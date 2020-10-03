/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/core.h>
#include <lm/renderer.h>
#include <lm/scene.h>
#include <lm/film.h>
#include <lm/scheduler.h>
#include <lm/path.h>
#include <lm/timer.h>
#include <lm/stats.h>
#include <lm/light.h>
#include <lm/accel.h>
#include "statstags.h"
#include "arepoloader.h"

LM_NAMESPACE_BEGIN(LM_NAMESPACE)



class Renderer_Arepo_VRL final : public Renderer {

protected:
    std::string strategy_;
    Float mis_power_;
    Scene* scene_;
    AccelKnn* vrl_accel_;
    Volume* volume_;
    Film* film_;
    int max_verts_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;
    Component::Ptr<lm::Light> reprlight_;


public:
    virtual void construct(const Json& prop) override {
        strategy_ = json::value<std::string>(prop, "strategy");
        mis_power_ = json::value<Float>(prop, "mis_power",2.0);
        scene_ = json::comp_ref<Scene>(prop, "scene");
        vrl_accel_ = json::comp_ref<AccelKnn>(prop, "vrl_accel");
        volume_= json::comp_ref<Volume>(prop, "volume");
        film_ = json::comp_ref<Film>(prop, "output");
        max_verts_ = json::value<int>(prop, "max_verts");
        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");


        
        
        //want to cast direct lighting on sensor, each sample treating one light
        
        auto copy = prop;
        copy["num_samples"] = scene_->num_lights();
        LM_INFO("num samples {}",scene_->num_lights());
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::" + sched_name, make_loc("scheduler"), copy);



    }

    LM_SERIALIZE_IMPL(ar) {
        ar(scene_, film_, max_verts_, rr_prob_, sched_);
    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        comp::visit(visit, scene_);
        comp::visit(visit, film_);
        comp::visit(visit, sched_);
        
    }


    virtual Json render() const override {
		scene_->require_renderable();

        //reconstruct scheduler with new sample count

        Json copy;
        copy["num_samples"] = scene_->num_lights();
        LM_INFO("num samples {}",scene_->num_lights());
        
        auto ontheflysched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::sample", make_loc("scheduler"), copy);


        stats::clearGlobal<lm::stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<lm::stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<lm::stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<lm::stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<lm::stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<lm::stats::TotalTetraTests,int,long long>( );

        stats::clearGlobal<stats::VRL,stats::TetraIndex,LightToCameraRaySegmentCDF>();
        stats::clearGlobal<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>();


        film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;


        //the scheduler gives one sample per light
        const auto processed = ontheflysched_->run([&](long long pixel_index, long long sample_index, int threadid) {

            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
            Vec2 raster_pos{};

            //store the sample id that this thread currently works on 
            stats::set<stats::CachedSampleId,int,long long>(0,sample_index);

            //LM_INFO("i am thread {} and work on sample {}", threadid, sample_index);
            if(sample_index >= scene_->num_lights())
                return;

            //my light
            auto lightPrimitiveIndex = scene_->light_primitive_index_at(sample_index);
            
            auto globalT = lightPrimitiveIndex.global_transform;
           // LM_INFO("t {}",globalT.M[0][0]);
            auto mylightSceneNode = scene_->node_at( lightPrimitiveIndex.index);
            //LM_INFO("l{}",mylightSceneNode.primitive.light != nullptr ? 1  :0 );
            assert(mylightSceneNode.primitive.light != nullptr);
            auto positionsample = mylightSceneNode.primitive.light->sample_position(rng.next<Light::PositionSampleU>(), globalT);
            //LM_INFO("p{}",positionsample.has_value());
            assert(positionsample);
            auto pos = positionsample->geom.p;
            auto sp = SceneInteraction::make_light_endpoint(lightPrimitiveIndex.index, positionsample->geom);

            auto cam_light_connection =
            path::sample_direct(rng,scene_, sp, TransDir::EL);

            assert(cam_light_connection);


            //now sample the contribution of the light, seen from camera endpoint (then need to evaluate transmittance, see further down)
            auto lightSample =  
                mylightSceneNode.primitive.light->sample_direct(rng.next<Light::RaySampleU>(),cam_light_connection->sp.geom, globalT);

            assert(lightSample);


            //hash map to save nearest lights
            //auto nearestLights = std::unordered_map<int,Float>();
            //nearestLights.reserve(scene_->num_lights() * 0.1);

            //save nearest lights along current ray 

            LightToCameraRaySegmentCDF segment;
            segment.weight = lightSample->weight;
            segment.cdfSoFar = 0.0;
            segment.localcdf = 0.0;
            segment.t = 0.0;
            segment.a = 0.0;
            segment.b = 0.0;
            
            std::function<void(Vec3,RaySegmentCDF const &, int)> raysegmentVisitor = [&] (lm::Vec3 boundarypos,lm::RaySegmentCDF const & tetrasegment, int tetraI) -> void {
                //add an entry for the current tetrahedron
                segment.localcdf = tetrasegment.localcdf;
                segment.t = tetrasegment.t;
                segment.a = tetrasegment.a;
                segment.b = tetrasegment.b;
                segment.p = boundarypos;
                segment.d = -cam_light_connection->wo;
                stats::enqueue<stats::VRL,stats::TetraIndex,LightToCameraRaySegmentCDF>(std::move(tetraI),std::move(segment));
                segment.tSoFar += tetrasegment.t;
                segment.cdfSoFar += tetrasegment.localcdf;
            };
            stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &,int)>>(0,raysegmentVisitor);
            //std::function<void(Vec3,RaySegmentCDF const &)> donothing = [](auto,auto) {}; 
            //stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &)>>(0,&donothing);

            //ask nearest light to current vertex
            //scene_->sample_light_selection_from_pos(0.0,sp.geom.p);

            //lm::Scene::LightPrimitiveIndex lightprim; 

            // Sample next scene interaction

            //importance sample distance following light source, strategy 0

            //importance sample distance following volume 
            std::optional<path::DistanceSample> sd = path::sample_distance(rng, scene_, sp, -cam_light_connection->wo);

            //this is very important, to unregister the function in subsequent render passes(!)
            stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &,int)>>(0,{});

            //splat result to pixel (test)
            const auto rp_ = path::raster_position(scene_, cam_light_connection->wo);
            auto Tr = glm::exp(- segment.cdfSoFar);
            if(rp_) { //it is part of the sensor?
                auto tosplat = Tr * segment.weight;
                //LM_INFO("splat  {},{},{} on sensor {},{}",tosplat[0],tosplat[1],tosplat[2],
                //rp_->x,rp_->y);
                film_->splat(*rp_,Tr * segment.weight);
            }
            

        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {

            stats::flushQueue<stats::VRL,stats::TetraIndex,LightToCameraRaySegmentCDF>([](auto & k, auto & v) {});
            
            //to be enqueued in the render loop
            stats::clear<stats::VRL,stats::TetraIndex,LightToCameraRaySegmentCDF>();
            //to be filled during flush queue after render loop
            stats::clear<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>();


            stats::clear<lm::stats::SampleIdCacheHits,int,long long>( );
            stats::clear<lm::stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<lm::stats::UsedCachedTetra,int,long long>( );
            stats::clear<lm::stats::UsedNeighborTetra,int,long long>( );
            stats::clear<lm::stats::ResampleAccel,int,long long>( );
            stats::clear<lm::stats::TotalTetraTests,int,long long>( );

        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
            
            //flush the queue to thread local per-tetrahedron vectors of raysegments
            stats::flushQueue<stats::VRL,stats::TetraIndex,LightToCameraRaySegmentCDF>(
                [](auto & k, auto & v) {
                    stats::update<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>(
//                        k,
                        0, //dont use tetra index, just add to 0
                        [&](auto & vec) {
                            vec.push_back(v);
                        }
                    );
                }
            );
            //merge the per tetrahedron vectors of ray segments to the global
            stats::mergeToGlobal<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>( 
                [](auto & vector1, auto & vector2 ) { vector1.insert(vector1.begin(),vector2.begin(), vector2.end());return vector1;}
            );


            stats::mergeToGlobal<lm::stats::SampleIdCacheHits,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::SampleIdCacheMisses,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::UsedCachedTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::UsedNeighborTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::ResampleAccel,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::TotalTetraTests,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );

        }
        );

        auto smplhits = stats::getGlobal<lm::stats::SampleIdCacheHits,int,long long>(0 );
        auto smplmisses = stats::getGlobal<lm::stats::SampleIdCacheMisses,int,long long>(0 );
        auto tetrahits =  stats::getGlobal<lm::stats::UsedCachedTetra,int,long long>( 0);
        auto tetraneighborhits = stats::getGlobal<lm::stats::UsedNeighborTetra,int,long long>(0 );
        auto accelsmpls = stats::getGlobal<lm::stats::ResampleAccel,int,long long>( 0);
        auto totaltetratests = stats::getGlobal<lm::stats::TotalTetraTests,int,long long>(0 );


        

     
        auto & vrls = stats::getGlobalRef<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>( )
            [0];


        /*stats::getGlobalRef<stats::KNNLineComp,int,std::function<void(int,Vec3&,Vec3&)>>()[0]
        = [&] (int vrlI,Vec3&vrlstart,Vec3&vrlend) {
            auto & vrl = stats::getGlobalRef<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>( )[0]
            [vrlI];
            vrlstart = vrl.p;
            vrlend = vrl.p + vrl.d * vrl.t;
        };*/
        
        int numvrls = vrls.size();
        Json prop;
        Mat4 transform;

        
        //vrl_temp_scene->add_child(vrl_temp_scene->root_node(), t);
        LM_INFO("vrl count: {}",numvrls);
       
        int currentIndex = 0;
        //assume everything is stored in one vector
#ifdef USE_KNN

        if(numvrls > 0) {
            std::function<bool(Mat4&,int&)> nextObject = [&](Mat4 & out_global_transform,int & out_someindex) {
                auto & vrl = vrls[currentIndex];
                transform = glm::mat4(1.0);
                transform[3] = glm::vec4(vrl.p + vrl.d * vrl.t * 0.5,1);
                out_someindex = currentIndex;
                currentIndex++;
                //LM_INFO("vrl {}",currentIndex);
                if(currentIndex <= numvrls)
                    return true;
                return false;
            };

            vrl_accel_->build(numvrls,nextObject);


        }
#endif
        


        LM_INFO("sample hits: {}, misses : {}, tetra hits {}, tetra neighbor hits {}, accel smpls {} . total tetra probes {}", 
         smplhits,
         smplmisses,
         tetrahits,
         tetraneighborhits,
         accelsmpls,
         totaltetratests
         );

        // Rescale film
        #if VOLPT_IMAGE_SAMPLING
        film_->rescale(Float(size.w* size.h) / processed);
        #else
        film_->rescale(1_f / processed);
        #endif

        return { {"processed", processed}, {"elapsed", st.now()} };
    }
};

LM_COMP_REG_IMPL(Renderer_Arepo_VRL, "renderer::volpt_arepo_vrl");

LM_NAMESPACE_END(LM_NAMESPACE)
