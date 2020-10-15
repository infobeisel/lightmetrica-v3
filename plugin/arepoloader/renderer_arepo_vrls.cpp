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
using namespace ArepoLoaderInternals;

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
    Float impact_threshold_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;
    Component::Ptr<Light> reprlight_;
    bool save_vrls_;

    std::vector<Float> star_coords_;
    std::vector<Float> star_ugriz_;
    Float model_scale_;
    Vec3 camera_pos_;

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
        save_vrls_ = json::value<bool>(prop, "save_vrls",true);
        impact_threshold_ = json::value<Float>(prop, "impact_threshold");

        star_coords_ = json::value<std::vector<Float>>(prop,"star_coords");
        star_ugriz_ = json::value<std::vector<Float>>(prop,"star_ugriz");
        model_scale_ = json::value<Float>(prop,"model_scale");
        camera_pos_ = json::value<Vec3>(prop,"camera_pos");

        
        
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

        stats::clearGlobal<stats::CachedSampleId,int,long long>( );

        stats::clearGlobal<stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<stats::TotalTetraTests,int,long long>( );

        stats::clearGlobal<stats::VRL,stats::TetraIndex,std::deque<LightToCameraRaySegmentCDF>>();
        stats::clearGlobal<stats::TetsPerLight,int,Float>();


        film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;

        //to sRGB
        auto m0 = Vec3(3.2404542,-0.9692660,0.0556434);     
        auto m1 = Vec3( -1.5371385, 1.8760108,-0.2040259);  
        auto m2 = Vec3(  -0.4985314, 0.0415560, 1.0572252);
        auto M = Mat3(m0,m1,m2);
      


        //the scheduler gives one sample per light
        const auto processed = ontheflysched_->run([&](long long pixel_index, long long sample_index, int threadid) {

            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
            Vec2 raster_pos{};

            //store the sample id that this thread currently works on 
            stats::set<stats::CachedSampleId,int,long long>(0,sample_index);

            //parse light
            Json lightprop;
            auto ind = sample_index;
            //for(int ind = 0; ind < star_coords_.size() / 3; ind++) {
            auto ugriz_i = ind * 8;
            auto coord_i = ind * 3;
            auto u_mag = star_ugriz_[ugriz_i + 0];
            auto g_mag = star_ugriz_[ugriz_i + 4];
            auto r_mag = star_ugriz_[ugriz_i + 5];
            auto i_mag = star_ugriz_[ugriz_i + 6];
            auto z_mag = star_ugriz_[ugriz_i + 7];

            auto u_flux = u_mag;
            auto g_flux = g_mag;
            auto r_flux = r_mag;
            auto i_flux = i_mag;
            auto z_flux = z_mag;

            //LM_INFO(" g mag {}",g_mag);
            mag_to_flux(u_flux,g_flux,r_flux,i_flux,z_flux);

            auto sum = g_flux +  r_flux +  i_flux +  z_flux;
            auto g_ratio = g_flux / sum;
            auto temperature = poly10th(g_ratio);

            auto g_response_flux = sun_solid_angle*applySdssGBand(temperature);
            auto scaling = g_flux / g_response_flux;

            auto reproduced_mag = scaling * g_response_flux;
            flux_to_mag_g(reproduced_mag);


            Float x_,y_,z_;
            tempToXYZ(temperature,x_,y_,z_);
            //LM_INFO("star xyz {},{},{}",x_,y_,z_);
            auto rgb = M * Vec3(x_,y_,z_); //if it was the sun 
            rgb *= scaling; //account for star/sun relative size
            rgb *=  6.09e+18;//convert from area light source to point light source: sun's area in m^2 
            //LM_INFO("star rgb {},{},{}",rgb[0],rgb[1],rgb[2]);
            rgb = glm::max(Vec3(0.0),rgb); //CLAMP: for low temperatures there is negative colour from our model. 
            StarSource star;
            star.position = model_scale_ * Vec3(star_coords_[coord_i],star_coords_[coord_i+1],star_coords_[coord_i+2]);
            star.intensity = rgb;
            star.index = ind;

            auto cam_light_connection = glm::normalize(star.position - camera_pos_);
 

            auto cdfSoFar = 0.0;
            auto tSoFar = 0.0;
            

            std::function<void(Vec3,RaySegmentCDF const &, int)> raysegmentVisitor = [&] (Vec3 boundarypos,RaySegmentCDF const & tetrasegment, int tetraI) -> void {
                //add an entry for the current tetrahedron
                LightToCameraRaySegmentCDF segment;
                segment.weight = star.intensity;
                segment.localcdf = tetrasegment.localcdf;
                segment.t = tetrasegment.t;
                segment.a = tetrasegment.a;
                segment.b = tetrasegment.b;
                segment.p = boundarypos;
                segment.d = -cam_light_connection;
                segment.tSoFar = tSoFar;
                segment.cdfSoFar = cdfSoFar;

                Float intensity = glm::max(segment.weight[0],glm::max(segment.weight[1],segment.weight[2]));

                bool comp = intensity * glm::exp(-segment.cdfSoFar) > impact_threshold_;//BIAS
                if(save_vrls_ && comp) {
                    
                    stats::update<stats::VRL,stats::TetraIndex,std::deque<LightToCameraRaySegmentCDF>>(
                        tetraI, 
                        [&](auto & vec) {
                            vec.push_back(segment);
                        }
                    );
                    stats::add<stats::TetsPerLight,int,Float>(tetraI,1);

                }
                cdfSoFar += tetrasegment.localcdf;
                tSoFar += tetrasegment.t;
            };
            stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &,int)>>(0,raysegmentVisitor);
          
            //importance sample distance following volume 
            const auto* medium = scene_->node_at(scene_->medium_node()).primitive.medium;
            const auto ds = medium->sample_distance(rng, {star.position, -cam_light_connection }, 0_f, Inf);


            //this is very important, to unregister the function in subsequent render passes(!)
            stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &,int)>>(0,{});

            //splat result to pixel (test)
            const auto rp_ = path::raster_position(scene_, cam_light_connection);
            auto Tr = glm::exp(- cdfSoFar);
            if(rp_) { //it is part of the sensor?
                auto tosplat = Tr * star.intensity;
                //LM_INFO("splat  {},{},{} on sensor {},{}",tosplat[0],tosplat[1],tosplat[2],
                //rp_->x,rp_->y);
                //LM_INFO("{},{},{}",sp.geom.p.x,sp.geom.p.y,sp.geom.p.z);

                film_->splat(*rp_,tosplat);
            }
            

        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {

            stats::clear<stats::VRL,stats::TetraIndex,std::deque<LightToCameraRaySegmentCDF>>();

            stats::clear<stats::CachedSampleId,int,long long>();
            stats::clear<stats::TetsPerLight,int,Float>();

            stats::clear<stats::SampleIdCacheHits,int,long long>( );
            stats::clear<stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<stats::UsedCachedTetra,int,long long>( );
            stats::clear<stats::UsedNeighborTetra,int,long long>( );
            stats::clear<stats::ResampleAccel,int,long long>( );
            stats::clear<stats::TotalTetraTests,int,long long>( );

        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
            
            
            //merge the per tetrahedron vectors of star light sources
            stats::mergeToGlobal<stats::VRL,stats::TetraIndex,std::deque<LightToCameraRaySegmentCDF>>( 
                [](auto & vector1, auto & vector2 ) { vector1.insert(vector1.begin(),vector2.begin(), vector2.end());return vector1;}
            );

            stats::clear<stats::VRL,stats::TetraIndex,std::deque<LightToCameraRaySegmentCDF>>();


            stats::mergeToGlobal<stats::TetsPerLight,int,Float>( 
                [](Float & v0,Float & v1 ) { return v0 + v1;} );

            stats::mergeToGlobal<stats::SampleIdCacheHits,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::SampleIdCacheMisses,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::UsedCachedTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::UsedNeighborTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::ResampleAccel,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::TotalTetraTests,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );

        }
        );

        auto smplhits = stats::getGlobal<stats::SampleIdCacheHits,int,long long>(0 );
        auto smplmisses = stats::getGlobal<stats::SampleIdCacheMisses,int,long long>(0 );
        auto tetrahits =  stats::getGlobal<stats::UsedCachedTetra,int,long long>( 0);
        auto tetraneighborhits = stats::getGlobal<stats::UsedNeighborTetra,int,long long>(0 );
        auto accelsmpls = stats::getGlobal<stats::ResampleAccel,int,long long>( 0);
        auto totaltetratests = stats::getGlobal<stats::TotalTetraTests,int,long long>(0 );
        
        auto & tetsPerLight = stats::getGlobalRef<stats::TetsPerLight,int,Float>();
        Float numlights = tetsPerLight.size();
        Float avg = 0.0;
        for(auto p : tetsPerLight) {
            avg += p.second / numlights;
        }
        LM_INFO("average number vrl per tetra: {}", avg);

#ifdef USE_KNN_EMBREE

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
        //film_->rescale(1_f / processed);
        #endif

        return { {"processed", processed}, {"elapsed", st.now()} };
    }
};

LM_COMP_REG_IMPL(Renderer_Arepo_VRL, "renderer::volpt_arepo_vrl");

LM_NAMESPACE_END(LM_NAMESPACE)
