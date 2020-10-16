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

        


void to_json(Json& j, const StarSource& p) {
    j = {
        {"intensity" , {
            {"r",p.intensity.x},
            {"g",p.intensity.y},
            {"b",p.intensity.z}
            }
        },
        {"position" , {
            {"x",p.position.x},
            {"y",p.position.y},
            {"z",p.position.z}
            }
        },
        {"index" , p.index}
        
    };
}


void from_json(const Json& j, StarSource& p) {
    //do nothing;
}


class Renderer_Arepo_Preprocess final : public Renderer {

protected:
    std::string strategy_;
    Float mis_power_;
    Scene* scene_;
    AccelKnn* light_accel_;
    Volume_Arepo* volume_;
    Film* film_;
    int max_verts_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;
    Component::Ptr<Light> reprlight_;
    Float impact_threshold_;

    std::vector<Float> star_coords_;
    std::vector<Float> star_ugriz_;
    Float model_scale_;

    std::string mode_;




public:
    virtual void construct(const Json& prop) override {
        scene_ = json::comp_ref<Scene>(prop, "scene");
        volume_= json::comp_ref<Volume_Arepo>(prop, "volume");

        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");
        impact_threshold_ = json::value<Float>(prop, "impact_threshold",1.0);

        star_coords_ = json::value<std::vector<Float>>(prop,"star_coords");
        star_ugriz_ = json::value<std::vector<Float>>(prop,"star_ugriz");
        model_scale_ = json::value<Float>(prop,"model_scale");

        mode_ = json::value<std::string>(prop, "mode"); 
        
        light_accel_ = json::comp_ref<AccelKnn>(prop, "light_accel");

        
        Json info;
        info["num_samples"] = star_coords_.size() / 3;

        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::sample", make_loc("scheduler"), info);


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

        
        stats::clearGlobal<stats::CachedSampleId,int,long long>( );

        stats::clearGlobal<stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<stats::TotalTetraTests,int,long long>( );

        //stats::clearGlobal<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>();
        stats::clearGlobal<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>();
        stats::clearGlobal<stats::LightSet,int,std::vector<StarSource>>();

        stats::clearGlobal<stats::TetsPerLight,int,Float>();



        timer::ScopedTimer st;



        //to sRGB
        auto m0 = Vec3(3.2404542,-0.9692660,0.0556434);     
        auto m1 = Vec3( -1.5371385, 1.8760108,-0.2040259);  
        auto m2 = Vec3(  -0.4985314, 0.0415560, 1.0572252);
        auto M = Mat3(m0,m1,m2);
      
        LM_INFO("star coords: {}, star photos {}",star_coords_.size()/3,star_ugriz_.size()/8);
        //createdLights_.reserve(star_coords_.size() / 3);


        //the scheduler gives one sample per light
        const auto processed = sched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            if(sample_index > star_coords_.size()/3) {
                LM_INFO("YOU THERE! STOOP!");
                return;
            }


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
            

            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
            Vec2 raster_pos{};

            if(mode_ != "all_lights") {

                //store the sample id that this thread currently works on 
                stats::set<stats::CachedSampleId,int,long long>(0,sample_index);

                int currentBFSLayer = 0;
                bool addedAny = true;
                int numLs = 0;
                volume_->visitBFS(star.position, [&] (int tetraI,glm::tmat4x3<Float> corners, int bfsLayer) -> bool {
                    //first try: store where this star is located in. 
                    //store information directly instead of indices. more cache efficient? but more storage.
                    //no store indices.
                    bool enteredNewLayer = currentBFSLayer != bfsLayer;
                    bool haveAddedAnyInLastLayer = addedAny;
                    
                    if(enteredNewLayer) //reset
                        addedAny = false;
                    currentBFSLayer = bfsLayer;

                    Vec3 starpos = star.position;
                    Float starintens = glm::max(star.intensity[0],glm::max(star.intensity[1],star.intensity[2]));
                    
                    bool isInside = false;
                    Float mainDeterminant = 0.0;
                    glm::tmat4x3<Float> pVs;
                    {
                        connectP(corners,starpos,pVs);
                        Vec4 dets;
                        computeDeterminants(pVs,dets);
                        mainDeterminant = det3x3(corners[0] - corners[3],corners[1] - corners[3],corners[2] - corners[3]);
                        isInside = inside(dets, mainDeterminant);
                    }

                    //take farthest tetra point as sphere center
                    //and longest edge as radius, test if inside radius

                    int farthestI;
                    Float dist = 0.0;
                    //compute farthest corner 
                    for(int i = 0; i < 4; i ++) {
                        auto d = glm::length2(pVs[i]);
                        if(d > dist) {
                            farthestI = i;
                            dist = d;
                        }
                    }
                    //compute longest edge in tetra
                    dist = 0.0;
                    for(int i = 1; i < 4; i++)
                        dist = glm::max( dist, 
                            glm::distance(corners[(farthestI + i) % 4], corners[farthestI]));
                    //compute distance between sphere and star
                    dist = glm::max(0.0,glm::distance(starpos, corners[farthestI]) - dist);


                    //smallest distance to tetra ? well, the above is something similar..?
                    
                    bool comp = isInside || (starintens / dist) > impact_threshold_;//BIAS


                    if(comp) {
                        //stats::update<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>(
                        stats::update<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>(
                            tetraI, 
                            [&](auto & vec) {
                            //    vec.push_back(star);
                                vec.push_back(star.index);
                            }
                        );
                        numLs++;
                        addedAny = true;
                    }
                
                    bool continu = !enteredNewLayer || (enteredNewLayer && haveAddedAnyInLastLayer);
                    if(!continu) {
                        stats::set<stats::TetsPerLight,int,Float>(ind,numLs);
                        //LM_INFO("gonna stop at bfs layer {}, add {} myself",bfsLayer,numLs);
                    }
                    return continu; //continue if didnt stop adding

                });
            }// else if(mode_ == "all_lights") {
                stats::update<stats::LightSet,int,std::vector<StarSource>>(
                    0, 
                    [&](auto & vec) {
                        vec.push_back(star);
                    }
                );

            //}
    
            

        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::clear<stats::CachedSampleId,int,long long>();

            //stats::clear<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>();
            stats::clear<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>();
            stats::clear<stats::LightSet,int,std::vector<StarSource>>();

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
            //stats::mergeToGlobal<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>( 
            stats::mergeToGlobal<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>( 
                [] (auto & vector1, auto & vector2 )   { vector1.insert(vector1.end(),vector2.begin(), vector2.end());return vector1;}
            );

            //if strategy is all lights, merge like this
            stats::mergeToGlobal<stats::LightSet,int,std::vector<StarSource>>( 
                [] (auto & vector1, auto & vector2 )   { vector1.insert(vector1.end(),vector2.begin(), vector2.end());return vector1;}
            );

            //stats::clear<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>();
            stats::clear<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>();

            stats::clear<stats::LightSet,int,std::vector<StarSource>>();


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

        LM_INFO("average number of tetrahedra a light is associated with: {}", avg);

        auto lights = stats::getGlobalRef<stats::LightSet,int,std::vector<StarSource>>()[0];
        LM_INFO("in case of strategy all lights: total number of lights: {}", lights.size());
#ifdef USE_KNN_EMBREE
        if(lights.size() > 0) {
            int currentIndex = 0;
            auto nextObject = 
            [&] (Vec3 & out_global_transform,int & out_someindex,Float & out_someradius) -> bool {
                auto & light = lights[currentIndex];
                out_global_transform = light.position;
                out_someindex = currentIndex;
                Float starintens = glm::max(light.intensity[0],glm::max(light.intensity[1],light.intensity[2]));

                out_someradius = starintens /impact_threshold_;
                currentIndex++;
                if(currentIndex <= lights.size())
                    return true;
                return false;
            };
            light_accel_->build(lights.size(),nextObject);
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

        

        return { {"processed", processed}, {"elapsed", st.now()} };
    }
};

LM_COMP_REG_IMPL(Renderer_Arepo_Preprocess, "renderer::volpt_arepo_preprocess");

LM_NAMESPACE_END(LM_NAMESPACE)
