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

        



void to_json(lm::Json& j, const StarSource& p) {
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


void from_json(const lm::Json& j, StarSource& p) {
    //do nothing;
}


class Renderer_Arepo_Preprocess final : public Renderer {

protected:
    std::string strategy_;
    Float mis_power_;
    Scene* scene_;
    AccelKnn* vrl_accel_;
    Volume_Arepo* volume_;
    Film* film_;
    int max_verts_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;
    Component::Ptr<lm::Light> reprlight_;
    Float impact_threshold_;

    std::vector<Float> star_coords_;
    std::vector<Float> star_ugriz_;
    Float model_scale_;


public:
    virtual void construct(const Json& prop) override {
        scene_ = json::comp_ref<Scene>(prop, "scene");
        volume_= json::comp_ref<Volume_Arepo>(prop, "volume");

        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");
        impact_threshold_ = json::value<Float>(prop, "impact_threshold",1.0);

        star_coords_ = lm::json::value<std::vector<Float>>(prop,"star_coords");
        star_ugriz_ = lm::json::value<std::vector<Float>>(prop,"star_ugriz");
        model_scale_ = lm::json::value<Float>(prop,"model_scale");
        
        
        //want to cast direct lighting on sensor, each sample treating one light
        


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
		//i am not actually rendering

        //reconstruct scheduler with new sample count
        Json info;
        info["num_samples"] = star_coords_.size() / 3;

        auto ontheflysched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::sample", make_loc("scheduler"), info);


        stats::clearGlobal<lm::stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<lm::stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<lm::stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<lm::stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<lm::stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<lm::stats::TotalTetraTests,int,long long>( );

        stats::clearGlobal<stats::LightsInTetra,stats::TetraIndex,std::deque<StarSource>>();


        timer::ScopedTimer st;



        //to sRGB
        auto m0 = Vec3(3.2404542,-0.9692660,0.0556434);     
        auto m1 = Vec3( -1.5371385, 1.8760108,-0.2040259);  
        auto m2 = Vec3(  -0.4985314, 0.0415560, 1.0572252);
        auto M = Mat3(m0,m1,m2);
      
        LM_INFO("star coords: {}, star photos {}",star_coords_.size()/3,star_ugriz_.size()/8);
        //createdLights_.reserve(star_coords_.size() / 3);


        //the scheduler gives one sample per light
        const auto processed = ontheflysched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            
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

            //LM_INFO("reproduced g mag {}",reproduced_mag);
            //LM_INFO("reproduced temp {}",temperature);
            //LM_INFO("incoming mag : {},{},{},{}, reproduced g mag: {}",
            //g_mag,r_mag,i_mag,z_mag,
        // reproduced_mag);

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
            

            /*lightprop["Le"] = rgb;

            lightprop["position"] = model_scale_ * Vec3(star_coords_[coord_i],star_coords_[coord_i+1],star_coords_[coord_i+2]);
            std::string name = "light" + std::to_string(ind);
            createdLights_.push_back(
                std::move(
                    comp::create<Light>(
                        "light::point", make_loc(name), lightprop)
                )
            );*/
            /*LM_INFO("star with jansky {},{},{},{},  with temp {}, rgb {},{},{}, at pos {},{},{}",
            g,r,i,z,
            temperature,
            rgb.r,
            rgb.g,
            rgb.b,
            star_coords_[coord_i],
            star_coords_[coord_i+1],star_coords_[coord_i+2]
            );*/

            //}
            



            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
            Vec2 raster_pos{};

            //store the sample id that this thread currently works on 
            stats::set<stats::CachedSampleId,int,long long>(0,sample_index);

            int currentBFSLayer = 0;
            bool addedAny = true;
            volume_->visitBFS(star.position, [&] (int tetraI,glm::tmat4x3<lm::Float> corners, int bfsLayer) -> bool {
                //first try: store where this star is located in. 
                //store information directly instead of indices. more cache efficient? but more storage.
                
                bool enteredNewLayer = currentBFSLayer != bfsLayer;
                bool haveAddedAnyInLastLayer = addedAny;
                
                if(enteredNewLayer) //reset
                    addedAny = false;
                currentBFSLayer = bfsLayer;

                Vec3 starpos = star.position;
                Float starintens = glm::max(star.intensity[0],glm::max(star.intensity[1],star.intensity[2]));
                
                glm::tmat4x3<lm::Float> pVs;
                connectP(corners,starpos,pVs);
                lm::Vec4 dets;
                computeDeterminants(pVs,dets);
                lm::Float mainDeterminant = det3x3(corners[0] - corners[3],corners[1] - corners[3],corners[2] - corners[3]);
                bool isInside = inside(dets, mainDeterminant);
                //smallest distance to tetra ? well something simila..?
                auto sth = glm::abs(mainDeterminant / (glm::abs(dets[0]) + glm::abs(dets[1]) + glm::abs(dets[2]) + glm::abs(dets[3])));
                
                //calculate tetra's minimum distance to star
                //auto d = glm::distance2(corners[0],starpos);
                //d = glm::min(d, glm::distance2(corners[1],starpos));
                //d = glm::min(d, glm::distance2(corners[2],starpos));
                //d = glm::min(d, glm::distance2(corners[3],starpos));

                //bool comp = (starintens / d) > impact_threshold_;

                bool comp = isInside || (sth < impact_threshold_);

                if(comp) {
                    stats::update<stats::LightsInTetra,stats::TetraIndex,std::deque<StarSource>>(
                        tetraI, 
                        [&](auto & vec) {
                            vec.push_back(star);
                        }
                    );
                    addedAny = true;
                }
             
                
                return !enteredNewLayer || (enteredNewLayer && haveAddedAnyInLastLayer); //continue if didnt stop adding

            });
    
            

        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {

            stats::clear<stats::LightsInTetra,stats::TetraIndex,std::deque<StarSource>>();

            stats::clear<lm::stats::SampleIdCacheHits,int,long long>( );
            stats::clear<lm::stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<lm::stats::UsedCachedTetra,int,long long>( );
            stats::clear<lm::stats::UsedNeighborTetra,int,long long>( );
            stats::clear<lm::stats::ResampleAccel,int,long long>( );
            stats::clear<lm::stats::TotalTetraTests,int,long long>( );

        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
            
        

            
            //merge the per tetrahedron vectors of star light sources
            stats::mergeToGlobal<stats::LightsInTetra,stats::TetraIndex,std::deque<StarSource>>( 
                [](auto & vector1, auto & vector2 ) { vector1.insert(vector1.begin(),vector2.begin(), vector2.end());return vector1;}
            );

            stats::clear<stats::LightsInTetra,stats::TetraIndex,std::deque<StarSource>>();

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
