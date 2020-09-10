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
#include "statstags.h"
#include "arepoloader.h"
#define VOLPT_IMAGE_SAMPLING 0

LM_NAMESPACE_BEGIN(LM_NAMESPACE)


class Renderer_VolPT_Base : public Renderer {
protected:
    std::string strategy_;
    Scene* scene_;
    Volume_Arepo* volume_;
    Film* film_;
    int max_verts_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;


public:
    LM_SERIALIZE_IMPL(ar) {
        ar(scene_, film_, max_verts_, rr_prob_, sched_);
    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        comp::visit(visit, scene_);
        comp::visit(visit, film_);
        comp::visit(visit, sched_);
    }

public:
    virtual void construct(const Json& prop) override {
        strategy_ = json::value<std::string>(prop, "strategy");
        scene_ = json::comp_ref<Scene>(prop, "scene");
        volume_= json::comp_ref<Volume_Arepo>(prop, "volume");
        film_ = json::comp_ref<Film>(prop, "output");
        max_verts_ = json::value<int>(prop, "max_verts");
        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");

        #if VOLPT_IMAGE_SAMPLING
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::" + sched_name, make_loc("scheduler"), prop);
        #else
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spp::" + sched_name, make_loc("scheduler"), prop);
        #endif

    }
};
class Renderer_VolPT_Arepo final : public Renderer_VolPT_Base {

public:
    virtual Json render() const override {
		scene_->require_renderable();

        stats::clearGlobal<lm::stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<lm::stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<lm::stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<lm::stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<lm::stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<lm::stats::TotalTetraTests,int,long long>( );


        film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;
        const auto processed = sched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            LM_KEEP_UNUSED(sample_index);

            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
             //LM_INFO("current sample id : {}", std::to_string(stats::get<stats::CurrentSampleId,long long>(0)));
            //
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_0,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_1,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_1_0,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_1_1,0.0);


            // ------------------------------------------------------------------------------------

            // Sample window
            Vec4 window(0,0,1,1);
            #if VOLPT_IMAGE_SAMPLING
            LM_UNUSED(pixel_index);
            #else
            {
                const int x = int(pixel_index % size.w);
                const int y = int(pixel_index / size.w);
                const auto dx = 1_f / size.w;
                const auto dy = 1_f / size.h;
                window = { dx * x, dy * y, dx, dy };
            }
            #endif

            // ------------------------------------------------------------------------------------

            // Sample initial vertex
            const auto sE = path::sample_position(rng, scene_, TransDir::EL);
            const auto sE_comp = path::sample_component(rng, scene_, sE->sp, {});
            auto sp = sE->sp;
            int comp = sE_comp.comp;
            auto throughput = sE->weight * sE_comp.weight;

            // ------------------------------------------------------------------------------------

            // Perform random walk
            Vec3 wi{};
            Vec2 raster_pos{};
            for (int num_verts = 1; num_verts < max_verts_; num_verts++) {


                //store the sample id that this thread currently works on 
                stats::set<stats::CachedSampleId,int,long long>(0,max_verts_ * spp_ * pixel_index + max_verts_ * sample_index + num_verts);
           

                // Sample a NEE edge
                #if VOLPT_IMAGE_SAMPLING
                const auto samplable_by_nee = !path::is_specular_component(scene_, sp, comp);
                #else
                const auto samplable_by_nee = num_verts > 1 && !path::is_specular_component(scene_, sp, comp);
                #endif
                if (samplable_by_nee) [&] {
                    // Sample a light
                    const auto sL = path::sample_direct(rng, scene_, sp, TransDir::LE);
                    if (!sL) {
                        return;
                    }

                    // Recompute raster position for the primary edge
                    Vec2 rp = raster_pos;
                    if (num_verts == 1) {
                        const auto rp_ = path::raster_position(scene_, -sL->wo);
                        if (!rp_) { return; }
                        rp = *rp_;
                    }

                    // Transmittance
                    const auto Tr = path::eval_transmittance(rng, scene_, sp, sL->sp);
                    if (math::is_zero(Tr)) {
                        return;
                    }

                    // Evaluate BSDF
                    const auto wo = -sL->wo;
                    const auto fs = path::eval_contrb_direction(scene_, sp, wi, wo, comp, TransDir::EL, true);
                    if (math::is_zero(fs)) {
                        return;
                    }

                    //int h = 0;
                    //lm::Float normFac = lm::stats::get<lm::stats::MaxTransmittance,int,lm::Float>(h);


                    // Evaluate and accumulate contribution
                    const auto C = throughput * Tr * fs * sL->weight;
                    film_->splat(rp, C);
                }();

                // --------------------------------------------------------------------------------

                // Sample direction
                const auto s = [&]() -> std::optional<path::DirectionSample> {
                    if (num_verts == 1) {
                        const auto [x, y, w, h] = window.data.data;
                        const auto ud = Vec2(x + w * rng.u(), y + h * rng.u());
                        return path::sample_direction({ ud, rng.next<Vec2>() }, scene_, sp, wi, comp, TransDir::EL);
                    }
                    else {
                        return path::sample_direction(rng, scene_, sp, wi, comp, TransDir::EL);
                    }
                }();
                if (!s) {
                    break;
                }

                // --------------------------------------------------------------------------------

                // Compute and cache raster position
                if (num_verts == 1) {
                    raster_pos = *path::raster_position(scene_, s->wo);
                }

                // --------------------------------------------------------------------------------

             
                // Sample next scene interaction
                
                //importance sample distance following light source, strategy 0

                
                auto lightDistanceSample = path::DistanceSample();
                lightDistanceSample.weight = Vec3(0);
                if(strategy_ != "delta" && strategy_ != "regular") {

                    auto mediumindex = scene_->medium_node();
                    auto & medium = scene_->node_at(mediumindex);
                    //medium.primitive.medium->
                    Ray ray = {
                    sp.geom.p,
                    s->wo};
                    
                    //sample next scene interaction following light sources situation
                    scene_->traverse_primitive_nodes(
                        [&] (const SceneNode& node, Mat4 global_transform) {
                            if(node.primitive.light != nullptr && !node.primitive.light->is_infinite()) {
                                auto & light = *node.primitive.light;
                                auto possample = rng.next<Light::PositionSampleU>();
                                auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;
                                //LM_INFO("light pos {},{},{}",lightPos.x,lightPos.y,lightPos.z);
                                auto shortestT = glm::dot((lightPos - ray.o), ray.d); 
                                auto d = ray.o + ray.d * shortestT - lightPos;
                                auto D = glm::length(d);
                                Float a = -10.0;//ray.o - lightPos;
                                Float b =  10.0;//ray.o + ray.d * 999999.0 - lightPos;
                                auto theta_a = glm::atan(a,D);
                                auto theta_b = glm::atan(b,D);
                                auto xi = rng.u();
                                auto equiT = D * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                auto pdf = D / ((theta_b - theta_a) * (D*D+equiT*equiT));
                                auto t = equiT + shortestT;
                                //store sample of strategy 0
                                int key = 0;
                                stats::set<stats::EquiangularStrategyDistanceSample,int,Float>(key,t);
                                //store pdf for sample of strategy 0
                                stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_0,pdf);

                                auto T = medium.primitive.medium->eval_transmittance(
                                 rng,  ray,0.0, t);
                                auto mu_s = volume_->eval_scalar(ray.o + ray.d * t,ray.d);
                                //assume mu_s with isotropic phase, i.e.,particle_dens *  1

                                auto albedo = Vec3( mu_s* T );



                                lightDistanceSample = path::DistanceSample {
                                    SceneInteraction::make_medium_interaction(
                                        scene_->medium_node(),
                                        PointGeometry::make_degenerated(ray.o + ray.d * t)
                                    ),
                                    albedo
                                };
                                
                            }
                    });
                }

                //importance sample distance following volume 
                auto sd = path::sample_distance(rng, scene_, sp, s->wo);
                
                

                
                //2 distance samples 0 and 1, 2 strategies 0 and 1,
                // we miss one pdf value: equiangular pdf 0 of regular distance sample 1
                if(strategy_ != "delta" && strategy_ != "regular") {

                    //auto mediumindex = scene_->medium_node();
                    //auto & medium = scene_->node_at(mediumindex);
                    //medium.primitive.medium->
                    Ray ray = {
                    sp.geom.p,
                    s->wo};
                    int key = 0;
                    auto regularT = stats::get<stats::RegularTrackingStrategyDistanceSample,int,Float>(key);
                                
                    //sample next scene interaction following light sources situation
                    scene_->traverse_primitive_nodes(
                        [&] (const SceneNode& node, Mat4 global_transform) {
                            if(node.primitive.light != nullptr && !node.primitive.light->is_infinite()) {
                                auto & light = *node.primitive.light;
                                auto possample = rng.next<Light::PositionSampleU>();
                                auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;
                                //LM_INFO("light pos {},{},{}",lightPos.x,lightPos.y,lightPos.z);
                                auto shortestT = glm::dot((lightPos - ray.o), ray.d); 
                                auto d = ray.o + ray.d * shortestT - lightPos;
                                auto D = glm::length(d);
                                Float a = -10.0;//ray.o - lightPos;
                                Float b =  10.0;//ray.o + ray.d * 999999.0 - lightPos;
                                auto theta_a = glm::atan(a,D);
                                auto theta_b = glm::atan(b,D);
                                //auto xi = rng.u();
                                auto equiT = regularT - shortestT;
                                auto pdf = D / ((theta_b - theta_a) * (D*D+equiT*equiT));
                                //store pdf for sample of strategy 0
                                stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_1,pdf);                                
                            }
                    });
                }
                auto _strat_smpl_key = stats::IJ::_0_0;
                auto equiStratEquiSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                
                _strat_smpl_key = stats::IJ::_0_1;
                auto equiStratRegularSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                            
                _strat_smpl_key = stats::IJ::_1_0;
                auto regularStratEquiSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                             
                _strat_smpl_key = stats::IJ::_1_1;
                auto regularStratRegularSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                          
                auto w_Equi = glm::pow(1.0 * equiStratEquiSmplPDF,1.0) / (glm::pow(1.0 * equiStratEquiSmplPDF,1.0) + glm::pow(1.0 * regularStratEquiSmplPDF,1.0));
                auto w_Regular = glm::pow(1.0 * regularStratRegularSmplPDF,1.0) / (glm::pow(1.0 * equiStratRegularSmplPDF,1.0) + glm::pow(1.0 * regularStratRegularSmplPDF,1.0));
                
                //first choice: leave sd as is

                //second choice: sample equiangular
                //sd = lightDistanceSample;


                //no medium interaction at all
                //todo evaluate this!!!
                if (!sd) {
                    break;
                }

                if(strategy_ == "experimental_mis") {

                    //third choice: balance heuristic
                    auto equiContr = lightDistanceSample.weight;
                    auto regularContr = sd->weight;
                    auto equiScatterPoint = lightDistanceSample.sp.geom.p;
                    auto regularScatterPoint = sd->sp.geom.p;
                    sd->sp.geom.p  =    w_Equi * lightDistanceSample.sp.geom.p / equiStratEquiSmplPDF
                                    +   w_Regular * sd->sp.geom.p  / regularStratRegularSmplPDF;
                                
                    sd->weight =  w_Equi * lightDistanceSample.weight / equiStratEquiSmplPDF
                                +   w_Regular * sd->weight  / regularStratRegularSmplPDF;
                } else if(strategy_ == "equiangular" ) {
                    sd = lightDistanceSample;
                    sd->weight = sd->weight / equiStratEquiSmplPDF;
                } else if(strategy_ == "regular") {
                    sd->weight = sd->weight / regularStratRegularSmplPDF;
                    //already done
                } else if(strategy_ == "onesample_mis") {
                    auto chooseStrategy =  rng.u();
                    if(chooseStrategy < 0.5) {
                        sd->weight = w_Regular * sd->weight / regularStratRegularSmplPDF / 0.5;
                    } else {
                        sd = lightDistanceSample;
                        sd->weight = w_Equi * sd->weight / equiStratEquiSmplPDF / 0.5;
                    }
                } else if(strategy_ == "delta") {
                    //everythin is as it should be
                }
                //LM_INFO("weights sum {}" ,Vec3(sd->weight  + lightDistanceSample.weight )[0]);
                //lightDistanceSample.weight * lightDistanceSample.sp.geom.p 
                //+ sd.weight * sd->sp.geom.p 
                //linear combine points
                //sd->sp.geom.p = lightDistanceSample.sp.geom.p * lightDistanceSample.weight * 0.5 +
                //                sd->sp.geom.p * sd->weight * 0.5;
                //draw another uniform value and take one
                //sd->sp.geom.p = rng.u() > 0.5 ? sd->sp.geom.p : lightDistanceSample.sp.geom.p;

                //always combine all weights together (right now with 0.5 weight)
                //sd->weight = 0.5 * lightDistanceSample.weight + 0.5 * sd->weight;

                

                //sd->weight = sd->weight / regularStratRegularSmplPDF;
                //sd->weight = sd->weight / equiStratEquiSmplPDF;
                
                
                // --------------------------------------------------------------------------------

                // Update throughput
                throughput *= s->weight * sd->weight;

                // --------------------------------------------------------------------------------

                // Accumulate contribution from emissive interaction
                if (!samplable_by_nee && scene_->is_light(sd->sp)) {
                    const auto spL = sd->sp.as_type(SceneInteraction::LightEndpoint);
                    const auto woL = -s->wo;
                    const auto Le = path::eval_contrb_direction(scene_, spL, {}, woL, comp, TransDir::LE, true);
                    const auto C = throughput * Le;
                    film_->splat(raster_pos, C);
                }
                
                // --------------------------------------------------------------------------------

                // Termination on a hit with environment
                if (sd->sp.geom.infinite) {
                    break;
                }

                // Russian roulette
                if (num_verts > 5) {
                    const auto q = glm::max(rr_prob_, 1_f - glm::compMax(throughput));
                    if (rng.u() < q) {
                        break;
                    }
                    throughput /= 1_f - q;
                }

                // --------------------------------------------------------------------------------

                // Sample component
                const auto s_comp = path::sample_component(rng, scene_, sd->sp, -s->wo);
                throughput *= s_comp.weight;

                // --------------------------------------------------------------------------------

                // Update information
                wi = -s->wo;
                sp = sd->sp;
                comp = s_comp.comp;
            }
        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::clear<lm::stats::SampleIdCacheHits,int,long long>( );
            stats::clear<lm::stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<lm::stats::UsedCachedTetra,int,long long>( );
            stats::clear<lm::stats::UsedNeighborTetra,int,long long>( );
            stats::clear<lm::stats::ResampleAccel,int,long long>( );
            stats::clear<lm::stats::TotalTetraTests,int,long long>( );

        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
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

        // Rescale film
        #if VOLPT_IMAGE_SAMPLING
        film_->rescale(Float(size.w* size.h) / processed);
        #else
        film_->rescale(1_f / processed);
        #endif

        return { {"processed", processed}, {"elapsed", st.now()} };
    }
};

LM_COMP_REG_IMPL(Renderer_VolPT_Arepo, "renderer::volpt_arepo");

LM_NAMESPACE_END(LM_NAMESPACE)
