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
    Scene* scene_;
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
        scene_ = json::comp_ref<Scene>(prop, "scene");
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
            
             //store the sample id that this thread currently works on 
            stats::set<stats::CachedSampleId,int,long long>(0,spp_ * pixel_index + sample_index);
           //LM_INFO("current sample id : {}", std::to_string(stats::get<stats::CurrentSampleId,long long>(0)));

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
                const auto sd = path::sample_distance(rng, scene_, sp, s->wo);
                if (!sd) {
                    break;
                }

                //sample next scene interaction following light sources situation
                scene_->traverse_primitive_nodes(
                    [&] (const SceneNode& node, Mat4 global_transform) {
                        if(node.primitive.light != nullptr) {
                            auto & light = *node.primitive.light;
                            LM_INFO("found a light {}",node.index);
                        }
                });


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
