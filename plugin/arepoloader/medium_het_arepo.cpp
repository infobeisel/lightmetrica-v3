/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/core.h>
#include <lm/medium.h>
#include <lm/surface.h>
#include <lm/phase.h>
#include <lm/volume.h>

#include "arepoloader.h"

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

/*
\rst
.. function:: medium::heterogeneous

    Heterogeneous medium.

    :param str density: Locator to ``volume`` asset representing density of the medium.
    :param str albedo: Locator to ``volume`` asset representing albedo of the medium.
    :param str phase: Locator to ``phase`` asset.
\endrst
*/
class Medium_Het_Arepo final : public Medium {
private:
    const Volume_Arepo* volume_density_;  // Density volume. density := \mu_t = \mu_a + \mu_s
    const Volume* volme_albedo_;	// Albedo volume. albedo := \mu_s / \mu_t
    const Phase* phase_;            // Underlying phase function.
    bool deltatracking = false;
    bool extrapolation = false;

public:
    LM_SERIALIZE_IMPL(ar) {
        ar(volume_density_, volme_albedo_, phase_);
    }

public:
    virtual void construct(const Json& prop) override {
        volume_density_ = json::comp_ref<Volume_Arepo>(prop, "volume_density");
        volme_albedo_ = json::comp_ref<Volume>(prop, "volume_albedo");
        phase_ = json::comp_ref<Phase>(prop, "phase");
        deltatracking = lm::json::value<int>(prop,"deltatracking",0) == 1 ? true : false;
        extrapolation = lm::json::value<int>(prop,"extrapolation",0) == 1 ? true : false;
        
    }

    virtual std::optional<DistanceSample> sample_distance(Rng& rng, Ray ray, Float tmin, Float tmax) const override {
        // Compute overlapping range between the ray and the volume
        if (!volume_density_->bound().isect_range(ray, tmin, tmax)) {
            // No intersection with the volume, use surface interaction
            return {};
        }
        
        if(deltatracking) {
            Float validMaxT = -1.0;
            Float t = tmin;
            Float t_it = 0.0;
            Float t_acc = 0.0;
            Float inv_max_density;
            Float aCoef, bCoef; // scalar = aCoef * t + bCoef 
            while (true) {
                if(t_acc >= validMaxT) { // update max scalar
                    inv_max_density = 1_f / volume_density_->max_scalar({ray.o + ray.d * t,ray.d},validMaxT, aCoef, bCoef);
                    t_acc = 0.0;
                }
                t_it = -glm::log(1_f-rng.u()) * inv_max_density;
                // Sample a distance from the 'homogenized' volume
                t += t_it;
                t_acc += t_it;
                if (t >= tmax) {
                    // Hit with boundary, use surface interaction
                    return {};
                }

                // Density at the sampled point
                const auto p = ray.o + ray.d*t;
                //const auto density = volume_density_->eval_scalar(p, ray.d);
                //const auto density = bCoef + volume_density_->eval_scalar(p, ray.d);
                const auto density = extrapolation ? 
                    bCoef + t_acc * aCoef:
                    volume_density_->eval_scalar(p, ray.d);

                // Determine scattering collision or null collision
                // Continue tracking if null collusion is seleced
                if (density * inv_max_density > rng.u()) {
                    // Scattering collision
                    const auto albedo = volme_albedo_->eval_color(p);
                    return DistanceSample{
                        ray.o + ray.d*t,
                        albedo,     // T_{\bar{\mu}}(t) / p_{\bar{\mu}}(t) * \mu_s(t)
                                    // = 1/\mu_t(t) * \mu_s(t) = albedo(t)
                        true
                    };
                }
            }
        } else {
            
            Float out_maxTransmittance = 0.0;
            Float t = volume_density_->sample_distance(ray,tmin,tmax, rng, out_maxTransmittance);
            //LM_INFO("result {} ----------------------------------------------", t);
            
            if(t >= tmax) {
                // Hit with boundary, use surface interaction
                return {};
            }

            t = glm::min(t, tmax);

            auto p = ray.o + ray.d*t;
            const auto albedo = volme_albedo_->eval_color(p);
            return DistanceSample{
                p,
                albedo,     // T_{\bar{\mu}}(t) / p_{\bar{\mu}}(t) * \mu_s(t)
                            // = 1/\mu_t(t) * \mu_s(t) = albedo(t)
                true
            };
        }
        LM_UNREACHABLE_RETURN();
    }
    
    virtual Vec3 eval_transmittance(Rng& rng, Ray ray, Float tmin, Float tmax) const override {
        // Compute overlapping range
        if (!volume_density_->bound().isect_range(ray, tmin, tmax)) {
            // No intersection with the volume, no attenuation
            return Vec3(1_f);
        }

        if(!deltatracking) {
            //LM_INFO("eval transmittance from {} to {} ---------------------------------------", tmin, tmax);

            auto ret = volume_density_->eval_transmittance(ray, tmin, tmax);


            //LM_INFO("result {} ----------------------------------------------", ret);

            return Vec3(ret);
        } else {
            // Perform ratio tracking [Novak et al. 2014]
            Float Tr = 1_f;
            Float validMaxT = -1.0;
            Float t = tmin;
            Float t_it = 0.0;
            Float t_acc = 0.0;
            Float inv_max_density;
            Float aCoef, bCoef; // scalar = aCoef * t + bCoef 

            while (true) {
                if(t_acc >= validMaxT) {// update max scalar
                    inv_max_density = 1_f / volume_density_->max_scalar({ray.o + ray.d * t,ray.d},validMaxT,aCoef, bCoef);
                    t_acc = 0.0;
                }
                t_it = -glm::log(1_f - rng.u()) * inv_max_density;
                t += t_it;
                t_acc += t_it;
                if (t >= tmax) {
                    break;
                }
                const auto p = ray.o + ray.d*t;
                //const auto density = volume_density_->eval_scalar(p, ray.d);
                const auto density = extrapolation ? 
                    bCoef + t_acc * aCoef:
                    volume_density_->eval_scalar(p, ray.d);

                Tr *= 1_f - density * inv_max_density;
            }

            return Vec3(Tr);
        }
    }

    virtual bool is_emitter() const override {
        return false;
    }

    virtual const Phase* phase() const override {
        return phase_;
    }
};

LM_COMP_REG_IMPL(Medium_Het_Arepo, "medium::het_arepo");

LM_NAMESPACE_END(LM_NAMESPACE)
