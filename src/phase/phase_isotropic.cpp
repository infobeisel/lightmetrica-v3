/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/core.h>
#include <lm/phase.h>
#include <lm/surface.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

/*
\rst
.. function:: phase::isotropic

    Isotropic phase function.
\endrst
*/
class Phase_Isotropic final : public Phase {
public:
    virtual std::optional<DirectionSample> sample_direction(const DirectionSampleU& us, const PointGeometry& geom, Vec3) const override {
        LM_UNUSED(geom);
        assert(geom.degenerated);
        return DirectionSample{
            math::sample_uniform_sphere(us.ud),
            Vec3(1_f)
        };
    }

    virtual Float pdf_direction(const PointGeometry& geom, Vec3, Vec3) const override {
        LM_UNUSED(geom);
        assert(geom.degenerated);
        return math::pdf_uniform_sphere();
    }

    virtual Vec3 eval(const PointGeometry&, Vec3, Vec3) const override {
        // Normalization constant = 1/(4*pi)
        return Vec3(math::pdf_uniform_sphere());
    }
};

LM_COMP_REG_IMPL(Phase_Isotropic, "phase::isotropic");

LM_NAMESPACE_END(LM_NAMESPACE)
