/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#pragma once

#include "component.h"
#include "math.h"
#include "exception.h"

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

/*!
    \addtogroup volume
    @{
*/

/*!
    \brief Volume data.
*/
class Volume : public Component {
public:
    /*!
        \brief Get bound of the volume.
        \return Bound of the volume.
    */
    virtual Bound bound() const = 0;

    // --------------------------------------------------------------------------------------------

    /*!
        \brief Check if the volume has scalar component.
        \return True if the volume has scalar component.
    */
    virtual bool has_scalar() const = 0;

    /*!
        \brief Evaluate maximum scalar value.
    */
    virtual Float max_scalar() const = 0;

    /*!
        \brief Evaluate scalar value.
        \param p Position in volume coordinates.
        \return Evaluated scalar value.
    */
    virtual Float eval_scalar(Vec3 p) const {
        LM_UNUSED(p);
        LM_THROW_EXCEPTION_DEFAULT(Error::Unimplemented);
    }

    // --------------------------------------------------------------------------------------------

    /*!
        \brief Check if the volume has color component.
        \return True if the volume has color component.
    */
    virtual bool has_color() const = 0;

    /*!
        \brief Evaluate color.
        \param p Position in volume coordinates.
        \return Evaluated color value.
    */
    virtual Vec3 eval_color(Vec3 p) const {
        LM_UNUSED(p);
        LM_THROW_EXCEPTION_DEFAULT(Error::Unimplemented);
    }

    // --------------------------------------------------------------------------------------------

    /*!
        \brief Callback function called for ray marching.
        \param t Distance from the ray origin.
        \retval true Continue raymarching.
        \retval false Abort raymarching.
    */
    using RaymarchFunc = std::function<bool(Float t)>;

    /*!
        \brief March the volume along with the ray.
        \param ray Ray.
        \param tmin Lower bound of the valid range of the ray.
        \param tmax Upper bound of the valid range of the ray.
        \param march_step Ray marching step in world space.
        \param raymarch_func Raymarching function.

        \rst
        This function performs volume-specific raymarching operations.
        \endrst
    */
    virtual void march(Ray ray, Float tmin, Float tmax, Float march_step, const RaymarchFunc& raymarch_func) const {
        LM_UNUSED(ray, tmin, tmax, march_step, raymarch_func);
        LM_THROW_EXCEPTION_DEFAULT(Error::Unimplemented);
    }
};

/*!
    @}
*/

LM_NAMESPACE_END(LM_NAMESPACE)
