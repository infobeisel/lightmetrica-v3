/*
	Lightmetrica - Copyright (c) 2019 Hisanari Otsu
	Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/core.h>
#include <lm/volume.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

/*
\rst
.. function:: volume::gaussian

    Sphere volume.

    :param color color: Stored color value.
    :param float scalar: Stored max scalar value.
    :param Vec3 pos: center position of the volume
    :param Vec3 sigma: gaussian parameter for x,y,z
\endrst
*/
class Volume_Sphere : public Volume {
private:
	Bound bound_;//is computed from 0.001% of the max of gaussian
	std::optional<Vec3> color_;
	std::optional<Float> scalar_;
    Vec3 pos_;
    Float radius_;


public:
	LM_SERIALIZE_IMPL(ar) {
		ar(bound_);
	}

public:
	virtual void construct(const Json& prop) override {
		color_ = json::value_or_none<Vec3>(prop, "color");
        pos_ = json::value<Vec3>(prop, "pos", Vec3(.0_f,.0_f,.0_f));
		scalar_ = json::value_or_none<Float>(prop, "scalar");
        radius_ = json::value<Float>(prop, "radius", 1._f);
		if (!color_ && !scalar_) {
            LM_THROW_EXCEPTION(Error::InvalidArgument,
			    "Either 'color' or 'scalar' property is necessary.");
		}
        computeBound();
	}

    void computeBound(){
        
        bound_.max = pos_+sqrt(3)*radius_;
        bound_.min = pos_-sqrt(3)*radius_;

        LM_INFO("min {}, {}, {}",bound_.min.x,bound_.min.y,bound_.min.z);
        LM_INFO("max {}, {}, {}",bound_.max.x,bound_.max.y,bound_.max.z);
         
    }
	
	virtual Bound bound() const override {
		return bound_;
	}

	virtual bool has_scalar() const override {
		return scalar_.has_value();
	}

	virtual Float max_scalar() const override {
		return *scalar_;
	}

    bool inSphere(Vec3 p, Float r) const{
        return (glm::length(p)<r)?true:false;
    }

	virtual Float eval_scalar(Vec3 p) const override {
		return max_scalar()*inSphere(pos_-p, radius_);
	}

	virtual bool has_color() const override {
		return color_.has_value();
	}

	virtual Vec3 eval_color(Vec3) const override {
		return *color_;
	}
};

LM_COMP_REG_IMPL(Volume_Sphere, "volume::sphere");

LM_NAMESPACE_END(LM_NAMESPACE)
