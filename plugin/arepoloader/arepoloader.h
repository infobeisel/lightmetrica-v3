#include <lm/core.h>
#include <lm/volume.h>
LM_NAMESPACE_BEGIN(LM_NAMESPACE)

class Volume_Arepo : public Volume {
private:
    Float scale_;
    Bound bound_;
    Float max_scalar_;

public:
    Volume_Arepo() ;
    ~Volume_Arepo();

    

public:
    virtual void construct(const Json& prop) override ;

    virtual Bound bound() const override ;
    virtual Float max_scalar() const override ;

    virtual bool has_scalar() const override ;

    virtual Float eval_scalar(Vec3 p) const override ;

    virtual bool has_color() const override ;

    virtual void march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const override ;
};

LM_NAMESPACE_END(LM_NAMESPACE)
