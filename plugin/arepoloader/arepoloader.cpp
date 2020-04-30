
#include <pybind11/pybind11.h>
#include <lm/pylm.h>

#include <snapio.h>
#include <unistd.h>
#include <png.h>
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"
#include <functional>
#include <vector>
#include "math.h"
#include <experimental/random>

#include "arepoloader.h"
LM_NAMESPACE_BEGIN(LM_NAMESPACE)



class Volume_Arepo_Impl final : public Volume_Arepo {

    private:
    Float scale_;
    Bound bound_;
    Float max_scalar_;

public:
    Volume_Arepo_Impl() ;
    ~Volume_Arepo_Impl();

    


    virtual void construct(const Json& prop) override {

        const auto configPath = json::value<std::string>(prop, "configpath");
        const auto cutoutPath = json::value<std::string>(prop, "cutoutpath");
        std::cout << "hello Arepo2OpenVDB" << std::endl;
        ConfigSet Config;
        Config.ReadFile( configPath );
        Arepo arepo = Arepo(cutoutPath, Config.paramFilename);
        int argc = 0;

        arepo.Init(&argc,nullptr);
        arepo.LoadSnapshot();
        arepo.ComputeQuantityBounds();
        //only need ArepoMesh implementation (Spectrum and TransferFunction aren't used) 
        const Spectrum s = Spectrum::FromRGB(Config.rgbAbsorb);
        TransferFunction tf(s);
        ArepoMesh * arepoMesh = new ArepoMesh(&tf);
        arepoMesh->ComputeVoronoiEdges();
        
        std::cout << " loaded snapshot " << std::endl;
        // Density scale
        scale_ = json::value<Float>(prop, "scale", 1_f);

        // Bound

        // Maximum density
        LM_INFO("constructed Volume Arepo");
    }

    virtual Bound bound() const override {
return bound_;
}
    virtual Float max_scalar() const override {
return max_scalar_;
}
    virtual bool has_scalar() const override{
return true;
}

    virtual Float eval_scalar(Vec3 p) const override {
//TODO
return 0.0f;
}

    virtual bool has_color() const override {
return false;
}

    virtual void march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const override {

}




};


LM_COMP_REG_IMPL(Volume_Arepo_Impl, "volume::arepo");


//LM_COMP_REG_IMPL(Volume_OpenVDBScalar, "volume::openvdb_scalar");

LM_NAMESPACE_END(LM_NAMESPACE)
