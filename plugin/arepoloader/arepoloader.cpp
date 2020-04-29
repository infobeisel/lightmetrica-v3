
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

Volume_Arepo::Volume_Arepo() {

}   

Volume_Arepo::~Volume_Arepo() {

}
void Volume_Arepo::construct(const Json& prop)  {

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

Bound Volume_Arepo::bound() const  {
return bound_;
}

Float Volume_Arepo::max_scalar()  const {
return max_scalar_;
}

bool has_scalar()   {
return true;
}

Float Volume_Arepo::eval_scalar(Vec3 p) const  {
//TODO
return 0.0f;
}

bool Volume_Arepo::has_color()  const {
return false;
}

void Volume_Arepo::march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const   {

}

LM_COMP_REG_IMPL(Volume_Arepo, "volume::arepo");


//LM_COMP_REG_IMPL(Volume_OpenVDBScalar, "volume::openvdb_scalar");

LM_NAMESPACE_END(LM_NAMESPACE)
