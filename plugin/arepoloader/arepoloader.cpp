
#include <pybind11/pybind11.h>
#include <lm/pylm.h>
#include "arepoconfig.h"

#include "arepoloader.h"
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"
#include <memory>
//extern ConfigSet Config;
//extern char ParameterFile[MAXLEN_PATH];
ConfigSet Config;
//LM_NAMESPACE_BEGIN(LM_NAMESPACE)




class  ArepoTempQuantities {
public:
    std::vector<float> vals;
    ArepoTempQuantities () {
        vals.resize(9, 0.0f);
    }
    void clear() {
        for (int i = 0; i < 9; i++)
            vals[i] = 0.0f;
    }
};

class Volume_Arepo_Impl final : public lm::Volume_Arepo {

protected:
    lm::Float scale_;
    lm::Bound bound_;
    lm::Float max_scalar_;

    std::unique_ptr<ArepoMesh>  arepoMesh;
    std::unique_ptr<Arepo> arepo;

    Spectrum s;
    TransferFunction tf;

private:
    void gatherValsAtPoint(lm::Vec3 p, std::vector<float> & toVals) const {
        Point point;
        point.x = p.x;
        point.y = p.y;
        point.z = p.z;

        double minDist = 0.0;
        Vector dir;
        dir.x = 1.0f;
        
        int found = -1;
        int found2 = -1;
        int found3 = -1;
        
        float foundDist = std::numeric_limits<float>::max();
        float foundDist2 = foundDist;
        float foundDist3 = foundDist;
        for(int i = 0; i < arepoMesh->Ndp; i++) {
            lm::Vec3 a = 
            lm::Vec3(arepoMesh->DP[i].x,
            arepoMesh->DP[i].y,
            arepoMesh->DP[i].z);//lm::Vec3(arepoMesh->DTC[i].cx,arepoMesh->DTC[i].cy,arepoMesh->DTC[i].cz);
            float curDist = glm::distance(a, p);
            if ( curDist < foundDist) {
                foundDist = curDist;
                found = arepoMesh->DP[i].index;
            } 
            if (curDist < foundDist2 && curDist > foundDist ) {
                found2 = i;
                foundDist2 = curDist;
            }
            if (curDist < foundDist3 && curDist > foundDist2 ) {
                found3 = i;
                foundDist3 = curDist;
            }
        }
       
        //LM_INFO(std::to_string(found));

        //LM_INFO("{}", partInd );
        if(found > -1 && NumGas > 0) {
            //int index = arepoMesh->getSphPID(found);
            if(found >= NumGas)
                found -= NumGas;
            addValsContribution(toVals,found,1.0 );
        }
        //if(found2 > -1)
        //    addValsContribution(toVals,found2,0.33);
        //if(found3 > -1)
         //   addValsContribution(toVals,found3,0.33);
        

    }
public:
    Volume_Arepo_Impl() : arepoMesh(nullptr), arepo(nullptr), s(0.0f), tf(s) {

    }
    ~Volume_Arepo_Impl() {
        
        
    }

    


    virtual void construct(const lm::Json& prop) override {
        const auto configPath = lm::json::value<std::string>(prop, "configpath");
        auto cutoutPath = lm::json::value<std::string>(prop, "cutoutpath");
        auto pos = cutoutPath.find(".hdf5");
        cutoutPath = cutoutPath.substr(0, pos);
        Config.ReadFile( configPath );
        arepo = std::make_unique<Arepo>(cutoutPath, Config.paramFilename);
        int argc = 0;
        std::string argv0 = "lm?";
        std::string argv1 = configPath;
        std::vector<char*> b = {argv0.data(),argv1.data()};
        char ** arg = &b[0];
        char *** argv = &arg;
        arepo->Init(&argc,argv);
        arepo->LoadSnapshot();
        arepo->ComputeQuantityBounds();
        //only need ArepoMesh implementation (Spectrum and TransferFunction aren't used) 
        s = Spectrum::FromRGB(Config.rgbAbsorb);
        tf = TransferFunction(s);
        arepoMesh = std::make_unique<ArepoMesh>(&tf);
        arepoMesh->ComputeVoronoiEdges();
        
        LM_INFO("constructed Volume Arepo");
         std::cout << " loaded snapshot " << std::endl;
        // Density scale
        scale_ = lm::json::value<lm::Float>(prop, "scale", 1.0f);
        
        LM_INFO("constructed Volume Arepo");
        auto arepoBound = arepoMesh->WorldBound();
        bound_.max = lm::Vec3(arepoBound.pMin.x,arepoBound.pMin.y,arepoBound.pMin.z);
        bound_.min = lm::Vec3(arepoBound.pMax.x,arepoBound.pMax.y,arepoBound.pMax.z);
        //bound
        for (int i = 0; i < arepoMesh->Ndp; i++) {
            auto currentP = lm::Vec3(arepoMesh->DP[i].x,arepoMesh->DP[i].y,arepoMesh->DP[i].z);
            bound_.max = glm::max(bound_.max,currentP);
            bound_.min = glm::min(bound_.min,currentP);
        }
        LM_INFO("spatial bounds {},{},{} ; {},{},{}",bound_.min.x,bound_.min.y,bound_.min.z, bound_.max.x, bound_.max.y, bound_.max.z);
        
        LM_INFO("max scalar  {}",max_scalar());
        LM_INFO("mean scalar  {}",arepo->valBounds[TF_VAL_DENS*3 + 2]);
        LM_INFO("num gas  {}",NumGas);
        

    }

    virtual lm::Bound bound() const override {
        return bound_;        
    }
    virtual lm::Float max_scalar() const override {
        return arepo->valBounds[TF_VAL_DENS*3 + 1];
    }
    virtual bool has_scalar() const override{
        return true;
    }

    virtual lm::Float eval_scalar(lm::Vec3 p) const override {
        
        thread_local ArepoTempQuantities tmpVals1;
        tmpVals1.clear();
        gatherValsAtPoint(p, tmpVals1.vals);
        //if(tmpVals1.vals[TF_VAL_DENS] > 0.0f)
        //   LM_INFO("{}", tmpVals1.vals[TF_VAL_DENS]);
        return scale_ * tmpVals1.vals[TF_VAL_DENS];
    }

    virtual bool has_color() const override {
        return true;
    }
    virtual lm::Vec3 eval_color(lm::Vec3 p) const override {
        thread_local ArepoTempQuantities tmpVals2;
        tmpVals2.clear();
        gatherValsAtPoint(p, tmpVals2.vals);
        return lm::Vec3(1.0,1.0,1.0); //some random stuff 
        //return Vec3(tmpVals2.vals[TF_VAL_BMAG],tmpVals2.vals[TF_VAL_METAL],tmpVals2.vals[TF_VAL_VMAG]); //some random stuff 
        

    }

    //virtual void march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const override {
      //  
    //}




};


LM_COMP_REG_IMPL(Volume_Arepo_Impl, "volume::arepo");


/*


class Mesh_Arepo : public OBJLoaderContext {
public:
    virtual bool load(
        const std::string& path,
        OBJSurfaceGeometry& geo,
        const ProcessMeshFunc& process_mesh,
        const ProcessMaterialFunc& process_material) override
    {
      
        return true;
    }
};

LM_COMP_REG_IMPL(Mesh_Arepo, "objloader::Mesh_Arepo");

*/

//LM_COMP_REG_IMPL(Volume_OpenVDBScalar, "volume::openvdb_scalar");

//LM_NAMESPACE_END(LM_NAMESPACE)
