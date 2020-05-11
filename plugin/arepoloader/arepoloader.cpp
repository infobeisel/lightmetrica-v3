
#include <pybind11/pybind11.h>
#include <lm/pylm.h>
#include "arepoconfig.h"

#include "arepoloader.h"
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"
#include <memory>


#include <lm/lm.h>
#include <lm/core.h>
#include <lm/scene.h>
#include <lm/json.h>
#include <lm/jsontype.h>
#include <lm/mesh.h>
#include <lm/accel.h>
//extern ConfigSet Config;
//extern char ParameterFile[MAXLEN_PATH];
ConfigSet Config;



namespace ArepoLoaderInternals {

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
}

//LM_NAMESPACE_BEGIN(LM_NAMESPACE)

class ArepoMeshImpl final : public lm::Mesh {
private:
    point const * dps;
    size_t ndp;
    tetra const * dts;
    size_t ndt;
  
public:
    //TODO!
    //LM_SERIALIZE_IMPL(ar) {
    //    ar(dps, ndp, dts, ndt);
   // }

public:
    virtual void construct(const lm::Json& prop) override {

        
        //std::uintptr_t ptrToPoints = prop["ps_addr"];
        //dps = reinterpret_cast<point*>(ptrToPoints);
        dps = prop["ps_addr"].get<const point*>();
        ndp = prop["ps_count"];
        
        //std::uintptr_t ptrToTetras = prop["ts_addr"];
        //dts = reinterpret_cast<tetra*>(ptrToPoints);;
        dts = prop["ts_addr"].get<const tetra*>();
        ndt = prop["ts_count"];

        
        
    }

    virtual void foreach_triangle(const ProcessTriangleFunc& process_triangle) const override {
        for(size_t fi = 0; fi < ndt * 4; fi++) {
            process_triangle(fi, triangle_at(fi));
        }
    }

    virtual lm::Mesh::Tri triangle_at(int face) const override {
        size_t tetra = face / 4;
        size_t faceIndex = face % 4;// (face - tetra * 4);

        //"oriented" tetrahedron points...?
        size_t vertexIndex0 = (0 + faceIndex) % 4;
        size_t vertexIndex1 = (1 + faceIndex) % 4;
        size_t vertexIndex2 = (2 + faceIndex) % 4;
        
        auto v0 = lm::Vec3(dps[dts[tetra].p[vertexIndex0]].x, dps[dts[tetra].p[vertexIndex0]].y, dps[dts[tetra].p[vertexIndex0]].z);
        auto v1 = lm::Vec3(dps[dts[tetra].p[vertexIndex1]].x, dps[dts[tetra].p[vertexIndex1]].y, dps[dts[tetra].p[vertexIndex1]].z);
        auto v2 = lm::Vec3(dps[dts[tetra].p[vertexIndex2]].x, dps[dts[tetra].p[vertexIndex2]].y, dps[dts[tetra].p[vertexIndex2]].z);
        
        auto triangleNormal = glm::cross (
            glm::normalize(v1 - v0),
            glm::normalize(v2 - v0));
        auto uv = lm::Vec2(0); 

        return { 
            {v0,triangleNormal, uv},
            {v1,triangleNormal, uv},
            {v2,triangleNormal, uv},
        };
    }

    virtual lm::Mesh::InterpolatedPoint surface_point(int face, lm::Vec2 uv) const override {

        auto triangle = triangle_at(face);
        auto triangleNormal = triangle.p1.n;
        return {
            lm::math::mix_barycentric(triangle.p1.p, triangle.p2.p, triangle.p3.p, uv),
            glm::normalize(lm::math::mix_barycentric(triangle.p1.n, triangle.p2.n, triangle.p3.n, uv)),
            lm::math::geometry_normal(triangle.p1.p, triangle.p2.p, triangle.p3.p),
            lm::math::mix_barycentric(triangle.p1.t, triangle.p2.t, triangle.p3.t, uv)
        };
    }

    virtual int num_triangles() const override {
        return ndt * 4;
    }
};

LM_COMP_REG_IMPL(ArepoMeshImpl, "mesh::arepo");



class Volume_Arepo_Impl final : public lm::Volume_Arepo {

    
    public:
    
    Volume_Arepo_Impl() : arepoMesh(nullptr), arepo(nullptr), s(0.0f), tf(s),scene(nullptr) {

    }
    ~Volume_Arepo_Impl() {
        
        
    }

    virtual Component* underlying(const std::string& name) const override {
        if (name == "tetramesh") {
            return meshAdapter.get();
        } else if(name == "tetraaccel") {
            return accel.get();
        } else if(name == "tetrascene") {
            return scene.get();
        } else if(name == "dummymat") {
            return dummyMat.get();
        }
        return nullptr;

    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        lm::comp::visit(visit, meshAdapter);
        lm::comp::visit(visit, accel);
        lm::comp::visit(visit, scene);
        lm::comp::visit(visit, dummyMat);
        
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
        

        meshAdapter = lm::comp::create<lm::Mesh>( //lm::load<lm::Mesh>( 
        "mesh::arepo", make_loc("tetramesh"), {
            {"ps_addr", (const point*)arepoMesh->DP},
            //{"ps_addr", reinterpret_cast<std::uintptr_t>(arepoMesh->DP)},
            {"ps_count", arepoMesh->Ndp},
            {"ts_addr", (const tetra*)arepoMesh->DT},
            //{"ts_addr", reinterpret_cast<std::uintptr_t>(arepoMesh->DT)},
            {"ts_count", arepoMesh->Ndt}
        });
        
        
        /*meshAdapter->foreach_triangle([&] (int face, const lm::Mesh::Tri& t) {
            LM_INFO("tri nr {} at {}{}{}",face,t.p1.p[0],t.p1.p[1],t.p1.p[2]);
        });*/
        
        dummyMat = lm::comp::create<lm::Material>(
        "material::diffuse", make_loc("dummymat"), {
            {"Kd", lm::Vec3(0.0f)}
        });

        accel = lm::comp::create<lm::Accel>("accel::sahbvh", make_loc("tetraaccel"), {});
        LM_INFO( accel->loc());
        scene = lm::comp::create<lm::Scene>("scene::default", make_loc("tetrascene"), {
            {"accel", accel->loc()}
        });
        
        
        auto c = scene->create_primitive_node({
             {"mesh" , meshAdapter->loc()},
             {"material" , dummyMat->loc()}
        });
        scene->add_child(scene->root_node(), c);
        
        scene->build();
        accel->build(*scene.get());
        
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
        
        thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals1;
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
        thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals2;
        tmpVals2.clear();
        gatherValsAtPoint(p, tmpVals2.vals);
        return lm::Vec3(1.0,1.0,1.0); //some random stuff 
        //return Vec3(tmpVals2.vals[TF_VAL_BMAG],tmpVals2.vals[TF_VAL_METAL],tmpVals2.vals[TF_VAL_VMAG]); //some random stuff 
        

    }

    //virtual void march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const override {
      //  
    //}

    protected:
    lm::Float scale_;
    lm::Bound bound_;
    lm::Float max_scalar_;

    std::unique_ptr<ArepoMesh>  arepoMesh;
    std::unique_ptr<Arepo> arepo;

    Spectrum s;
    TransferFunction tf;

    private:    

    lm::Component::Ptr<lm::Mesh> meshAdapter;
    lm::Component::Ptr<lm::Accel> accel;
    lm::Component::Ptr<lm::Scene> scene;
    lm::Component::Ptr<lm::Material> dummyMat;

    float det3x3(lm::Vec3 b0,lm::Vec3 b1,lm::Vec3 b2) const {
        glm::determinant(glm::mat3(b0,b1,b2));
    }

    bool insideTetra(lm::Vec3 p, size_t tetraIndex,lm::Vec3 const & v0,lm::Vec3 const & v1,lm::Vec3 const & v2,lm::Vec3 const & v3) const {
        
        auto pa = p - v0;
        auto pb = p - v1;
        auto pc = p - v2;
        auto pd = p - v3;
        
        auto da = det3x3(pb,pc,pd);
        auto db = det3x3(pa,pc,pd);
        auto dc = det3x3(pb,pa,pd);
        auto dd = det3x3(pb,pc,pa);
        return (da < 0.0f && db > 0.0f && dc < 0.0f && dd > 0.0f)
            || (da > 0.0f && db < 0.0f && dc > 0.0f && dd < 0.0f);
    }

    void gatherValsAtPoint(lm::Vec3 p, std::vector<float> & toVals) const {
       
        lm::Ray r;
        r.o = p;
        r.d = lm::Vec3(1.0f,0.0f,0.0f); //arbitrary
        
        auto hit = accel->intersect(r,0.0f,std::numeric_limits<float>::max());
        if(hit.has_value()) { //check if inside the tetra of hit triangle
            size_t tetraIndex = hit.value().face / 4;
            glm::tmat4x3<lm::Float> pVs;//point to vertex connections
            glm::tmat4x3<lm::Float> verts;
            glm::ivec4 vertInds;
            for(int i = 0; i < 4; i++) {
                vertInds[i] = arepoMesh->DT[tetraIndex].p[0];
                auto av = arepoMesh->DP[vertInds[i]];
                verts[i] = lm::Vec3(av.x,av.y,av.z);
                pVs[i] = p - verts[i];
            }
            

            if (insideTetra(p,  tetraIndex, verts[0], verts[1], verts[2], verts[3])) {
                auto lengths =  lm::Vec4(
                    glm::length(pVs[0]),
                    glm::length(pVs[1]),
                    glm::length(pVs[2]),
                    glm::length(pVs[3]));
                auto totalD = lengths.x + lengths.y + lengths.z + lengths.w;
                //TODO weird check, saw it in arepo vtk project
                for(int i = 0; i < 4; i ++) {
                    if(vertInds[i] >= NumGas && NumGas > 0)
                        vertInds[i] -= NumGas;
                    addValsContribution(toVals,vertInds[i],lengths[i] / totalD);
                }
            }
        }
    }

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
