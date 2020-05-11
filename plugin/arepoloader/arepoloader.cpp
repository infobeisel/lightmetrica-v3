
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
    point * dps;
    size_t ndp;
    tetra * dts;
    size_t ndt;

    std::vector<lm::Mesh::Tri> triangles;
  
public:
    //TODO!
    //LM_SERIALIZE_IMPL(ar) {
    //    ar(dps, ndp, dts, ndt);
   // }

public:
    virtual void construct(const lm::Json& prop) override {

        
        //std::uintptr_t ptrToPoints = prop["ps_addr"];
        //dps = reinterpret_cast<point*>(ptrToPoints);
        
        dps = reinterpret_cast<point*> ( prop["ps_addr"].get<uintptr_t>() );
        ndp = prop["ps_count"];
        
        //std::uintptr_t ptrToTetras = prop["ts_addr"];
        //dts = reinterpret_cast<tetra*>(ptrToPoints);;
        dts = reinterpret_cast<tetra*> ( prop["ts_addr"].get<uintptr_t>() );
        ndt = prop["ts_count"];


        triangles.reserve(ndt * 4);
        for(size_t tetrai = 0; tetrai < ndt; tetrai++) {

            //"oriented" tetrahedron points...?
           // size_t vertexIndex0 = (0 + faceIndex) % 4;
            //ize_t vertexIndex1 = (1 + faceIndex) % 4;
           // size_t vertexIndex2 = (2 + faceIndex) % 4;
           // size_t vertexIndex3 = (2 + faceIndex) % 4;
            
           //auto v3 = lm::Vec3(dps[dts[tetrai].p[2]].x, dps[dts[tetrai].p[2]].y, dps[dts[tetrai].p[2]].z);
            
            for(int i = 0; i < 4; i++) {
                int a = (i + 0) % 4;
                int b = (i + 1) % 4;
                int c = (i + 2) % 4;
               
                auto v0 = lm::Vec3(dps[dts[tetrai].p[a]].x, dps[dts[tetrai].p[a]].y, dps[dts[tetrai].p[a]].z);
                auto v1 = lm::Vec3(dps[dts[tetrai].p[b]].x, dps[dts[tetrai].p[b]].y, dps[dts[tetrai].p[b]].z);
                auto v2 = lm::Vec3(dps[dts[tetrai].p[c]].x, dps[dts[tetrai].p[c]].y, dps[dts[tetrai].p[c]].z);
                auto triangleNormal = glm::cross (glm::normalize(v1 - v0),glm::normalize(v2 - v0));
                triangles.push_back({ 
                    {v0,triangleNormal,lm::Vec2(v0.x)},
                    {v1,triangleNormal, lm::Vec2(v0.y)},
                    {v2,triangleNormal, lm::Vec2(v0.z)},
                });
            
            /*LM_INFO("add tetrahedron {} ? {} . triangle {},{},{} with v0x {}" ,std::to_string(tetrai), 
            std::to_string(dts[tetrai].t[0] == -1),
            std::to_string(a),std::to_string(b),std::to_string(c),
            std::to_string(dps[dts[tetrai].p[a]].x) );*/
            }
        }

        
        
    }

    virtual void foreach_triangle(const ProcessTriangleFunc& process_triangle) const override {
        /*for(size_t fi = 0; fi < ndt * 4; fi++) {
            process_triangle(fi, triangle_at(fi));
        }*/
        for(int t = 0; t <num_triangles();t++)
            process_triangle(t,triangles[t]);
    }

    virtual lm::Mesh::Tri triangle_at(int face) const override {
        return triangles[face];
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


/*
class Mesh_Arepo_Impl final : public lm::Mesh {
    private:

    std::unique_ptr<ArepoMesh>  arepoMesh;
    std::unique_ptr<Arepo> arepo;

    Spectrum s;
    TransferFunction tf;

    lm::Component::Ptr<lm::Mesh> mesh;

    public:
    Mesh_Arepo_Impl() : arepoMesh(nullptr), arepo(nullptr),s(0.0f), tf(s){};
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
    

        mesh = lm::comp::create<lm::Mesh>( //lm::load<lm::Mesh>( 
        "mesh::arepo", make_loc("tetramesh"), {
            {"ps_addr", reinterpret_cast<uintptr_t> (arepoMesh->DP)},
            //{"ps_addr", reinterpret_cast<std::uintptr_t>(arepoMesh->DP)},
            {"ps_count", arepoMesh->Ndp},
            {"ts_addr", reinterpret_cast<uintptr_t> (arepoMesh->DT)},
            //{"ts_addr", reinterpret_cast<std::uintptr_t>(arepoMesh->DT)},
            {"ts_count", arepoMesh->Ndt}
        });

        
    }

    virtual void foreach_triangle(const ProcessTriangleFunc& process_triangle) const override {
        mesh->foreach_triangle(process_triangle);
    }

    virtual lm::Mesh::Tri triangle_at(int face) const override {
        return mesh->triangle_at(face);
    }

    virtual lm::Mesh::InterpolatedPoint surface_point(int face, lm::Vec2 uv) const override {

        return mesh->surface_point(face,uv);
    }

    virtual int num_triangles() const override {
        return mesh->num_triangles();
    }
};

LM_COMP_REG_IMPL(Mesh_Arepo_Impl, "mesh::fromarepo");
*/

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
        naive = lm::json::value<int>(prop,"naive",0) == 1 ? true : false;
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
        scale_ = lm::json::value<lm::Float>(prop, "scale", 1.0);
        
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
            //{"ps_addr", (const point*)arepoMesh->DP},
            {"ps_addr", reinterpret_cast<uintptr_t>(arepoMesh->DP)},
            {"ps_count", arepoMesh->Ndp},
            //{"ts_addr", (const tetra*)arepoMesh->DT},
            {"ts_addr", reinterpret_cast<uintptr_t>(arepoMesh->DT)},
            {"ts_count", arepoMesh->Ndt}
        });
        
        
        /*meshAdapter->foreach_triangle([&] (int face, const lm::Mesh::Tri& t) {
            LM_INFO("tri nr {} at {}{}{}",face,t.p1.p[0],t.p1.p[1],t.p1.p[2]);
        });*/
        
        dummyMat = lm::comp::create<lm::Material>(
        "material::diffuse", make_loc("dummymat"), {
            {"Kd", lm::Vec3(0.0f)}
        });

        
        accel = lm::comp::create<lm::Accel>("accel::embree", make_loc("tetraaccel"), {});
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
        //accel->build(*scene.get());
        
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
        if(naive)
            gatherValsAtPointNaive(p, tmpVals1.vals);
        else
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
        //gatherValsAtPoint(p, tmpVals2.vals);
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
    bool naive;

    double det3x3(lm::Vec3 b,lm::Vec3 c,lm::Vec3 d) const {
        //glm::determinant(lm::Mat3(b0,b1,b2));
        return b[0]*c[1]*d[2] + c[0]*d[1]*b[2] + d[0]*b[1]*c[2] - d[0]*c[1]*b[2] - c[0]*b[1]*d[2] - b[0]*d[1]*c[2];
    }

    bool insideTetra(glm::vec3 v0,glm::vec3 v1,glm::vec3 v2,glm::vec3 v3,glm::vec3 p) const {
        auto a = v0 - p;
        auto b = v1 - p;
        auto c = v2 - p;
        auto d = v3 - p;
        auto detA = det3x3(b,c,d);
        auto detB = det3x3(a,c,d);
        auto detC = det3x3(a,b,d);
        auto detD = det3x3(a,b,c);
        auto ret0 = detA > 0.0 && detB < 0.0 && detC > 0.0 && detD < 0.0;
        auto ret1 = detA < 0.0 && detB > 0.0 && detC < 0.0 && detD > 0.0;
        return ret0 || ret1;

    }

    void gatherValsAtPoint(lm::Vec3 p, std::vector<float> & toVals) const {
       
        lm::Ray r;
        r.o = p;
        r.d = lm::Vec3(1.0f,0.0f,0.0f); //arbitrary
        
        auto hit = accel->intersect(r,0.0f,std::numeric_limits<float>::max());
        if(hit.has_value()) { //check if inside the tetra of hit triangle
            int tetraIndex = hit.value().face / 4;
            glm::tmat4x3<lm::Float> pVs;//point to vertex connections
            glm::tmat4x3<lm::Float> verts;
            glm::ivec4 vertInds;
            for(int i = 0; i < 4; i++) {
                vertInds[i] = arepoMesh->DT[tetraIndex].p[i];
                auto av = arepoMesh->DP[vertInds[i]];
                verts[i] = lm::Vec3(av.x,av.y,av.z);
                pVs[i] = verts[i] - p;
            }
            //skip points 
            if  (
            (arepoMesh->DT[tetraIndex].t[0] < 0 || arepoMesh->DT[tetraIndex].p[0] == DPinfinity || arepoMesh->DT[tetraIndex].p[1] == DPinfinity
            || arepoMesh->DT[tetraIndex].p[2] == DPinfinity || arepoMesh->DT[tetraIndex].p[3] == DPinfinity)
            || arepoMesh->DT[tetraIndex].t[0] == -1)
            {
                //LM_INFO("skip");
                return;
            }


            if (insideTetra( verts[0], verts[1], verts[2], verts[3], p)) {
                //toVals[TF_VAL_DENS] = 1.0f;
                //LM_INFO("inside");

                auto lengths =  lm::Vec4(
                    glm::length(pVs[0]),
                    glm::length(pVs[1]),
                    glm::length(pVs[2]),
                    glm::length(pVs[3]));

                int minDistIndex = lengths[0] < lengths[1] ? 0 : 1;
                minDistIndex = lengths[minDistIndex] > lengths[2] ?  2 : minDistIndex;
                minDistIndex = lengths[minDistIndex] > lengths[3] ?  3 : minDistIndex;
                
                
                //auto totalD = lengths.x + lengths.y + lengths.z + lengths.w;
                //TODO weird check, saw it in arepo vtk project
                /*float foundDist = std::numeric_limits<float>::max();
                int found = -1;
                for(int i = 0; i < 4; i++) {
                    float curDist = glm::distance(p, verts[i]);
                    if ( curDist < foundDist) {
                        foundDist = curDist;
                        found = arepoMesh->DP[vertInds[i]].index;
                    } 
                
                }
                if(found > 0 && hydroIndex >= NumGas && NumGas > 0)
                        hydroIndex -= NumGas;
                addValsContribution(toVals,hydroIndex,lengths[i] / totalD);
                */
                //for(int i = 0; i < 4; i ++) {
                // skip those with initial points outside the box or connecting to DPinfinity
                

                int hydroIndex = arepoMesh->DP[vertInds[minDistIndex]].index;
                
                if(hydroIndex > -1 && hydroIndex  < NumGas &&  NumGas > 0) {
                    addValsContribution(toVals,hydroIndex,1.0);//lengths[minDistIndex] / totalD );
                } 
                //}
            }
        }
    }

    void gatherValsAtPointNaive(lm::Vec3 p, std::vector<float> & toVals) const {
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
        /*meshAdapter->foreach_triangle( [&] ( int face, const lm::Mesh::Tri & tr) {
            float curDist = glm::min(
                glm::min(glm::distance(p,tr.p1.p),glm::distance(p,tr.p2.p)),
                glm::distance(p,tr.p3.p));
            if ( curDist < foundDist) {
                foundDist = curDist;
                found = face;
            } 
        });*/
        for(int i = 0; i < arepoMesh->Ndp; i++) {
            lm::Vec3 a = 
            lm::Vec3(arepoMesh->DP[i].x,
            arepoMesh->DP[i].y,
            arepoMesh->DP[i].z);//lm::Vec3(arepoMesh->DTC[i].cx,arepoMesh->DTC[i].cy,arepoMesh->DTC[i].cz);
            float curDist = glm::distance(a, p);
            if ( curDist < foundDist && i != DPinfinity) {
                foundDist = curDist;
                found = arepoMesh->DP[i].index;
            } 
          
        }
            
        //}
       
        //LM_INFO(std::to_string(found));

      
        if(found > -1 && NumGas > 0 && found < NumGas) {
            //int index = arepoMesh->getSphPID(found);
            //if(found >= NumGas)
             //   found -= NumGas;
            addValsContribution(toVals,found,1.0 );
        }
        //if(found2 > -1)
        //    addValsContribution(toVals,found2,0.33);
        //if(found3 > -1)
         //   addValsContribution(toVals,found3,0.33);
        

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
