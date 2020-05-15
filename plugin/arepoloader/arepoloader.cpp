
#include <pybind11/pybind11.h>
#include <lm/pylm.h>
#include "arepoconfig.h"

#include "arepoloader.h"
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"
#include <memory>

#include "lm/stats.h"
#include "statstags.h"

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

    class CachedSample {
        public:
        long long sampleIndex;
        int tetraI;
        int hydroI;
        int minDistI;
        glm::tmat4x3<lm::Float> tetraVs;
        glm::ivec4 tetraInds;
        std::vector<float> values;

        glm::tmat4x3<lm::Float> tmpPVs; //some point to vertex connections
        lm::Vec4 tmpDets; //determinants

        CachedSample() : 
        sampleIndex(std::numeric_limits<long long>::max()),
        tetraI(-1),
        hydroI(-1),
        minDistI(-1),
        values(9,0.0f)
        {
            
        }

    };
}


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
        
        Volume_Arepo_Impl::cachedSample().sampleIndex = std::numeric_limits<long long>::max();
        
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

    static ArepoLoaderInternals::CachedSample & cachedSample() {
        thread_local ArepoLoaderInternals::CachedSample c;
        return c;
    }

    double det3x3(lm::Vec3 b,lm::Vec3 c,lm::Vec3 d) const {
        //glm::determinant(lm::Mat3(b0,b1,b2));
        return b[0]*c[1]*d[2] + c[0]*d[1]*b[2] + d[0]*b[1]*c[2] - d[0]*c[1]*b[2] - c[0]*b[1]*d[2] - b[0]*d[1]*c[2];
    }

    void connectP(glm::tmat4x3<lm::Float> & verts, lm::Vec3 & p, glm::tmat4x3<lm::Float> & pToVerts) const {

        for(int i = 0; i < 4; i++)
            pToVerts[i] = verts[i] - p;
    }

    void computeDeterminants(glm::tmat4x3<lm::Float> & pToVerts, lm::Vec4 & results) const {
        results[0] = det3x3(pToVerts[1],pToVerts[2],pToVerts[3]);
        results[1] = det3x3(pToVerts[0],pToVerts[2],pToVerts[3]);
        results[2] = det3x3(pToVerts[0],pToVerts[1],pToVerts[3]);
        results[3] = det3x3(pToVerts[0],pToVerts[1],pToVerts[2]); 
    }

    bool inside(lm::Vec4 & determinants) const {
        auto ret0 = determinants[0] > 0.0 && determinants[1] < 0.0 && determinants[2] > 0.0 && determinants[3] < 0.0;
        auto ret1 = determinants[0] < 0.0 && determinants[1] > 0.0 && determinants[2] < 0.0 && determinants[3] > 0.0;
        return ret0 || ret1;
    }


    bool insideCachedTetra(lm::Vec3 p) const {
        auto & c = cachedSample();
        connectP(c.tetraVs,p,c.tmpPVs);
        computeDeterminants(c.tmpPVs,c.tmpDets);
        return inside(c.tmpDets);
    }

    //assumes the cached sample has been updated,
    //returns the index to the neighbor tetraeder the queried point is in. returns -1 if it is
    //outside the current and not in any neighbor tetraeder (then it should be outside the arepoMesh)
    int checkCachedNeighbors(lm::Vec3 p) const {
        auto & cached = cachedSample();
        bool checkLessZero = true;
        glm::ivec4 indices = {-1,-1,-1,-1};
        int group0 = 0;
        int group1 = 3; 
        for(int i = 0; i < 4; i++) {
            bool test = checkLessZero ? cached.tmpDets[i] < 0.0 : cached.tmpDets[i] > 0.0;
            int storeTo;
            if(test) {
                storeTo = group0;
                group0++;
            } else {
                storeTo = group1;
                group1--;
            }
            indices[storeTo] = i;
            checkLessZero = !checkLessZero;//switch test every iteration
        }

        //indices now contain the  vertex indices (ranging from 0 to 3) of the tetrahedron
        //it contains two groups, one group where the point lies within, and one group where it lies outside a face
        //the face of vertex i is defined by the triangle of the tetrahedron that doesnt have the vertex i as a corner (i.e. the opposing triangle of vertex i)

        //now indices are in two groups, one group of insides, one of outsides. which is which is unclear, need to find out manually
        assert(group0 != 4 && group0 != 0); //this method wouldnt be called if p was inside the original tetra

        auto copyCached = cachedSample(); //save a copy
        bool foundNeighbor = false;
        bool foundBoundary = false;

        int checkOppositeTetraOfPoint = -1;
        if(group0 == 1 || group0 == 3) { 
            //case 1 or 3 , we have a group with one element only, check this first
            checkOppositeTetraOfPoint = group0 == 1 ? indices[0] : indices[3];
            int neighborTetI = arepoMesh->DT[copyCached.tetraI].t[checkOppositeTetraOfPoint];
            foundBoundary |= neighborTetI < 0; 
            if(neighborTetI >= 0) { //valid tet
                foundNeighbor = insideTetra(neighborTetI,p); //writes into cachedSample
                if(foundNeighbor)
                    return neighborTetI;
            }
            //found: simplest and most efficient case, it is in an opposing tetra
            //not found: then it is in some other neighbor tetra hovering "over a corner" or "over an edge"
        }

        if(!foundNeighbor) {
            std::vector<int> neighbors = {copyCached.tetraI};
            //first test all opposing faces to the current tetra
            for(int i = 0; i < 4  && !foundNeighbor; i++) {
                int ind = indices[i];
                if(ind == checkOppositeTetraOfPoint)
                    continue;//already checked this one
                int neighborTetI = arepoMesh->DT[copyCached.tetraI].t[ind];
                foundBoundary |= neighborTetI < 0; 
                if(neighborTetI >= 0)
                    foundNeighbor = insideTetra(neighborTetI,p);

                if(foundNeighbor)
                    return neighborTetI;

                neighbors.push_back(neighborTetI);
            }
            //still not found, worst case: test neighbors of neighbors
            if(!foundNeighbor) {
                auto alreadyContains = [&] (int tet) {
                    for(auto n : neighbors)
                        if(n == tet)
                            return true;
                    return false;
                };
                auto isNeighbor = [&] (int tet, int of_tet) {
                    int * p1s = arepoMesh->DT[tet].p;
                    int * p2s = arepoMesh->DT[of_tet].p;
                    for(int i = 0; i < 4; i++)
                        for(int j = 0; j < 4; j++)
                            if(p1s[i] == p2s[j])
                                return true;
                    return false;
                };
                std::function<void(int,int)> add_neighbors = [&] (int of_tet, int original_tet) -> void {
                    for(int i = 0; i < 4; i++) {
                        int tetInd = arepoMesh->DT[of_tet].t[i];
                        foundBoundary |= tetInd < 0; 
                        //if the tet is valid, not already added, and a neighbor of our original cached tet
                        if(tetInd >= 0 && !alreadyContains(tetInd) && isNeighbor(tetInd,original_tet)) {
                            neighbors.push_back(tetInd);
                            //also add the neighbors
                            add_neighbors(tetInd, original_tet);
                        }
                    }
                };

                //add neighbors of neighbors
                add_neighbors(neighbors[1], copyCached.tetraI);
                add_neighbors(neighbors[2], copyCached.tetraI);
                add_neighbors(neighbors[3], copyCached.tetraI);
                add_neighbors(neighbors[4], copyCached.tetraI);

                //then perform inside tests 
                //(beginning from 6th element because original and its 4 neighbors were already tested) and cancel when found
                for(int i = 5; i < neighbors.size(); i++) {
                    if(insideTetra(neighbors[i],p)) {
                        return neighbors[i];
                    }
                }

            }

        } 

        //still not found !?
        //then we are outside of the arepo mesh! 
        assert(foundBoundary);
        return -1; 

    }

    //performs point in tetrahedron test, returns result. 
    //always stores all relevant test data into thread local cachedS struct
    bool insideTetra(int tetraIndex, lm::Vec3 p) const {
        assert(tetraIndex >= 0);
        glm::tmat4x3<lm::Float> pVs;//point to vertex connections
        glm::tmat4x3<lm::Float> verts;
        glm::ivec4 vertInds;
        lm::Vec4 determinants;
        
        for(int i = 0; i < 4; i++) {
            vertInds[i] = arepoMesh->DT[tetraIndex].p[i];
            auto av = arepoMesh->DP[vertInds[i]];
            verts[i] = lm::Vec3(av.x,av.y,av.z);
            pVs[i] = verts[i] - p;
        }

        auto & cachedS = Volume_Arepo_Impl::cachedSample();
        cachedS.tetraI = tetraIndex;
        cachedS.tetraVs = verts;
        cachedS.tetraInds = vertInds;

        //skip points 
        if  (
        (arepoMesh->DT[tetraIndex].t[0] < 0 || arepoMesh->DT[tetraIndex].p[0] == DPinfinity || arepoMesh->DT[tetraIndex].p[1] == DPinfinity
        || arepoMesh->DT[tetraIndex].p[2] == DPinfinity || arepoMesh->DT[tetraIndex].p[3] == DPinfinity)
        || arepoMesh->DT[tetraIndex].t[0] == -1)
        {
            //LM_INFO("skip");
            return false;
        }
        
        connectP(verts,p,pVs);
        computeDeterminants(pVs,determinants);
        bool insideTet = inside(determinants);
        
        cachedS.tmpPVs = pVs;
        cachedS.tmpDets = determinants;
        return insideTet;
        
    }

    void evaluateDensityCached( std::vector<float> & toVals) const {
        auto & cachedS = Volume_Arepo_Impl::cachedSample();
        auto lengths =  lm::Vec4(
            glm::length2(cachedS.tmpPVs[0]),
            glm::length2(cachedS.tmpPVs[1]),
            glm::length2(cachedS.tmpPVs[2]),
            glm::length2(cachedS.tmpPVs[3]));

        int minDistIndex = lengths[0] < lengths[1] ? 0 : 1;
        minDistIndex = lengths[minDistIndex] > lengths[2] ?  2 : minDistIndex;
        minDistIndex = lengths[minDistIndex] > lengths[3] ?  3 : minDistIndex;
        int hydroIndex = 0;
        if(cachedS.minDistI == minDistIndex) {
            hydroIndex = cachedS.hydroI;   
            toVals[TF_VAL_DENS] = cachedS.values[TF_VAL_DENS];
        } else {//need to look up in arepomesh structure
            hydroIndex = arepoMesh->DP[cachedS.tetraInds[minDistIndex]].index;
            if(hydroIndex > -1 && hydroIndex  < NumGas &&  NumGas > 0) {
                for (int i = 0; i < 9; i++)
                    cachedS.values[i] = 0.0f;
                addValsContribution(cachedS.values,hydroIndex,1.0);//lengths[minDistIndex] / totalD );
                toVals[TF_VAL_DENS] = cachedS.values[TF_VAL_DENS];
            } 
        }                
        cachedS.hydroI = hydroIndex;
        cachedS.minDistI = minDistIndex;
    }


    void gatherValsAtPoint(lm::Vec3 p, std::vector<float> & toVals) const {

        int h = 0;
        auto currentSample = lm::stats::get<int,int,long long>(h);
        
        auto & cachedS = Volume_Arepo_Impl::cachedSample();
        if(cachedS.sampleIndex == currentSample) { //is cached
            if (insideCachedTetra(p)) {
                evaluateDensityCached(toVals);
                return; //most efficient case, still in cached tetra
            }
            else { //check neighbors
                int tetraIndex = checkCachedNeighbors(p);
                bool inside = tetraIndex >= 0;
                if(inside) { //sample changed to a neighbor
                    cachedS.sampleIndex = currentSample;
                    //need to invalidate some information
                    cachedS.hydroI = -1;
                    cachedS.minDistI = -1;
                    //then evaluate
                    evaluateDensityCached(toVals);
                    return;
                } else {// not even in neighbors, "lost" track
                    cachedS.sampleIndex = std::numeric_limits<long long>::max();
                }
            }
        }
        //uncached, need to ray intersect with volume
        lm::Ray r; //TODO to be further developed, with real ray!
        r.o = p;
        r.d = lm::Vec3(1.0f,0.0f,0.0f); //arbitrary
        auto hit = accel->intersect(r,0.0f,std::numeric_limits<float>::max());
        if(hit.has_value()) { //check if inside the tetra of hit triangle
            int tetraIndex = hit.value().face / 4;
            bool inside = insideTetra(tetraIndex, p);
            //can be outside because we have double triangles in mesh, 
            //having exact same positions but belonging to two opposing tetrahedra,
            //need to check for both.
            if(!inside) {
                tetraIndex = checkCachedNeighbors(p);
                inside = tetraIndex >= 0;
            }
           
            if(inside) { // we found a tetrahedron where we are inside
                cachedS.sampleIndex = currentSample;
                //need to invalidate some information
                cachedS.hydroI = -1;
                cachedS.minDistI = -1;
                evaluateDensityCached(toVals);
                return;
            } else { //still not inside, then we are outside the whole arepoMesh at the moment, next sample will have to perform ray intersection test again.  
                cachedS.sampleIndex = std::numeric_limits<long long>::max();
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
