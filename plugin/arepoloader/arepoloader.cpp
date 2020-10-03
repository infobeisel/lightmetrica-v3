
#include <pybind11/pybind11.h>
#include <lm/pylm.h>
#include "arepoconfig.h"

#include "arepoloader.h"
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"

#include <memory>

#include "pluecker.h"

#include "lm/stats.h"
#include "statstags.h"

#include <lm/lm.h>
#include <lm/core.h>
#include <lm/scene.h>
#include <lm/json.h>
#include <lm/jsontype.h>
#include <lm/mesh.h>
#include <lm/accel.h>



#include <algorithm>
//extern ConfigSet Config;
//extern char ParameterFile[MAXLEN_PATH];
ConfigSet Config;





namespace ArepoLoaderInternals {

    #define DENSITY_CRANKUP 1.0
    #define REJECTION_SAMPLES_COUNT 10
    #define RAY_SEGMENT_ALLOC 200


    static lm::Float MODEL_SCALE = 1.0;

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

    struct ArepoMeshWrapper : public IArepoMeshMock {


        ArepoMesh * ref;
        std::vector<lm::Float> emptyDensities;
        ArepoMeshWrapper(ArepoMesh*m) : ref(m) {}

        virtual tetra *  getDT() override {return ref->DT;}
        virtual point * getDP() override {return ref->DP;}
        virtual std::vector<lm::Float> & getdensities() {LM_ERROR("getdensities: not implemented");return emptyDensities;}
        virtual int getNdt() {return ref->Ndt;}
        virtual int getNdp() {return ref->Ndp;}

        virtual BBox WorldBound() override {return ref->WorldBound();}
        virtual lm::Float getDensity(int index) {LM_ERROR("getDensity: not implemented");return -1;}
        virtual lm::Float max_density() {LM_ERROR("max_density: not implemented");return -1;}

    
    };
    struct ArepoMeshMock : public IArepoMeshMock {
        std::vector<tetra> DT;
        std::vector<point> DP;
        std::vector<lm::Float> densities;
        int Ndt;
        int Ndp;

        virtual tetra * getDT() override {return &DT[0];}
        virtual point * getDP() override {return &DP[0];}
        virtual std::vector<lm::Float> & getdensities() override {return densities;}
        virtual int getNdt() override {return Ndt;}
        virtual int getNdp() override {return Ndp;}

        virtual BBox WorldBound() override  {
            BBox b;
            b.pMin.x = b.pMin.y = b.pMin.z = std::numeric_limits<lm::Float>::max();
            b.pMax.x = b.pMax.y = b.pMax.z = - std::numeric_limits<lm::Float>::max();
            
            for (auto p : DP) {
                b.pMin.x = std::min(b.pMin.x, p.x);b.pMin.y = std::min(b.pMin.y, p.y);b.pMin.z = std::min(b.pMin.z, p.z);
                b.pMax.x = std::max(b.pMax.x, p.x);b.pMax.y = std::max(b.pMax.y, p.y);b.pMax.z = std::max(b.pMax.z, p.z);
            }
            return b;
        }
        //ArepoMeshMock(const lm::Json& prop)  {

        //    density_ = json::value<Float>(prop, "density");
        //    albedo_ = json::value<Vec3>(prop, "albedo");
        
        //}
        ArepoMeshMock()  {

             

            /*
            point A;A.x=-1;A.y=1;A.z=1;A.index=0;
            point B;B.x=-1;B.y=-1;B.z=1;B.index=1;
            point C;C.x=1;C.y=-1;C.z=1;C.index=2;
            point D;D.x=1;D.y=1;D.z=1;D.index=3;
            point E;E.x=-1;E.y=1;E.z=-1;E.index=4;
            point F;F.x=-1;F.y=-1;F.z=-1;F.index=5;
            point G;G.x=1;G.y=-1;G.z=-1;G.index=6;
            point H;H.x=1;H.y=1;H.z=-1;H.index=7;

            DP.resize(0);
            DP.reserve(8);
            DP.push_back(A);
            DP.push_back(B);
            DP.push_back(C);
            DP.push_back(D);
            DP.push_back(E);
            DP.push_back(F);
            DP.push_back(G);
            DP.push_back(H);
            
            densities.resize(0);
            densities.reserve(8);
            densities.push_back(0.00000001);
            densities.push_back(0.1);
            densities.push_back(0.00000001);
            densities.push_back(0.00000001);
            densities.push_back(0.00000001);
            densities.push_back(0.00000001);
            densities.push_back(0.00000001);
            densities.push_back(0.00000001);

            tetra ACFH;
            ACFH.p[0] = 0;ACFH.p[1] = 2;ACFH.p[2] = 5;ACFH.p[3] = 7;
            ACFH.t[0] = 4;ACFH.t[1] = 2;ACFH.t[2] = 3;ACFH.t[3] = 1; //has a neighbor in every direction
            //ACFH.t[0] = 5;ACFH.t[1] = 5;ACFH.t[2] = 5;ACFH.t[3] = 1; //has a neighbor in every direction
            //ACFH.t[0] = 5;ACFH.t[1] = 5;ACFH.t[2] = 5;ACFH.t[3] = 5;
             
            tetra ABCF;
            ABCF.p[0] = 0;ABCF.p[1] = 1;ABCF.p[2] = 2;ABCF.p[3] = 5;
            ABCF.t[0] = 5;ABCF.t[1] = 0;ABCF.t[2] = 5;ABCF.t[3] = 5;//5 leads to the deleted tetra->no neighbor

            tetra AEHF;
            AEHF.p[0] = 0;AEHF.p[1] = 4;AEHF.p[2] = 7;AEHF.p[3] = 5;
            AEHF.t[0] = 5;AEHF.t[1] = 0;AEHF.t[2] = 5;AEHF.t[3] = 5;

            tetra ACDH;
            ACDH.p[0] = 0;ACDH.p[1] = 2;ACDH.p[2] = 3;ACDH.p[3] = 7;
            ACDH.t[0] = 5;ACDH.t[1] = 5;ACDH.t[2] = 0;ACDH.t[3] = 5;

            tetra FCGH;
            FCGH.p[0] = 5;FCGH.p[1] = 2;FCGH.p[2] = 6;FCGH.p[3] = 7;
            FCGH.t[0] = 5;FCGH.t[1] = 5;FCGH.t[2] = 0;FCGH.t[3] = 5;
            
            
            tetra tDel;
            tDel.t[0] = -1;
            tDel.t[1] = -1;
            tDel.t[2] = -1;
            tDel.t[3] = -1;
            tDel.p[0] = 0;
            tDel.p[1] = 0;
            tDel.p[2] = 0;
            tDel.p[3] = 0;
          
            DT.resize(0);
            DT.reserve(6);
            DT.push_back(ACFH);
            DT.push_back(ABCF);
            DT.push_back(AEHF);
            DT.push_back(ACDH);
            DT.push_back(FCGH);
            DT.push_back(tDel);

            Ndt = 6;
            
            Ndp = 8;*/


            {//ACDG
                tetra t;
               
                t.t[0] = 5;
                t.t[1] = 4;
                t.t[2] = 5;
                t.t[3] = 5;


                t.p[0] = 0;
                t.p[1] = 2;
                t.p[2] = 3;
                t.p[3] = 6;
                t.s[0] = -1;
                t.s[1] = -1;
                t.s[2] = -1;
                t.s[3] = -1;

                DT.push_back(t);
            }

            {//DFGH
                tetra t;
               
                t.t[0] = 5;
                t.t[1] = 5;
                t.t[2] = 5;
                t.t[3] = 4;


                t.p[0] = 3;
                t.p[1] = 5;
                t.p[2] = 6;
                t.p[3] = 7;
                t.s[0] = -1;
                t.s[1] = -1;
                t.s[2] = -1;
                t.s[3] = -1;

                DT.push_back(t);
            }

            {//AEFG
                tetra t;
               
                t.t[0] = 5;
                t.t[1] = 4;
                t.t[2] = 5;
                t.t[3] = 5;


                t.p[0] = 0;
                t.p[1] = 4;
                t.p[2] = 5;
                t.p[3] = 6;
                t.s[0] = -1;
                t.s[1] = -1;
                t.s[2] = -1;
                t.s[3] = -1;

                DT.push_back(t);
            }

            {//ABFG
                tetra t;
               
                t.t[0] = 5;
                t.t[1] = 4;
                t.t[2] = 5;
                t.t[3] = 5;


                t.p[0] = 0;
                t.p[1] = 1;
                t.p[2] = 5;
                t.p[3] = 3;
                t.s[0] = -1;
                t.s[1] = -1;
                t.s[2] = -1;
                t.s[3] = -1;

                DT.push_back(t);
            }

            {//AFDG
                tetra t;
               
                t.t[0] = 1;
                t.t[1] = 0;
                t.t[2] = 2;
                t.t[3] = 3;


                t.p[0] = 0;
                t.p[1] = 5;
                t.p[2] = 3;
                t.p[3] = 6;
                t.s[0] = -1;
                t.s[1] = -1;
                t.s[2] = -1;
                t.s[3] = -1;

                DT.push_back(t);
            }

            

       

            point p;

            p.x = -10;p.y = 10;p.z = 10;p.index = 0;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = 10;p.y = 10;p.z = 10;p.index=1;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = -10;p.y = -10;p.z = 10;p.index=2;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = 10;p.y = -10;p.z = 10;p.index=3;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = -10;p.y = 10;p.z = -10;p.index=4;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = 10;p.y = 10;p.z = -10;p.index = 5; //F
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = -10;p.y = -10;p.z = -10;p.index=6;
            DP.push_back(p);
            densities.push_back(0.0000000000003);

            p.x = 10;p.y = -10;p.z = -10;p.index=7;
            DP.push_back(p);
            densities.push_back(0.0000000000003);


           

            /*{
                tetra t;
                t.t[0] = 3;
                t.t[1] = 3;
                t.t[2] = 3;
                t.t[3] = 0;

                t.p[0] = 1;
                t.p[1] = 2;
                t.p[2] = 3;
                t.p[3] = 5;
                t.s[0] = 1;
                t.s[1] = 1;
                t.s[2] = 1;
                t.s[3] = 1;

                DT.push_back(t);
            }*/


            tetra tDel;
            tDel.t[0] = -1;
            tDel.t[1] = -1;
            tDel.t[2] = -1;
            tDel.t[3] = -1;
            tDel.p[0] = -1;
            tDel.p[1] = -1;
            tDel.p[2] = -1;
            tDel.p[3] = -1;
            tDel.s[0] = 0;
            tDel.s[1] = 0;
            tDel.s[2] = 0;
            tDel.s[3] = 0;
            DT.push_back(tDel);

            Ndt = 6;
            Ndp = 8;


        }

        virtual lm::Float getDensity(int index) override  {

            //for(int i = 0; i < densities.size(); i++) {
            //    LM_INFO("density {}: {}, ask {} ",i, densities[i], index);
           // }
            return densities[index] / MODEL_SCALE;
        }

        virtual lm::Float max_density() override  {
            return (*std::max_element(densities.begin(), densities.end())) / MODEL_SCALE;
        }
    };

#ifdef MOCK_AREPO
    static IArepoMeshMock * arepoMeshRef = nullptr;
#else
    static IArepoMeshMock * arepoMeshRef = nullptr;
#endif


    static lm::Accel *  accelRef = nullptr;


    class CachedSample {
        public:
        long long sampleIndex; //the id of the sample that is currently cached
        int tetraI;
        int hydroI;
        int minDistI;
        glm::tmat4x3<lm::Float> tetraVs;
        glm::ivec4 tetraInds;
        std::vector<float> values;

        glm::tmat4x3<lm::Float> tmpPVs; //some point to vertex connections
        lm::Vec4 tmpDets; //determinants

        //T^-1, a barycentric coordinates matrix
        lm::Float mainDeterminant; //the tetrahedron's determinant
        lm::Mat3 baryInvT;
        std::vector<std::vector<float>> cornerVals;

        std::vector<tetra> neighbors;
        std::vector<int> neighborInds;

        lm::Vec4 dirDets;

        int looksAtTriId;//the tetrahedron's triangle index (0-3) the current sample ray looks at.
        lm::Accel::Hit lastHit;
        lm::Mesh::Tri lastHitTri;

        //lm::Float max_scalar_;


        CachedSample() : 
        sampleIndex(std::numeric_limits<long long>::max()),
        tetraI(-1),
        hydroI(-1),
        minDistI(-1),
        values(9,0.0f),
        cornerVals(4,{0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f}),
        looksAtTriId(0),
        neighbors(100,tetra()),
        neighborInds(100,0)

        {
            
        }

    };



    
    static std::vector<lm::RaySegmentCDF> & raySegments() {
        int i = 0;
        auto & c =  lm::stats::getRef<int,int,std::vector<lm::RaySegmentCDF>>(i); 
        return c;
    }



    struct UniformSampleCache {
        std::vector<lm::Float> urs;//uniform random samples
        std::vector<lm::Float> acccdf;//accumulated cdfs
        std::vector<lm::Float> freet;//free paths
        //std::vector<bool> scattered; //the sample already scattered
        
        
        lm::Float minfreet() {
            lm::Float ret = std::numeric_limits<lm::Float>::max();
            for(auto & t : freet)
                ret = glm::min(ret,t);
            return ret;
        }

    };

    inline double det3x3(lm::Vec3 b,lm::Vec3 c,lm::Vec3 d) {
        //glm::determinant(lm::Mat3(b0,b1,b2));
        return b[0]*c[1]*d[2] + c[0]*d[1]*b[2] + d[0]*b[1]*c[2] - d[0]*c[1]*b[2] - c[0]*b[1]*d[2] - b[0]*d[1]*c[2];
    }
    
    inline void connectP(glm::tmat4x3<lm::Float> const  & verts, lm::Vec3 const & p, glm::tmat4x3<lm::Float> & pToVerts) {

        for(int i = 0; i < 4; i++)
            pToVerts[i] = verts[i] - p;
    }

    inline void computeDeterminants(glm::tmat4x3<lm::Float> & pToVerts, lm::Vec4 & results) {
        results[0] = det3x3(pToVerts[1],pToVerts[2],pToVerts[3]);
        results[1] = det3x3(pToVerts[0],pToVerts[2],pToVerts[3]);
        results[2] = det3x3(pToVerts[0],pToVerts[1],pToVerts[3]);
        results[3] = det3x3(pToVerts[0],pToVerts[1],pToVerts[2]); 
    }

    //this only works iff the tetrahedron is specified with a specific vertex order, namely, see mesh_arepo.cpp !!!
    inline bool insideOld(lm::Vec4 & determinants) {
        auto ret0 = determinants[0] > 0.0 && determinants[1] < 0.0 && determinants[2] > 0.0 && determinants[3] < 0.0;
        auto ret1 = determinants[0] < 0.0 && determinants[1] > 0.0 && determinants[2] < 0.0 && determinants[3] > 0.0;
        return ret0 || ret1;
    }

    //this test works irrespective of vertex order of the tetrahedron
    inline bool inside(lm::Vec4 & determinants, lm::Float tetraDeterminant) {
        return abs(abs(tetraDeterminant) - (abs(determinants[0])+abs(determinants[1])+abs(determinants[2])+abs(determinants[3]))) < INSIDE_TOLERANCE;       
    }


    inline bool insideCachedTetra(lm::Vec3 p, CachedSample & c) {
        connectP(c.tetraVs,p,c.tmpPVs);
        computeDeterminants(c.tmpPVs,c.tmpDets);
        return inside(c.tmpDets, c.mainDeterminant);
    }

    inline void sampleCachedScalarCoefficients(lm::Ray ray, lm::Float & a, lm::Float &  b, CachedSample const & cached) {
        //this method assumes that the ray resides inside the current tetrahedron
        //therefore it clamps the barycentric coordinates in the b term to 1
        
        auto lambda012_a = cached.baryInvT * ray.d;
        auto lambda012_b = cached.baryInvT * (ray.o - cached.tetraVs[3]);

        auto maxd = glm::max(
            cached.cornerVals[0][TF_VAL_DENS] ,glm::max(
            cached.cornerVals[1][TF_VAL_DENS] ,glm::max(
            cached.cornerVals[2][TF_VAL_DENS] ,
            cached.cornerVals[3][TF_VAL_DENS] )));
        
        auto densities = lm::Vec4(
            cached.cornerVals[0][TF_VAL_DENS] ,
            cached.cornerVals[1][TF_VAL_DENS] ,
            cached.cornerVals[2][TF_VAL_DENS] ,
            cached.cornerVals[3][TF_VAL_DENS] );

       

        //densities *= scale_;
        //auto lambda012_a = glm::transpose(cached.baryInvT) * lm::Vec3(densities);
        //lambda012_a = lambda012_a - glm::transpose(cached.baryInvT) * lm::Vec3(densities[3]);

        b = glm::dot( lm::Vec4(lambda012_b, 1.0 - lambda012_b.x - lambda012_b.y - lambda012_b.z), densities);
        a = glm::dot( lambda012_a, lm::Vec3(densities));
        a +=  densities[3] * (- lambda012_a.x - lambda012_a.y - lambda012_a.z);
        
        //a = 0.001;
        //b = 0.000001;

        //LM_INFO("{} and {}",a,b);
        if(b<=0.0) {
           // LM_ERROR("b is non positive?!");
           // LM_INFO("b is non positive?!");

             if(densities.x <= 0.0 &&
        densities.y <= 0.0 &&
        densities.z <= 0.0 && 
        densities.w <= 0.0 ) {
          //          LM_INFO("densities dont have positive value ?!");
            
           //     LM_INFO("dens {},{},{},{}",cached.cornerVals[0][TF_VAL_DENS],
           ///     cached.cornerVals[1][TF_VAL_DENS],
             //   cached.cornerVals[2][TF_VAL_DENS],
             //   cached.cornerVals[3][TF_VAL_DENS]);
            
        }

        }
        //if(a < 0.0) {
        //    LM_ERROR("a is negative?!");
        //    LM_INFO("a is negative?!");
        //}
        
        
    }


    inline lm::Float sampleTransmittance(lm::Ray ray, lm::Float fromT, lm::Float toT,lm::Float a, lm::Float b, CachedSample const & cached) {
        return glm::exp(- sampleCDF(toT,a,b));
    }




    inline lm::Float intersectCachedTetra(lm::Ray ray, CachedSample & cached, bool usepluecker)  {

        //assumes ray.o is within cached tetra!


        
        
        //find out in which of the the 4 sub tetras (spanned by ray.o and the current tetrahedron )
        //the ray direction lies, determine vertex indices that represent the triangle the ray dir looks at.
        //choose tetra vertex 3 as "roof" R, test certain determinants (dir,X,R),
        // where negative means dir is "right of", positive "left of" axis XR
        glm::ivec3 indices; 

        cached.dirDets[0] = det3x3(ray.d,cached.tmpPVs[1],cached.tmpPVs[3]);
        bool det_B_D =  cached.mainDeterminant < 0.0 ? cached.dirDets[0] > 0.0 : cached.dirDets[0] < 0.0;
        indices.x = det_B_D ? 0 : 2; //A or C is guaranteed first vertex of final triangle

        cached.dirDets[1] = det3x3(ray.d,cached.tmpPVs[indices.x],cached.tmpPVs[3]);
        bool det_AorC_D = cached.mainDeterminant < 0.0 ? cached.dirDets[1] < 0.0 : cached.dirDets[1] > 0.0;
        indices.y = det_B_D   ? 
        (det_AorC_D ? 1 : 2) : 
        (det_AorC_D ? 0 : 1);


        cached.dirDets[2] = det3x3(ray.d,cached.tmpPVs[indices.x],cached.tmpPVs[indices.y]);
        //left of horizontal x to y ?
        bool det_horizontal = cached.mainDeterminant < 0.0 ? cached.dirDets[2] < 0.0 :cached.dirDets[2] > 0.0;
        //the top or the bottom?
        indices.z = 
        det_B_D ? 
            (det_AorC_D ?
                (det_horizontal ? 2 : 3) : //A->B->?
                (det_horizontal ? 3 : 1)   //A->C->?
            ) :
            (det_AorC_D ?
                (det_horizontal ? 1 : 3) : //C->A->?
                (det_horizontal ? 3 : 0)   //C->B->?
            )
        ;
        
        auto faceid = det_B_D ? 
            (det_AorC_D ?
                (det_horizontal ? 0 : 3) :
                (det_horizontal ? 2 : 0) 
            ) :
            (det_AorC_D ?
                (det_horizontal ? 0 : 2) : 
                (det_horizontal ? 1 : 0) 
            )
        ;

        //indices now contains 3 indices into tetra vertices forming the triangle the ray dir is looking at.
        //can compute point on triangle (intersection) using some determinants (->barycentric coordinates)
        //main determinant, assuming tmpDets contain the correct values
        auto mainD = det_B_D ? 
            (det_AorC_D ?
                (det_horizontal ? cached.tmpDets[3] : cached.tmpDets[2]) :
                (det_horizontal ? cached.tmpDets[1] : cached.tmpDets[3]) 
            ) :
            (det_AorC_D ?
                (det_horizontal ? cached.tmpDets[3] : cached.tmpDets[1]) : 
                (det_horizontal ? cached.tmpDets[0] : cached.tmpDets[3])
            )
        ;

        //TODO fix this intersection method to work with all kinds of vertex orders :/

        
        
        //construct determinants : ray direction (assume it is normalized) and PVs
        auto p = ray.o + ray.d;
        glm::tmat4x3<lm::Float> pVs;//point to vertex connections
        glm::tmat4x3<lm::Float> verts;
        connectP(cached.tetraVs,p, pVs);


        //TODO already better but still having t issues O:O
        lm::Float volumeRelation = std::abs(mainD) + std::abs(det3x3(pVs[indices[0]],pVs[indices[1]], pVs[indices[2]]));
        volumeRelation =  volumeRelation / std::abs(mainD); 
        //ray.d *= abs(mainD);
        lm::Float bary0 = std::abs( det3x3(cached.tmpPVs[indices[1]],cached.tmpPVs[indices[2]], ray.d) );
        lm::Float bary1 = std::abs( det3x3(cached.tmpPVs[indices[0]],cached.tmpPVs[indices[2]], ray.d) );
        lm::Float bary2 = std::abs( det3x3(cached.tmpPVs[indices[0]],cached.tmpPVs[indices[1]], ray.d) );
        lm::Float baryOrigin = 1.0 - bary0 - bary1 - bary2;
        
        auto t =  abs(mainD) / (bary0 + bary1 + bary2); 
        //auto t = (1.0 / volumeRelation); 
        //so this is the t the ray dir can be multiplied to, without leaving the tetra
        //ray.o + ray.d * t is intersection point.
        
        //LM_INFO("{} until intersection with tetra bounds ", t);


/*
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j< 3; j++) {
                
               LM_INFO("tmppvs {}",cached.tmpPVs[i][j]);
            }

            
        }

               LM_INFO("maindet {}",mainD);*/
        cached.looksAtTriId = faceid;

        if(!usepluecker) {
            //test plucker instead
            return t;

        } else {
        
        
            lm::Float plueckert = 0.0;
            std::vector<lm::Vec3> vertList = {cached.tetraVs[0],cached.tetraVs[1],cached.tetraVs[2],cached.tetraVs[3]};
            int enterface, leaveface;
            lm::Vec3 enterpoint, leavepoint;
            lm::Float uenter1, uenter2, uleave1, uleave2, tenter, tleave;
            bool retval = RayTetraPluecker(ray.o, ray.d,
                vertList.data(),
                enterface, leaveface,
                enterpoint,leavepoint,
                uenter1, uenter2, uleave1,uleave2,
                tenter, tleave);
                // LM_INFO(" {} t enter and leave {}, {} vs {}",retval,tenter,tleave, accelT  );
            if(retval) {
                if(tleave < tenter) { 
                    if(tleave < 0.0 && tenter > 0.0)  // tenter is valid
                        plueckert = tenter;
                    if(tleave > 0.0) // tleave is valud
                        plueckert = tleave;
                }
            }
            return plueckert;
        }

        
    }


    //copied from ArepoVTK
    // get primary hydro IND - handle local ghosts
    inline int getCorrectedHydroInd(int dpInd)
    {
        int SphP_ID = -1;
  
        if (dpInd >= 0 && dpInd < NumGas)
            SphP_ID = dpInd;
        else if (dpInd >= NumGas)
            SphP_ID = dpInd - NumGas;
        
        if (SphP_ID < 0)
            LM_ERROR("1135");
            
        return SphP_ID;
    }



    inline void updateCachedBaryInvT(CachedSample & cachedS)  {
         //need to calculate a barycentric coordinates matrix 
        lm::Float a = cachedS.tetraVs[0].x - cachedS.tetraVs[3].x;
        lm::Float d = cachedS.tetraVs[0].y - cachedS.tetraVs[3].y;
        lm::Float g = cachedS.tetraVs[0].z - cachedS.tetraVs[3].z;

        lm::Float b = cachedS.tetraVs[1].x - cachedS.tetraVs[3].x;
        lm::Float e = cachedS.tetraVs[1].y - cachedS.tetraVs[3].y;
        lm::Float h = cachedS.tetraVs[1].z - cachedS.tetraVs[3].z;
        
        lm::Float c = cachedS.tetraVs[2].x - cachedS.tetraVs[3].x;
        lm::Float f = cachedS.tetraVs[2].y - cachedS.tetraVs[3].y;
        lm::Float i = cachedS.tetraVs[2].z - cachedS.tetraVs[3].z;
        //https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices

        lm::Float A =   e * i - f * h; 
        lm::Float B =   -(d * i - f * g);
        lm::Float C =   d * h - e * g;//m0[1] * m1[2] - m1[1] * m0[2]; 

        lm::Float D =   -(b * i - c * h); //- m1[0] * m2[2] + m2[0] * m1[2]; 
        lm::Float E =   a * i - c * g;//m0[0] * m2[2] - m2[0] * m0[2]; 
        lm::Float F =   -(a * h - b * g);//- m0[0] * m1[2] + m1[0] * m0[2]; 
        
        lm::Float G =   b * f - c * e;//m1[0] * m2[1] - m2[0] * m1[1]; 
        lm::Float H =   -(a * f - c * d); //- m0[0] * m2[1] + m2[0] * m0[1]; 
        lm::Float I =   a * e - b * d;//m0[0] * m1[1] - m1[0] * m0[1];


        //calculate the inverse matrix that delivers the barycoordinates (from wikipedia)
        //already computed at this point cachedS.mainDeterminant =  (a * A + b * B + c * C);
        cachedS.baryInvT = 
        lm::Mat3(lm::Vec3(A,B,C),lm::Vec3(D,E,F),lm::Vec3(G,H,I))
        / (a*A+b*B+c*C); //cachedS.mainDeterminant;

    }

    
    inline void cacheCornerValues(CachedSample & cachedS)  {

        for(int i = 0; i < 4; i++) {
            for (int j = 0; j < 9; j++)
                cachedS.cornerVals[i][j] = 0.0f;
            int num = arepoMeshRef->getDP()[ cachedS.tetraInds[i]].index;
#ifdef MOCK_AREPO
            cachedS.cornerVals[i][TF_VAL_DENS] = arepoMeshRef->getDensity(num);
#else
            if(num >= 0) {
                int hydroIndex = getCorrectedHydroInd(num);
                
                if(hydroIndex > -1 && num  < NumGas) {
                    addValsContribution(cachedS.cornerVals[i],hydroIndex,DENSITY_CRANKUP/MODEL_SCALE);//lengths[minDistIndex] / totalD );
                }  
            }
#endif
        }
        
         
    }

    inline bool insideTetra(int tetraIndex,lm::Vec3 const & loc, tetra const & tetra,
        glm::tmat4x3<lm::Float> & verts, CachedSample & cachedS)  {

            glm::tmat4x3<lm::Float> pVs;//point to vertex connections
            //glm::tmat4x3<lm::Float> verts;
            glm::ivec4 vertInds;
            lm::Vec4 determinants;


            for(int i = 0; i < 4; i++) {
                vertInds[i] = tetra.p[i];
                //auto av =  arepoMeshRef->getDP()[vertInds[i]];
                //verts[i] = MODEL_SCALE * lm::Vec3(av.x,av.y,av.z);
                pVs[i] = verts[i] - loc;
            }

            //also transports sign (choose vertex 3 as "roof"): negative means ccw, positive cw tetrahedron definition
            lm::Float mainDeterminant = det3x3(verts[0] - verts[3],verts[1] - verts[3],verts[2] - verts[3]);
            //if(mainDeterminant > 0.0)
                //   LM_INFO( "cw");


            //skip points 
            //wtf see arepo vtk
            DPinfinity = -1;
            if  (
            (tetra.t[0] < 0 ||tetra.p[0] <= DPinfinity || tetra.p[1] <= DPinfinity
            || tetra.p[2] <= DPinfinity || tetra.p[3] <= DPinfinity)
            || tetra.t[0] == -1)
            {
                //LM_INFO("skip");
                return false;
            }



            connectP(verts,loc,pVs);
            computeDeterminants(pVs,determinants);
            bool insideTet = inside(determinants,mainDeterminant);
            if(insideTet) {
                cachedS.tetraI = tetraIndex;
                cachedS.tetraVs = verts;
                cachedS.tetraInds = vertInds;
                cachedS.tmpPVs = pVs;
                cachedS.tmpDets = determinants;
                cachedS.mainDeterminant = mainDeterminant;
                updateCachedBaryInvT(cachedS);
                cachedS.hydroI = -1;
                cachedS.minDistI = -1;
                cacheCornerValues(cachedS);
            }

            return insideTet;

    }
    

    //performs point in tetrahedron test, returns result. 
    //stores all relevant test data into cachedS if p is in the tetra
    inline bool insideTetra(int tetraIndex, tetra const & tetra, lm::Vec3 const & p, CachedSample & cachedS)  {
        glm::tmat4x3<lm::Float> verts;
        for(int i = 0; i < 4; i++) {
            auto av =  arepoMeshRef->getDP()[tetra.p[i]];
            verts[i] = MODEL_SCALE * lm::Vec3(av.x,av.y,av.z);
        }
        return insideTetra(tetraIndex,p,tetra,verts,cachedS);
        
    }

    //overload
    inline bool insideTetra(int tetra, lm::Vec3 const & p, CachedSample & cachedS)  {
        return insideTetra(tetra, arepoMeshRef->getDT()[tetra], p, cachedS);
    }

    inline void updateCachedNeighbors(CachedSample & cached,  lm::ArepoLMMesh * toQueryTetraId) {
        auto ofTetra = cached.tetraI;
        bool foundNeighbor = false;
        bool foundBoundary = false;
        auto & neighbors = cached.neighbors; 
        auto & neighborInds = cached.neighborInds; 
        neighbors.clear();
        neighborInds.clear();

        
        

        auto alreadyContains = [&] (int tet) {
            for(auto n : neighborInds)
                if(n == tet)
                    return true;
            return false;
        };
        auto isNeighbor = [&] (int tet, int of_tet) {
            int * p1s = arepoMeshRef->getDT()[tet].p;
            int * p2s = arepoMeshRef->getDT()[of_tet].p;
            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++)
                    if(p1s[i] == p2s[j])
                        return true;
            return false;
        };

        auto isNeighbor2 = [&] (int tet, int of_tet) {
            bool foundN2 = isNeighbor(tet, of_tet);
            for(int i = 0;i < 4; i++) {
                int firstNeighbor = arepoMeshRef->getDT()[tet].t[i];
                if(firstNeighbor >= 0)
                    foundN2 = foundN2 || isNeighbor(firstNeighbor,tet);
            }
            return foundN2;
            
        };

        std::function<void(int,int)> addNeighbors = [&] (int of_tet, int original_tet) -> void {
            auto tetStrct = arepoMeshRef->getDT()[of_tet]; 
            int n0 = tetStrct.t[0];
            foundBoundary = foundBoundary || n0 < 0; 
            if(n0 >= 0) {
                neighbors.push_back(tetStrct);
                neighborInds.push_back(of_tet);

                for(int i = 0;i < 4; i++) {
                    int tetInd = tetStrct.t[i];
                    if(tetInd >= 0 && !alreadyContains(tetInd) && isNeighbor(tetInd,original_tet)) {
                        //also add the neighbors
                        addNeighbors(tetInd, original_tet);
                    }
                }
            }
        };
        auto * ps = arepoMeshRef->getDT()[ofTetra].p;
        auto * ts = arepoMeshRef->getDT();
        //auto mostprobable = arepoMeshRef->getDT()[ofTetra].t[(cached.looksAtTriId - 1) % 4];

        for(int p = 0; p < 4; p++) {
            auto & adjacents = toQueryTetraId->adjacentTs(ps[p]);
            for(auto & i : adjacents) {
                neighborInds.push_back(i);
                neighbors.push_back(ts[i]);
            }
        }

        

        //neighbors.push_back(mostprobable);
        //neighborInds.push_back(mostprobable.t[(cached.looksAtTriId - 1) % 3]);
        /*if(mostprobable >= 0 && mostprobable < arepoMeshRef->getNdt() && arepoMeshRef->getDT()[mostprobable].t[0] >= 0) {
            neighbors.push_back(arepoMeshRef->getDT()[mostprobable]);
            neighborInds.push_back(mostprobable);
            addNeighbors(mostprobable,ofTetra);
        }*/
            
        //else
        //    addNeighbors(ofTetra,ofTetra);
        
    
    }

    //assumes the cached sample has been updated,
    //returns the index to the neighbor tetraeder the queried point is in. returns -1 if it is
    //outside the current and not in any neighbor tetraeder (then it should be outside the arepoMeshRef)
    inline int checkCachedNeighbors(lm::Vec3 p, CachedSample & cached)  {
        auto ofTetra = cached.tetraI;
        bool foundNeighbor = false;
        bool foundBoundary = false;
        //double check if in original tetra
        auto & neighbors = cached.neighbors;     
        for(int i = 0; i < neighbors.size(); i ++) {
            foundNeighbor = insideTetra(cached.neighborInds[i], neighbors[i],p, cached);
            if(foundNeighbor) {
                return cached.neighborInds[i];
            }
        }


        //still not found !?
        //then we are outside of the arepo mesh! 
        //if(!foundBoundary)
        //   LM_ERROR("tetra mesh traversal: didnt find neighbor but was not at boundary");
        return -1; 

    }

   

    inline bool findAndCacheTetra( CachedSample & cachedS, lm::Vec3 p, lm::Vec3 dir, lm::ArepoLMMesh * toQueryTetraId)  {
        int h = 0;

        //gather a guess that might have been set before calling this method!
        auto guess = lm::stats::get<lm::stats::TetraIdGuess,int,int>(h);


        auto currentSample = lm::stats::get<lm::stats::CachedSampleId,int,long long>(h);
        lm::stats::add<lm::stats::TotalTetraTests,int,long long>(h,1);

        bool returnValue = false;

        if(cachedS.sampleIndex == currentSample) { //is cached
            lm::stats::add<lm::stats::SampleIdCacheHits,int,long long>(h,1);
            if (insideCachedTetra(p,cachedS)) {
                //evaluateDensityCached(toVals);
                lm::stats::add<lm::stats::UsedCachedTetra,int,long long>(h,1);
                return true; //most efficient case, still in cached tetra
            }
            else {
                //check guess
                bool inside = false;
                int tetraIndex = -1;
                //guess serves for jumped query positions. 
                if(guess >= 0) { 
                    inside = insideTetra(guess, p, cachedS) ;//|| insideTetra(ni, p, cachedS);
                    tetraIndex = guess;
                }

                //for continuous query positions, probably it changed to one of the neighbours!
                if(!inside) {
                    //check neighbors
                    tetraIndex = checkCachedNeighbors(p,cachedS);
                    inside = tetraIndex >= 0;
                }
                if(inside) { //sample changed to a neighbor
                    updateCachedNeighbors(cachedS,toQueryTetraId);
                    lm::stats::add<lm::stats::UsedNeighborTetra,int,long long>(h,1);
                    //need to invalidate some information
                    cachedS.hydroI = -1;
                    cachedS.tetraI = tetraIndex;
                    cachedS.minDistI = -1;
                    //also store density at corners
                    cacheCornerValues(cachedS);
                    //add it to the queue of traversed tetras
                    return true;

                } else {
                    //we totally lost the tetrahedron, need to make request to acceleration structure
                }
            }
        } else {
            cachedS.sampleIndex = std::numeric_limits<long long>::max();

            if(guess >= 0) { //try guess
                if(insideTetra(guess, p, cachedS)){
                    updateCachedNeighbors(cachedS,toQueryTetraId);
                    //need to invalidate some information
                    cachedS.hydroI = -1;
                    //also store density at corners
                    cacheCornerValues(cachedS);
                    return true;
                }
            }
            //empty cached traversal queue
            //cachedTraversal().clear();
        }
        
        lm::stats::add<lm::stats::ResampleAccel,int,long long>(h,1);
        //uncached, need to ray intersect with volume
        lm::Ray r; 
        r.o = p;
        r.d = dir;
        auto hit = accelRef->intersect(r,0.0,std::numeric_limits<lm::Float>::max());

        if (hit != std::nullopt && hit.value().face >= 0 && hit.value().face < toQueryTetraId->num_triangles()) { //check if inside the tetra of hit triangle
            //LM_INFO("hit sth");
            cachedS.lastHit = hit.value(); //save hit value
            int localFaceIndex = (hit.value().face ) % 4;
            cachedS.looksAtTriId = localFaceIndex;
            int tetraIndex =  toQueryTetraId->correspondingTetra(hit.value().face);
            cachedS.lastHitTri =  toQueryTetraId->triangle_at(hit.value().face);
            //auto ni = arepoMeshRef->DT[tetraIndex].t[localFaceIndex];
            bool inside = insideTetra(tetraIndex, p, cachedS) ;//|| insideTetra(ni, p, cachedS);
            //LM_INFO("inside {} : {}", tetraIndex, inside);
            //can be outside because we have double triangles in mesh, 
            //having exact same positions but belonging to two opposing tetrahedra,
            //need to check for both.
            //check neighbors
           // if(!inside) {
            //    LM_INFO("test shit?");
            //}
           for(int i = 0; i < 4; i++) {
                auto ni = arepoMeshRef->getDT()[tetraIndex].t[i];
                if(ni >= 0 && arepoMeshRef->getDT()[ni].t[0] >= 0 ) { //only check if the neighbor was not deleted
                    inside = inside ||  insideTetra(ni, p, cachedS);
                   // LM_INFO("inside {} : {}", tetraIndex, inside);
                }
            }
            //if(!inside) {
            //    LM_INFO("we really should be inside somewhere now!");
            //}
            
            if(inside) { // we found a tetrahedron where we are inside
                updateCachedNeighbors(cachedS,toQueryTetraId);
                //LM_INFO("found tetra inside");
                cachedS.sampleIndex = currentSample;
                //need to invalidate some information
                cachedS.hydroI = -1;
                cachedS.minDistI = -1;
               //LM_ERROR("this doesnt happend unfortunately?!");
                cacheCornerValues(cachedS);
                //add traversed tetra
                //cachedTraversal().push_back(cachedS);

                returnValue = true;
            }
            if(!inside) { 
                //LM_INFO("found intersection but not corresponding tetra");
                //cachedS.lastHit.t = std::numeric_limits<lm::Float>::infinity(); 
                returnValue = false;
                //this is a weird case: we have hit sth with the intersection test but 
                //didnt find the corresponding tetrahedron...!?
            }
        } else {
            //cachedTraversal().clear();

            returnValue = false;
            cachedS.lastHit.t = std::numeric_limits<lm::Float>::infinity(); 
        }
        
        return returnValue;
    }


}

using namespace ArepoLoaderInternals;

class Volume_Arepo_Impl final : public lm::Volume_Arepo {

    
    public:
    
    Volume_Arepo_Impl() : arepoMesh(nullptr), arepo(nullptr), s(0.0f), tf(s),gastetrascene(nullptr),lightscene(nullptr) {

    }
    ~Volume_Arepo_Impl() {
        
        
    }

    virtual Component* underlying(const std::string& name) const override {
        if (name == "tetramesh") {
            return  dynamic_cast<lm::Mesh*>(meshAdapter.get());
            
        } else if(name == "tetraaccel") {
            return tetraaccel.get();
        } else if(name == "tetrascene") {
            return gastetrascene.get();
        } else if(name == "dummymat") {
            return dummyMat.get();
        } else if(name == "lightscene") {
            return lightscene.get();
        } else if(name == "lighthierarchy") {
            return lighthierarchy.get();
        }


        return nullptr;

    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        lm::comp::visit(visit, meshAdapter);
        lm::comp::visit(visit, tetraaccel);
        lm::comp::visit(visit, gastetrascene);
        lm::comp::visit(visit, lightscene);
        lm::comp::visit(visit, lighthierarchy);
        lm::comp::visit(visit, dummyMat);
        
    }

    


    virtual void construct(const lm::Json& prop) override {
        const auto configPath = lm::json::value<std::string>(prop, "configpath");
        tetraTestMargin = lm::json::value<lm::Float>(prop, "tetrahedronTestMargin", 0.1);
        auto cutoutPath = lm::json::value<std::string>(prop, "cutoutpath");

        // scale
        scale_ = lm::json::value<lm::Float>(prop, "scale", 1.0);
        MODEL_SCALE = scale_;
        usepluecker_ = lm::json::value<bool>(prop, "usepluecker", false);
        auto pos = cutoutPath.find(".hdf5");
        cutoutPath = cutoutPath.substr(0, pos);
        
       //Scene * am = reinterpret_cast<ArepoMesh*> ( prop["arepoMesh_addr"].get<uintptr_t>() );

        //only need ArepoMesh implementation (Spectrum and TransferFunction aren't used) 
#ifdef MOCK_AREPO
        arepoMesh = std::make_unique<ArepoMeshMock>();
        arepoMeshRef = arepoMesh.get();
        auto arepoBound = arepoMesh->WorldBound();

        auto Ndp = arepoMesh->getNdp();
        auto DP = arepoMesh->getDP();
#else 
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

        arepoMesh = std::make_unique<ArepoMesh>(&tf);
        arepoMesh->ComputeVoronoiEdges();
        arepoMeshWrapper = std::make_unique<ArepoMeshWrapper>(arepoMesh.get());
        arepoMeshRef = arepoMeshWrapper.get();
        s = Spectrum::FromRGB(Config.rgbAbsorb);
        tf = TransferFunction(s);

         LM_INFO("constructed Volume Arepo");
         std::cout << " loaded snapshot " << std::endl;
        
        
        LM_INFO("constructed Volume Arepo");
        
        max_scalar_ = max_scalar();
        LM_INFO("max scalar  {}",max_scalar_);
        LM_INFO("mean scalar  {}", arepo->valBounds[TF_VAL_DENS*3 + 2] / MODEL_SCALE);
        LM_INFO("num gas  {}",NumGas);
        LM_INFO("test delaunay mesh");
        auto arepoBound = arepoMeshWrapper->WorldBound();
        int Ndp = arepoMeshWrapper->getNdp();
        point * DP = arepoMeshWrapper->getDP();
#endif
        bound_.max = lm::Vec3(arepoBound.pMin.x,arepoBound.pMin.y,arepoBound.pMin.z);
        bound_.min = lm::Vec3(arepoBound.pMax.x,arepoBound.pMax.y,arepoBound.pMax.z);
        //bound
        for (int i = 0; i < Ndp; i++) {
            auto currentP = lm::Vec3(DP[i].x,DP[i].y,DP[i].z);
            bound_.max = glm::max(bound_.max,currentP);
            bound_.min = glm::min(bound_.min,currentP);
        }
        bound_.min *= MODEL_SCALE;
        bound_.max *= MODEL_SCALE;
        LM_INFO("spatial bounds {},{},{} ; {},{},{}",bound_.min.x,bound_.min.y,bound_.min.z, bound_.max.x, bound_.max.y, bound_.max.z);

       

        //test if it is a delaunay mesh, i.e. there are no other points within a tetrahedron
   /*     for(size_t tetrai = 0; tetrai < arepoMesh->Ndt; tetrai++) {
        
            int a = arepoMesh->DT[tetrai].p[0];
            int b = arepoMesh->DT[tetrai].p[1];
            int c = arepoMesh->DT[tetrai].p[2];
            int d = arepoMesh->DT[tetrai].p[3];
            
            auto v0 = lm::Vec3(arepoMesh->DP[a].x, arepoMesh->DP[a].y, arepoMesh->DP[a].z);
            auto v1 = lm::Vec3(arepoMesh->DP[b].x, arepoMesh->DP[b].y, arepoMesh->DP[b].z);
            auto v2 = lm::Vec3(arepoMesh->DP[c].x, arepoMesh->DP[c].y, arepoMesh->DP[c].z);
            auto v3 = lm::Vec3(arepoMesh->DP[d].x, arepoMesh->DP[d].y, arepoMesh->DP[d].z);
            auto tetraVs = glm::tmat4x3<lm::Float>(v0,v1,v2,v3);
            auto ctr = lm::Vec3(arepoMesh->DTC[tetrai].cx,arepoMesh->DTC[tetrai].cy,arepoMesh->DTC[tetrai].cz);




            for(size_t tetraj = 0; tetraj < arepoMesh->Ndt; tetraj++) {
                glm::tmat4x3<lm::Float> tmpPVs;

                lm::Vec4 p0Dets;
                auto p0 = lm::Vec3(arepoMesh->DP[0].x, arepoMesh->DP[0].y, arepoMesh->DP[0].z);
                connectP(tetraVs,p0,tmpPVs);
                computeDeterminants(tmpPVs,p0Dets);

                lm::Vec4 p1Dets;
                auto p1 = lm::Vec3(arepoMesh->DP[1].x, arepoMesh->DP[1].y, arepoMesh->DP[1].z);
                connectP(tetraVs,p1,tmpPVs);
                computeDeterminants(tmpPVs,p1Dets);

                lm::Vec4 p2Dets;
                auto p2 = lm::Vec3(arepoMesh->DP[2].x, arepoMesh->DP[2].y, arepoMesh->DP[2].z);
                connectP(tetraVs,p2,tmpPVs);
                computeDeterminants(tmpPVs,p2Dets);

                lm::Vec4 p3Dets;
                auto p3 = lm::Vec3(arepoMesh->DP[3].x, arepoMesh->DP[3].y, arepoMesh->DP[3].z);
                connectP(tetraVs,p3,tmpPVs);
                computeDeterminants(tmpPVs,p3Dets);

                auto sameSign0 = 
                p0Dets[0] * p1Dets[0] >= 0.0f && 
                 p0Dets[0] * p2Dets[0] >= 0.0f &&
                  p0Dets[0] * p3Dets[0] >= 0.0f;
                  
                auto sameSign1 = 
                p0Dets[1] * p1Dets[1] >= 0.0f && 
                 p0Dets[1] * p2Dets[1] >= 0.0f &&
                  p0Dets[1] * p3Dets[1] >= 0.0f;
                  
                auto sameSign2 = 
                p0Dets[2] * p1Dets[2] >= 0.0f && 
                 p0Dets[2] * p2Dets[2] >= 0.0f &&
                  p0Dets[2] * p3Dets[2] >= 0.0f;
                  
                auto sameSign3 = 
                p0Dets[3] * p1Dets[3] >= 0.0f && 
                 p0Dets[3] * p2Dets[3] >= 0.0f &&
                  p0Dets[3] * p3Dets[3] >= 0.0f;
                  
                if (!sameSign0 || !sameSign1 || !sameSign2 || ! sameSign3) {
                    LM_ERROR("ERROR: overlapping tetrahedra" );
                }

                for(int i =0; i < 4; i++) {
                    int index = arepoMesh->DT[tetraj].p[i];
                    if(index != a && index != b && index != c && index != d) { //not one of the points to check
                        auto p = lm::Vec3(arepoMesh->DP[index].x, arepoMesh->DP[index].y, arepoMesh->DP[index].z);
                        if( glm::distance(ctr , p) < glm::distance(ctr,v0) ) {
                            LM_ERROR("ERROR: NOT A DELAUNAY MESH {} < {}",glm::distance(ctr , p), glm::distance(ctr,v0) );
                        }
                        //if(insideTetra(tetrai,p,Volume_Arepo_Impl::cachedSample())) //require that no points are within the tetrahedron abcd!
                        //    LM_ERROR("ERROR: NOT A DELAUNAY MESH");

                    }
                }
            }
            
        }*/
        


        meshAdapter = lm::comp::create<lm::ArepoLMMesh>( //lm::load<lm::Mesh>( 
        "mesh::arepo", make_loc("tetramesh"), {
            //{"ps_addr", (const point*)arepoMesh->DP},
           // {"ps_addr", reinterpret_cast<uintptr_t>(arepoMesh->DP)},
            //{"ps_count", arepoMesh->Ndp},
            //{"ts_addr", (const tetra*)arepoMesh->DT},
            //{"ts_addr", reinterpret_cast<uintptr_t>(arepoMesh->DT)},
            {"arepoMesh_addr", reinterpret_cast<uintptr_t>(arepoMesh.get())},
            {"scale", scale_}
            //{"ts_count", arepoMesh->Ndt}
        });
        
        auto & cached = Volume_Arepo_Impl::cachedDistanceSample();
        cached.sampleIndex = std::numeric_limits<long long>::max();
        cached.hydroI = 0;
        cached.looksAtTriId = -1;
        cached.tetraI = 0;
        
        dummyMat = lm::comp::create<lm::Material>(
        "material::diffuse", make_loc("dummymat"), {
            {"Kd", lm::Vec3(0.0f)}
        });

        
        //tetraaccel = lm::comp::create<lm::Accel>("accel::nanort", make_loc("tetraaccel"), {});
        tetraaccel = lm::comp::create<lm::Accel>("accel::embree", make_loc("tetraaccel"), {});
        LM_INFO( tetraaccel->loc());
        gastetrascene = lm::comp::create<lm::Scene>("scene::default", make_loc("tetrascene"), {
            {"accel", tetraaccel->loc()}
        });
        
        accelRef = tetraaccel.get();
        //
        auto c = gastetrascene->create_primitive_node({
             {"mesh" , meshAdapter->loc()},
             {"material" , dummyMat->loc()}
        });
        gastetrascene->add_child(gastetrascene->root_node(), c);
        
        gastetrascene->build();
        //accel->build(*gastetrascene.get());


        lighthierarchy = lm::comp::create<lm::AccelKnn>("accel::embree_knn", make_loc("lighthierarchy"), {});
        //LM_INFO( accel->loc());
        lightscene = lm::comp::create<lm::Scene>("scene::default", make_loc("lightscene"), {
            {"accel", lighthierarchy->loc()}
        });
        //TODO fill light hierachy with star positions
        //accelRef = accel.get();
        //
        //auto c = lightscene->create_primitive_node({
        //     {"mesh" , meshAdapter->loc()},
        //     {"material" , dummyMat->loc()}
        //});
        //lightscene->add_child(lightscene->root_node(), c);

       // lightscene->build();


        
    }

    virtual lm::Bound bound() const override {
        return bound_;        
    }
    virtual lm::Float max_scalar() const override {
#ifdef MOCK_AREPO
        return arepoMesh->max_density();
#else
        return arepo->valBounds[TF_VAL_DENS*3 + 1] / MODEL_SCALE;
#endif
    }
    virtual lm::Float max_scalar(lm::Ray ray, lm::Float & out_t_forhowlong, lm::Float & out_a, lm::Float & out_b) const override {

        auto & info = cachedDistanceSample();
        bool inside = findAndCacheTetra(info, ray.o, ray.d, meshAdapter.get()); 
        lm::Float maxScalarInCurrentTetra = inside ? glm::max(info.cornerVals[0][TF_VAL_DENS],
        glm::max(info.cornerVals[1][TF_VAL_DENS],
        glm::max(info.cornerVals[2][TF_VAL_DENS],
                    info.cornerVals[3][TF_VAL_DENS]))) : 0.0;
        auto nextBound = inside ? intersectCachedTetra(ray,info,usepluecker_) : info.lastHit.t;
        out_t_forhowlong = nextBound;
        if(inside) {
            sampleCachedScalarCoefficients(ray, out_a, out_b, info);
        } else {
            out_a = 0.0;
            out_b = 0.0;
        }
        return glm::max(maxScalarInCurrentTetra , 1.0 / nextBound);
    }


    virtual bool has_scalar() const override {
        return true;
    }



    void travel(lm::Ray ray,  CachedSample & useCache,  std::function<bool(bool inside,lm::Ray currentRay, lm::Float t, CachedSample & info)> processor) const {
        bool inside = false;
        lm::Float t = 0.0;
        auto & info = useCache; 
        auto correctedRay = ray;
        ray.o += ray.d *  tetraTestMargin; //perturbate a bit 
        do {
            //step
            correctedRay.o = ray.o;
            ray.o += ray.d * (t + tetraTestMargin);
            /*if(inside) { //last iteration we were inside, can test neighbors
                auto ni = arepoMeshRef->DT[info.tetraI].t[(info.looksAtTriId - 1) % 4]; //next tetrahedron with common face
                if(ni >= 0 && ni < arepoMeshRef->Ndt) {
                    auto oppositeNPoint = arepoMeshRef->DT[info.tetraI].s[(info.looksAtTriId - 1) % 4]; //next tetrahedron's opposite point (not part of the common face)
                    auto tet = arepoMeshRef->DT[ni];
                    inside = insideTetra(ni, tet,ray.o,info);
                
                    if(!inside && tet.t[0] >= 0) { //not the "simple" case where the ray resides in the next tetra... what to do?
                        //check next after next? 
                        for(int i = 1; i < 4 && ! inside; i++) {
                            auto altI = (oppositeNPoint + i) % 4;
                            auto nni = tet.t[altI];
                            inside = insideTetra(nni,ray.o,info);
                            //TODO maybe make intersection tests with deleted tetras!?
                        }
                    }
                }

            } */
            //if(!inside) {
                inside = findAndCacheTetra(info,ray.o, ray.d, meshAdapter.get());
            //}

            //this seems to be numerically unstable (depending on the vertex order of the triangle O.o)
            //try to hide it by interpolating with barycentric coordinates 3 times (with all vertex combis)
            //and average the result
            auto accelT = info.lastHit.t;

            
            if(info.looksAtTriId >= 0 && !std::isinf(accelT) && usepluecker_) {
                std::vector<lm::Vec3> verts = {info.tetraVs[0],info.tetraVs[1],info.tetraVs[2],info.tetraVs[3]};
                int enterface, leaveface;
                lm::Vec3 enterpoint, leavepoint;
                lm::Float uenter1, uenter2, uleave1, uleave2, tenter, tleave;
                bool retval = RayTetraPluecker(ray.o, ray.d,
                    verts.data(),
                    enterface, leaveface,
                    enterpoint,leavepoint,
                    uenter1, uenter2, uleave1,uleave2,
                    tenter, tleave);
                   // LM_INFO(" {} t enter and leave {}, {} vs {}",retval,tenter,tleave, accelT  );
                if(retval) {
                    if(tleave < tenter) { 
                        if(tleave < 0.0 && tenter > 0.0)  // tenter is valid
                            accelT = tenter;
                        if(tleave > 0.0) // tleave is valud
                            accelT = tleave;
                    }
                }
                  //  accelT = tenter < 0.0 ? tleave : tenter;
                
            }
            //LM_INFO("hit {} on fac {}",info.lastHit.primitive, info.lastHit.face );
            //auto hitPoint = meshAdapter->surface_point(info.lastHit.face, info.lastHit.uv);
            //auto accelT = glm::distance(hitPoint.p,ray.o);

            t = (!inside ? accelT : 
                intersectCachedTetra(ray,info,usepluecker_));
            //t = info.lastHit.t;
            correctedRay.o = ray.o;//+= ray.d * t;
            //LM_INFO("inside: {}, t: {}",inside,t);
        } while(processor(inside,correctedRay,t + tetraTestMargin ,info));

    }

    virtual lm::Float sample_distance(lm::Ray originalRay,lm::Float tmin, lm::Float tmax, lm::Rng& rng, lm::Float & weight) const override {

        //retrieve distance sample of strategy 0 (equiangular sampling) 
        int key = 0;
        auto equiAngularT = lm::stats::get<lm::stats::EquiangularStrategyDistanceSample,int,lm::Float>(key);
        //store pdf for sample of strategy 0
        //auto & pdfs = lm::stats::get<lm::stats::DistanceSamplesPDFs,lm::stats::IJ::,lm::Float>(key);
        lm::stats::IJ pdfkey = lm::stats::IJ::_0_0;
        //auto equipdf = lm::stats::get<lm::stats::DistanceSamplesPDFs,lm::stats::IJ,lm::Float>(pdfkey);
                             

        int currentRegularRngU_I = lm::stats::get<lm::stats::RegularDistanceSampleRandomValueVertexIndex,int,int>(key);
        //receive current random sample
        auto rng_u_regular = lm::stats::get<lm::stats::DistanceSampleRandomValues,int,lm::Float>(currentRegularRngU_I);
        
        auto totalRayVisitor = lm::stats::get<lm::stats::BoundaryVisitor,int,std::function<void(lm::Vec3,lm::RaySegmentCDF const &, int tetraindex)>>(key);

        auto visitor = lm::stats::get<lm::stats::BoundaryVisitor,int,std::function<void(lm::Vec3,lm::RaySegmentCDF const &)>>(key);

        auto visitor2 = lm::stats::get<lm::stats::BoundaryVisitor,int,std::function<void(lm::Vec3,lm::RaySegmentCDF const &, int tetraindex,lm::Float,lm::Float)>>(key);

        int segmentCount = 0;
        std::vector<lm::RaySegmentCDF> & segments = raySegments();
        lm::stats::set<lm::stats::LastBoundarySequence,int,std::vector<lm::RaySegmentCDF>*>(0,&segments);
        

        lm::Float totalEffectiveT = 0.0;
        lm::Float totalacc = 0.0;
        lm::Float totalT = tmin;
        lm::Float minT = tmin;
        bool setMinT = false;
        {
            auto & cached = cachedDistanceSample(); 
            auto travelray = originalRay;
            travelray.o = travelray.o + travelray.d * tmin; 
            lm::Float traveledT;

            if(segments.size() == 0)
                segments.resize(RAY_SEGMENT_ALLOC);
                segments[0].localcdf = 0.0;//TODO VALID!?
                segments[0].t = tmin;
                segments[0].a = 0.0;
                segments[0].b = 0.0;
                segments[0].tetraI = cached.tetraI; //smells ugly, but is used for performance reasons only (see renderer_volpt_arepo, l. 1281), does not affect correctness
                segmentCount++;

                if(totalRayVisitor)
                    totalRayVisitor(originalRay.o , segments[0],segments[0].tetraI );

            
            
            travel(travelray, cached,
            [&] (bool inside,lm::Ray currentRay, lm::Float t, CachedSample & info) -> bool {
                lm::RaySegmentCDF toFill;
                if(segments.size() <= segmentCount)
                    segments.resize(segments.size() * 2);

                bool ret = false;
                if(!inside) {
                    if(std::isinf(t)) { // we wont hit any tetra again
                        toFill.localcdf = 0.0;
                        toFill.t = tmax;
                        toFill.a = 0;
                        toFill.b = 0;
                        toFill.tetraI = info.tetraI; //last valid tetra
                    } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                        toFill.localcdf = 0.0;
                        toFill.t = t;
                        toFill.a = 0;
                        toFill.b = 0;
                        toFill.tetraI = info.tetraI; //upcoming valid tetra
                        ret = true;
                        totalT += toFill.t;
                        
                        if(!setMinT) {
                            setMinT = true;
                            minT = totalT;
                        }
                    }
                } else {
                    lm::Float a,b;
                    sampleCachedScalarCoefficients(currentRay,  a, b, info);
                    toFill.localcdf = sampleCDF( t, a, b);
                    toFill.t = t;
                    toFill.a = a;
                    toFill.b = b;
                    toFill.tetraI = info.tetraI;
                    //only count segments that will participate in scatter.
                    if(toFill.localcdf > std::numeric_limits<lm::Float>::epsilon())
                        totalEffectiveT += toFill.t;

                    totalacc += toFill.localcdf;
                    ret = true;
                    totalT += toFill.t;
                }
                if(totalRayVisitor)
                    totalRayVisitor(currentRay.o , toFill,toFill.tetraI );

                segments[segmentCount] = std::move(toFill);
                segmentCount++;
                return ret;
            });
        }
        lm::stats::set<lm::stats::RegularTrackingStrategyTotalT,int,lm::Float>(0,totalT);
        lm::stats::set<lm::stats::RegularTrackingStrategyTotalEffT,int,lm::Float>(0,totalEffectiveT);
        lm::stats::set<lm::stats::RegularTrackingStrategyMinT,int,lm::Float>(0,minT);

       // LM_INFO("total {},",totalacc);
        lm::Float normFac =  1.0 - glm::exp(-totalacc);
        lm::stats::set<lm::stats::RegularTrackingStrategyNormFac,int,lm::Float>(0,normFac);





        lm::Float integratedNormFac = 1.0;//1.0 - glm::exp(- (totalacc));
        //have found cdf normalization, now can sample the volume exactly following 
        //its density
        lm::Float transmittance = 1.0;
        lm::Float xi =  rng.u() * normFac;

        lm::stats::set<lm::stats::RegularTrackingStrategyXi,int,lm::Float>(0,xi);
        lm::stats::set<lm::stats::RegularTrackingStrategyTotalTau,int,lm::Float>(0,totalacc);

        lm::Float zeta = xi * rng.u();//within the interval, 
        //warp a random uniform variable across the transmittance cdf

        //lm::Float xi =  rng_u_regular * normFac;
        
        //lm::Float logxi = -glm::log( 1.0 - xi);
        lm::Float logxi = -gsl_log1p(-xi);
        lm::Float logzeta = -gsl_log1p(-zeta);
        lm::Float acc = 0.0;
        lm::Vec3 contribution = lm::Vec3(0.0);
        lm::Float pdf = 0.0;
        lm::Float regularPF_of_equiSample= 0.0;
        auto freeT = 0.0;
        auto retFreeT = std::numeric_limits<lm::Float>::max();
        auto retAcc = 0.0;
        bool stopAccumulating = false;
        bool stopZeta = false;
        bool stopEquiangular = false;

        lm::stats::set<lm::stats::LastBoundarySequence,int,int>(key,segmentCount);
        
        for(int segmentI = 0; segmentI < segmentCount && !stopAccumulating; segmentI++) {
            auto & raySegment = segments[segmentI];


            if(visitor) visitor(originalRay.o + originalRay.d * freeT, raySegment);
            //by the way, calculate PDF for samples coming from other strategies:
            /*TODO works only if iterating over the complete ray, taken out for optimization 
            if(freeT + raySegment.t > equiAngularT && !stopEquiangular) {
                stopEquiangular = true;
                auto t = (equiAngularT - freeT);
                auto cdf = acc + sampleCDF(t,raySegment.a,raySegment.b );
                auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                auto particle_density = raySegment.b + raySegment.a * t;
                auto mu_a = crosssection * particle_density;
                auto phase_integrated = 1.0;//isotrope
                auto mu_s = phase_integrated* particle_density;
                auto mu_t = mu_a + mu_s;
                auto scattering_albedo = mu_s / mu_t; 
                regularPF_of_equiSample = mu_t * glm::exp(-cdf) ;/// normFac; 
                

            }*/
            if ( !stopAccumulating &&  (acc  + raySegment.localcdf ) > logxi) {
                auto normcdf =  acc ;
                lm::Float t = sampleCachedICDF_andCDF( logxi,xi , raySegment.t ,
                normcdf ,  raySegment.a ,   raySegment.b );
               // retFreeT = tmin + freeT + t; //fixed this problem by starting at 0 not at tmin:
                retFreeT =  freeT + t; 
                normcdf = sampleCDF(t,raySegment.a,raySegment.b );
                transmittance *= glm::exp(-  normcdf );
                //accumulate non normalized!

        

                retAcc = acc + normcdf;
                
                auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                auto particle_density = raySegment.b + raySegment.a * t;
                auto mu_a = crosssection * particle_density;
                auto phase_integrated = 1.0;//isotropic
                auto mu_s = phase_integrated* particle_density;
                auto mu_t =  mu_a + mu_s;
                contribution = lm::Vec3(
                    A_R_A_V_S * mu_t * glm::exp(-A_R_A_V_T*retAcc),
                    mu_t * glm::exp(-retAcc),
                    A_B_A_V_S * mu_t * glm::exp(-A_B_A_V_T*retAcc)
                );

                //auto scattering_albedo = mu_s / mu_t; 
                //contribution = mu_s * transmittance;

                pdf = mu_t * transmittance / normFac;             
                //pdf = (raySegment.b + raySegment.a * t) * glm::exp(-retAcc) ; //normfac already at scatter albedo

                
                lm::stats::set<lm::stats::RegularTrackingStrategyTetraIndex,int,int>(0, raySegment.tetraI);
        
                stopAccumulating = true;
            }

            if (!stopAccumulating) {
                //pdf *= (raySegment.b + raySegment.a * raySegment.t) * glm::exp(-raySegment.localcdf) ;
                transmittance *= glm::exp(- raySegment.localcdf);
            }
            
            //visitor that needs to know the tetrahedron index corresponding to the ray segment as well
            if(visitor2)
                visitor2(originalRay.o + originalRay.d * freeT, raySegment,raySegment.tetraI,
                stopAccumulating ? retFreeT : freeT , totalT);

           
            acc += raySegment.localcdf ;
            freeT += raySegment.t;
            
        }
        

        //lm::stats::set<lm::stats::DistanceSamplesPDFs,lm::stats::IJ,lm::Float>(
        //            lm::stats::IJ::_1_0,regularPF_of_equiSample);
                
                
        lm::stats::set<lm::stats::ScatteringAlbedo,int,lm::Vec3>(0, contribution);


        //acccdf contains NON-normalized cdf which is exactly what we want
        //acccdf *= cdfNorm;
        lm::stats::set<lm::stats::FreePathTransmittance,int,lm::Float>(0,transmittance );

        lm::stats::set<lm::stats::RegularTrackingStrategyTauUntilScatter,int,lm::Float>(0,retAcc);

        //TODO multiply with normFAC?!?!
        lm::stats::set<lm::stats::DistanceSamplesPDFs,lm::stats::IJ,lm::Float>(lm::stats::IJ::_1_1, pdf);
        lm::stats::set<lm::stats::RegularTrackingStrategyDistanceSample,int,lm::Float>(key,retFreeT);
        //weight =  1.0/pdf; 
        
        return retFreeT;

   }

    virtual lm::Float eval_transmittance(lm::Ray originalRay, lm::Float tmin, lm::Float tmax) const override {
         
        //auto cached = cachedDistanceSample(); //work directly on cached, as the next request will start from this request's result
        auto cached = cachedTransmittanceSample();
        //check if we have information from sample_distance available
        int h = 0;        
        auto currentSample = lm::stats::get<lm::stats::CachedSampleId,int,long long>(h);
        if(cached.sampleIndex == currentSample) { //is cached
            //auto maxTransmittance = lm::stats::get<lm::stats::MaxTransmittance,int,lm::Float>(h);
            auto freePathTransmittance = lm::stats::get<lm::stats::FreePathTransmittance,int,lm::Float>(h);
          //  return freePathTransmittance;
        } 

        auto transmittance = 1.0;
        auto accT = tmin;
        auto accCDF = 0.0;
        originalRay.o = originalRay.o + originalRay.d * tmin; 
        lm::Float negligibleTransmittance = +0.0;
        travel(originalRay, cached,
        [&] (bool inside,lm::Ray currentRay, lm::Float nextT, CachedSample & info) -> bool {
            bool ret = false;
            lm::Float a,b;
            sampleCachedScalarCoefficients(currentRay,  a, b, info);

            if(!inside) {
                if(std::isinf(nextT)) { // we wont hit any tetra again
                } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                    accT += nextT;
                    accCDF += sampleCDF(nextT,a,b);
                    ret = true;
                }
            } else { 
                auto transmittanceT = glm::min(tmax - accT,  nextT);
                accCDF += sampleCDF(transmittanceT,a,b);
                transmittance *= sampleTransmittance(currentRay, 0.0,0.0 + transmittanceT  ,a,b, info);
                if( nextT + accT < tmax && transmittance >= negligibleTransmittance)
                    ret = true;
                accT += nextT;
            }
            return ret;
        });
        lm::stats::set<lm::stats::OpticalThickness,int,lm::Float>(0,accCDF);
      
        return transmittance;

    }
    


    virtual lm::Float eval_scalar(lm::Vec3 p) const override {
        
        thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals1;
        tmpVals1.clear();
        gatherValsAtPoint(p,glm::normalize(lm::Vec3(1)), tmpVals1.vals);


        //cache scattering albedo
        auto particle_density = tmpVals1.vals[TF_VAL_DENS];
        auto crosssection = 1.0;//327.0/1000.0; //barn ...? Hydrogen
        auto mu_a = crosssection * particle_density;
        auto phase_integrated = 1.0;//isotrope
        auto mu_s = phase_integrated * particle_density;
        auto mu_t = mu_a + mu_s;
        auto scattering_albedo = mu_s / mu_t; 
        lm::stats::set<lm::stats::ScatteringAlbedo,int,lm::Float>(0, scattering_albedo);

        return  tmpVals1.vals[TF_VAL_DENS];
    }

    virtual lm::Float eval_scalar(lm::Vec3 p , lm::Vec3 dir) const override {
        thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals1;
        tmpVals1.clear();
        gatherValsAtPoint(p, dir, tmpVals1.vals);

        //cache scattering albedo
        auto particle_density = tmpVals1.vals[TF_VAL_DENS];
        auto crosssection = 1.0;//327.0/1000.0; //barn ...? Hydrogen
        auto mu_a = crosssection * particle_density;
        auto phase_integrated = 1.0;//isotrope
        auto mu_s = phase_integrated* particle_density;
        auto mu_t = mu_a + mu_s;
        auto scattering_albedo = mu_s / mu_t; 
        lm::stats::set<lm::stats::ScatteringAlbedo,int,lm::Float>(0, scattering_albedo);


        return  tmpVals1.vals[TF_VAL_DENS];
    }

    virtual bool has_color() const override {
        return true;
    }
    virtual lm::Vec3 eval_color(lm::Vec3 p) const override {
        //thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals2;
        auto cached = cachedDistanceSample();
        //check if we have information from sample_distance available
        int h = 0;        
        auto currentSample = lm::stats::get<lm::stats::CachedSampleId,int,long long>(h);
        if(cached.sampleIndex == currentSample) { //is cached
            auto scatteringAlbedo = lm::stats::get<lm::stats::ScatteringAlbedo,int,lm::Float>(h);
            return lm::Vec3(scatteringAlbedo);
        } 

        //else search for it 
        if(! findAndCacheTetra(cached,p,lm::Vec3(1,0,0),meshAdapter.get())) {
            return lm::Vec3(1.0);
        } else {
            auto crosssection = 1.0;//327.0/1000.0; //barn ...?

            lm::Ray r;
            r.o = p;
            r.d = lm::Vec3(1,0,0);
            lm::Float a,b;
            sampleCachedScalarCoefficients(r,a,b,cached);
            auto particle_density = b;
            auto mu_a = crosssection * particle_density;
            auto phase_integrated = 1.0;//isotrope
            auto mu_s = phase_integrated* particle_density;
            auto mu_t = mu_a + mu_s;
            auto scattering_albedo = mu_s / mu_t; 
            return lm::Vec3(scattering_albedo);
        }


        //TODOOOOOOOOOOOOOOOOOOOOOOO
        /*auto copy = cachedSample();
        if(! findAndCacheTetra(cachedSample(),p)) {
            cachedSample() = copy ;
            return lm::Vec3(0.0); //wtf just add sth to land in tetra
        }
        cachedSample() = copy ;
        return lm::Vec3(1.0);
        int r = cachedSample().hydroI & 0x000000FF;
        int g = (cachedSample().hydroI & 0x0000FF00) >> 8;
        int b = (cachedSample().hydroI & 0x00FF0000) >> 16;

        tmpVals2.clear();
        //gatherValsAtPoint(p, tmpVals2.vals);
        return lm::Vec3(
            (lm::Float) r / 256.0,
            (lm::Float) g / 256.0,
            (lm::Float) b / 256.0); //some random stuff 
        //return Vec3(tmpVals2.vals[TF_VAL_BMAG],tmpVals2.vals[TF_VAL_METAL],tmpVals2.vals[TF_VAL_VMAG]); //some random stuff 
        */

    }

    //virtual void march(Ray ray, Float tmin, Float tmax, Float marchStep, const RaymarchFunc& raymarchFunc) const override {
      //  
    //}

    protected:
    lm::Float scale_;
    lm::Bound bound_;
    lm::Float max_scalar_;
    lm::Float tetraTestMargin;

#ifdef MOCK_AREPO
    std::unique_ptr<IArepoMeshMock>  arepoMesh;
#else
    std::unique_ptr<ArepoMesh>  arepoMesh;
    std::unique_ptr<IArepoMeshMock>  arepoMeshWrapper;
    
#endif
    std::unique_ptr<Arepo> arepo;

    Spectrum s;
    TransferFunction tf;

    private:    

    lm::Component::Ptr<lm::ArepoLMMesh> meshAdapter;
    lm::Component::Ptr<lm::Accel> tetraaccel;
    lm::Component::Ptr<lm::Scene> gastetrascene;


    lm::Component::Ptr<lm::AccelKnn> lighthierarchy;
    lm::Component::Ptr<lm::Scene> lightscene;

    lm::Component::Ptr<lm::Material> dummyMat;

    bool usepluecker_;



    static ArepoLoaderInternals::CachedSample & cachedDistanceSample() {
        thread_local ArepoLoaderInternals::CachedSample c;
        return c;
    }

    

    
    static ArepoLoaderInternals::CachedSample & cachedTransmittanceSample() {
        thread_local ArepoLoaderInternals::CachedSample cm;
        return cm;
    }

    

    void gatherValsAtPoint(lm::Vec3 p, lm::Vec3 dir, std::vector<float> & toVals) const {
        
        lm::Ray r;
        r.o = p;
        r.d = dir;
        auto & cached = cachedDistanceSample();
        if (! findAndCacheTetra(cached,r.o,r.d, meshAdapter.get())) {
           // LM_INFO("return nothing");
            toVals[TF_VAL_DENS] = 0.0; 

        }
        else {
            lm::Float a,b;
            sampleCachedScalarCoefficients(r,a,b,cached);
            //only need b
           // LM_INFO("return {}", b/invNorm);
            toVals[TF_VAL_DENS] = b;
        }

    }


};


LM_COMP_REG_IMPL(Volume_Arepo_Impl, "volume::arepo");


