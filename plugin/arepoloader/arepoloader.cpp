
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

#include <algorithm>
//extern ConfigSet Config;
//extern char ParameterFile[MAXLEN_PATH];
ConfigSet Config;



namespace ArepoLoaderInternals {

    #define INSIDE_TOLERANCE 10.0 * std::numeric_limits<lm::Float>::epsilon()
    #define TRAVEL_BIAS 0.01
    #define DENSITY_CRANKUP 100.0

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


            {
                tetra t;
                t.t[0] = 2;
                t.t[1] = 3;
                t.t[2] = 1;
                t.t[3] = 3;

                t.p[0] = 0;
                t.p[1] = 1;
                t.p[2] = 2;
                t.p[3] = 3;
                t.s[0] = 1;
                t.s[1] = 1;
                t.s[2] = 1;
                t.s[3] = 1;

                DT.push_back(t);
            }

            

            point p;
            p.x = -1;p.y = 0;p.z = -2;p.index = 0;
            DP.push_back(p);
            densities.push_back(0.5);
            p.x = 0;p.y = 0;p.z = 0;p.index=1;
            DP.push_back(p);
            densities.push_back(0.000000);
            p.x = 1;p.y = 0;p.z = -2;p.index=2;
            DP.push_back(p);
            densities.push_back(0.5);
            p.x = 0;p.y = 2;p.z = -1;p.index=3;
            DP.push_back(p);
            densities.push_back(0.000000);

            p.x = -1;p.y = 2;p.z = 0;p.index=4;
            DP.push_back(p);
            densities.push_back(0.000000);

            p.x = 1;p.y = 2;p.z = 0;p.index=5;
            DP.push_back(p);
            densities.push_back(0.0);


            {
                tetra t;
                t.t[0] = 3;
                t.t[1] = 3;
                t.t[2] = 3;
                t.t[3] = 0;

                t.p[0] = 0;
                t.p[1] = 1;
                t.p[2] = 3;
                t.p[3] = 4;
                t.s[0] = 1;
                t.s[1] = 1;
                t.s[2] = 1;
                t.s[3] = 1;

                DT.push_back(t);
            }

            {
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
            }


            tetra tDel;
            tDel.t[0] = -1;
            tDel.t[1] = -1;
            tDel.t[2] = -1;
            tDel.t[3] = -1;
            tDel.p[0] = 0;
            tDel.p[1] = 0;
            tDel.p[2] = 0;
            tDel.p[3] = 0;
            tDel.s[0] = 0;
            tDel.s[1] = 0;
            tDel.s[2] = 0;
            tDel.s[3] = 0;
            DT.push_back(tDel);

            Ndt = 4;
            Ndp = 6;


        }

        virtual lm::Float getDensity(int index) override  {
            return densities[index];
        }

        virtual lm::Float max_density() override  {
            return *std::max_element(densities.begin(), densities.end());
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



    inline double det3x3(lm::Vec3 b,lm::Vec3 c,lm::Vec3 d) {
        //glm::determinant(lm::Mat3(b0,b1,b2));
        return b[0]*c[1]*d[2] + c[0]*d[1]*b[2] + d[0]*b[1]*c[2] - d[0]*c[1]*b[2] - c[0]*b[1]*d[2] - b[0]*d[1]*c[2];
    }
    
    inline void connectP(glm::tmat4x3<lm::Float> & verts, lm::Vec3 & p, glm::tmat4x3<lm::Float> & pToVerts) {

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

    inline void sampleCachedCDFCoefficients(lm::Ray ray, lm::Float & a, lm::Float &  b, CachedSample const & cached, lm::Float & out_invNorm) {
        //this method assumes that the ray resides inside the current tetrahedron
        //therefore it clamps the barycentric coordinates in the b term to 1
        
        auto lambda012_a = cached.baryInvT * ray.d;
        auto lambda012_b = cached.baryInvT * (ray.o - cached.tetraVs[3]);

        auto maxd = glm::max(
            cached.cornerVals[0][TF_VAL_DENS] ,glm::max(
            cached.cornerVals[1][TF_VAL_DENS] ,glm::max(
            cached.cornerVals[2][TF_VAL_DENS] ,
            cached.cornerVals[3][TF_VAL_DENS] )));
        out_invNorm = 1.0 / maxd;
        auto densities = out_invNorm * lm::Vec4(
            cached.cornerVals[0][TF_VAL_DENS] ,
            cached.cornerVals[1][TF_VAL_DENS] ,
            cached.cornerVals[2][TF_VAL_DENS] ,
            cached.cornerVals[3][TF_VAL_DENS] );
        //densities *= scale_;
        b = glm::dot( lm::Vec4(lambda012_b, 1.0 - lambda012_b.x - lambda012_b.y - lambda012_b.z), densities);
        a = glm::dot( lm::Vec4(lambda012_a, 1.0 - lambda012_a.x - lambda012_a.y - lambda012_a.z), densities);
        
    }

    inline lm::Float sampleCDF(lm::Ray ray, lm::Float fromT, lm::Float toT, CachedSample const & cached) {
        lm::Float a,b, invNorm;
        sampleCachedCDFCoefficients(ray, a, b, cached,invNorm);
        return ( b * (toT-fromT) / invNorm +  a * 0.5 *(toT * toT - fromT * fromT) / invNorm );
    }

    inline lm::Float sampleTransmittance(lm::Ray ray, lm::Float fromT, lm::Float toT, CachedSample const & cached) {
        return glm::exp(- sampleCDF(ray,fromT,toT,cached));
    }



    inline lm::Float sampleCachedICDF_andCDF(lm::Ray ray, lm::Float logxi, lm::Float tmin, lm::Float tmax, lm::Float & out_cdf, CachedSample const & cached) {
        
        lm::Float a,b,invNorm;
        sampleCachedCDFCoefficients(ray, a, b, cached,invNorm);
        
        //use tau*_t1 (t) which is the integral from t1 to t minus the integral from 0 to t1
        auto y = logxi + a *  0.5 * tmin*tmin  / invNorm + b * tmin / invNorm - out_cdf;
        
        //lm::Float freeT = glm::sqrt(
        //    ((b*b)/(4.0*a*a)) + ( y / (a / invNorm )) ) - b/(2.0*a);
        lm::Float freeT;
        


        freeT = glm::sqrt(
                ((b*b)/(a*a)) + ( 2.0 * y / (a / invNorm )) ) - b/a;
        
        freeT = isnan(freeT) ? std::numeric_limits<lm::Float>::max() : freeT;
        //for evaluating tau, limit free path to tmax
        lm::Float t = glm::min(freeT + tmin, tmax);
        
            out_cdf += 
            ((t - tmin) * b / invNorm +  a  * 0.5 * (t * t - tmin * tmin) / invNorm) ; //cdf within tmin and  min of (t , tmax) 
        return freeT;//returns sth between tmin - tmin (so 0) and tmax - tmin 
    }

    inline lm::Float intersectCachedTetra(lm::Ray ray, CachedSample & cached)  {

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
        lm::Float bary0 = std::abs( det3x3(cached.tmpPVs[indices[1]],cached.tmpPVs[indices[2]], ray.d) / mainD);
        lm::Float bary1 = std::abs( det3x3(cached.tmpPVs[indices[0]],cached.tmpPVs[indices[2]], ray.d) / mainD);
        lm::Float bary2 = std::abs( det3x3(cached.tmpPVs[indices[0]],cached.tmpPVs[indices[1]], ray.d) / mainD);
        lm::Float baryOrigin = 1.0 - bary0 - bary1 - bary2;
        
        auto t = 1.0 / (bary0 + bary1 + bary2); 
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
        return t;
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
        * (1.0 / cachedS.mainDeterminant);

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
                    addValsContribution(cachedS.cornerVals[i],hydroIndex,DENSITY_CRANKUP);//lengths[minDistIndex] / totalD );
                }  
            }
#endif
        }
        
         
    }

    //performs point in tetrahedron test, returns result. 
    //stores all relevant test data into cachedS if p is in the tetra
    inline bool insideTetra(int tetraIndex, tetra tetra, lm::Vec3 p, CachedSample & cachedS)  {
        glm::tmat4x3<lm::Float> pVs;//point to vertex connections
        glm::tmat4x3<lm::Float> verts;
        glm::ivec4 vertInds;
        lm::Vec4 determinants;
        
        
        for(int i = 0; i < 4; i++) {
            vertInds[i] = tetra.p[i];
            auto av = arepoMeshRef->getDP()[vertInds[i]];
            verts[i] = lm::Vec3(av.x,av.y,av.z);
            pVs[i] = verts[i] - p;
        }

        //also transports sign (choose vertex 3 as "roof"): negative means ccw, positive cw tetrahedron definition
        lm::Float mainDeterminant = det3x3(verts[0] - verts[3],verts[1] - verts[3],verts[2] - verts[3]);
        //if(mainDeterminant > 0.0)
         //   LM_INFO( "cw");
        
        
        //skip points 
        //wtf see arepo vtk
        DPinfinity = -5;
        if  (
        (tetra.t[0] < 0 ||tetra.p[0] == DPinfinity || tetra.p[1] == DPinfinity
        || tetra.p[2] == DPinfinity || tetra.p[3] == DPinfinity)
        || tetra.t[0] == -1)
        {
            //LM_INFO("skip");
            return false;
        }

        

        connectP(verts,p,pVs);
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

    //overload
    inline bool insideTetra(int tetra, lm::Vec3 p, CachedSample & cachedS)  {
        return insideTetra(tetra, arepoMeshRef->getDT()[tetra], p, cachedS);
    }

    inline void updateCachedNeighbors(CachedSample & cached) {
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
        auto mostprobable = arepoMeshRef->getDT()[ofTetra].t[(cached.looksAtTriId - 1) % 4];
        //neighbors.push_back(mostprobable);
        //neighborInds.push_back(mostprobable.t[(cached.looksAtTriId - 1) % 3]);
        if(mostprobable >= 0 && mostprobable < arepoMeshRef->getNdt() && arepoMeshRef->getDT()[mostprobable].t[0] >= 0) {
            neighbors.push_back(arepoMeshRef->getDT()[mostprobable]);
            neighborInds.push_back(mostprobable);
            //addNeighbors(mostprobable,ofTetra);
        }
            
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
            if(foundNeighbor)
                return cached.neighborInds[i];
        }


        //still not found !?
        //then we are outside of the arepo mesh! 
        //if(!foundBoundary)
        //   LM_ERROR("tetra mesh traversal: didnt find neighbor but was not at boundary");
        return -1; 

    }

   

    inline bool findAndCacheTetra( CachedSample & cachedS, lm::Vec3 p, lm::Vec3 dir, lm::ArepoLMMesh * toQueryTetraId)  {
        int h = 0;
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
            else { //check neighbors
                updateCachedNeighbors(cachedS);
                int tetraIndex = checkCachedNeighbors(p,cachedS);
                bool inside = tetraIndex >= 0;
                if(inside) { //sample changed to a neighbor
                    lm::stats::add<lm::stats::UsedNeighborTetra,int,long long>(h,1);
                    //need to invalidate some information
                    cachedS.hydroI = -1;
                    cachedS.tetraI = tetraIndex;
                    cachedS.minDistI = -1;
                    //also store density at corners
                    cacheCornerValues(cachedS);
                    return true;
                } else {
                    //we totally lost the tetrahedron, need to make request to acceleration structure
                }
            }
        } else {
            cachedS.sampleIndex = std::numeric_limits<long long>::max();
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
            int localFaceIndex = (hit.value().face - 1) % 4;
            int tetraIndex =  toQueryTetraId->correspondingTetra(hit.value().face);
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
                updateCachedNeighbors(cachedS);
                //LM_INFO("found tetra inside");
                cachedS.sampleIndex = currentSample;
                //need to invalidate some information
                cachedS.hydroI = -1;
                cachedS.minDistI = -1;
               //LM_ERROR("this doesnt happend unfortunately?!");
                cacheCornerValues(cachedS);
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
            returnValue = false;
            cachedS.lastHit.t = std::numeric_limits<lm::Float>::infinity(); 
        }
        
        return returnValue;
    }


}

using namespace ArepoLoaderInternals;

class Volume_Arepo_Impl final : public lm::Volume_Arepo {

    
    public:
    
    Volume_Arepo_Impl() : arepoMesh(nullptr), arepo(nullptr), s(0.0f), tf(s),scene(nullptr) {

    }
    ~Volume_Arepo_Impl() {
        
        
    }

    virtual Component* underlying(const std::string& name) const override {
        if (name == "tetramesh") {
            return  dynamic_cast<lm::Mesh*>(meshAdapter.get());
            
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
        // Density scale
        scale_ = lm::json::value<lm::Float>(prop, "scale", 1.0);
        
        LM_INFO("constructed Volume Arepo");
        
        max_scalar_ = max_scalar();
        LM_INFO("max scalar  {}",max_scalar_);
        LM_INFO("mean scalar  {}",arepo->valBounds[TF_VAL_DENS*3 + 2]);
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

        
        accel = lm::comp::create<lm::Accel>("accel::embree", make_loc("tetraaccel"), {});
        LM_INFO( accel->loc());
        scene = lm::comp::create<lm::Scene>("scene::default", make_loc("tetrascene"), {
            {"accel", accel->loc()}
        });
        
        accelRef = accel.get();
        //
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
#ifdef MOCK_AREPO
        return arepoMesh->max_density();
#else
        return DENSITY_CRANKUP * arepo->valBounds[TF_VAL_DENS*3 + 1];
#endif
    }
    virtual bool has_scalar() const override {
        return true;
    }



    void travel(lm::Ray ray,  CachedSample & useCache, std::function<bool(bool inside,lm::Ray currentRay, lm::Float t, CachedSample & info)> processor) const {
        bool inside = false;
        lm::Float t = 0.0;
        auto & info = useCache; 
        auto correctedRay = ray;
        do {
            //step
            ray.o += ray.d * (t + TRAVEL_BIAS);
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
            t = (!inside ? info.lastHit.t : 
                intersectCachedTetra(ray,info));
            //t = info.lastHit.t;
            correctedRay.o = ray.o - ray.d * TRAVEL_BIAS;
            //LM_INFO("inside: {}, t: {}",inside,t);
        } while(processor(inside,correctedRay,t + TRAVEL_BIAS,info));

    }

    virtual lm::Float sample_distance(lm::Ray originalRay,lm::Float tmin, lm::Float tmax, lm::Rng& rng) const override {
        //auto xi = rng.u();
        bool sampleWasUpdated = false;
        auto cached = cachedDistanceSample(); //work directly on cached, as the next request will start from this request's result
        auto freeT = tmin;

        originalRay.o = originalRay.o + originalRay.d * tmin; 
        lm::Float inoutcdfValue = 0.0;
        //TODO: cdf ergibt nicht 1
        auto logxi = -glm::log(1.0-rng.u()) ;

        travel(originalRay, cached,
        [&] (bool inside,lm::Ray currentRay, lm::Float t, CachedSample & info) -> bool {
            bool ret = false;
            if(!inside) {
                if(std::isinf(t)) { // we wont hit any tetra again
                    //do nothing, return because the sample managed to travel through everything without scattering
                    freeT += tmax;//cached.lastHit.t;
                } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                    freeT += t;
                    ret = true;
                }
            } else {
                auto freeTCandidate = sampleCachedICDF_andCDF(currentRay, logxi, 0.0, 0.0 + t , inoutcdfValue, info);
                if(freeTCandidate > t) {  //we have to continue with the next tetra
                    ret = true;
                    freeT += t;
                    
                } else {
                    freeT += freeTCandidate; //we stop inside 
                    //freeT += t; //we stop inside 
                }
            }
            return ret;
        });
        //LM_INFO("return freeT {}", freeT);

        return freeT;

   }

    virtual lm::Float eval_transmittance(lm::Ray originalRay, lm::Float tmin, lm::Float tmax) const override {
        
        auto accT = tmin;
        originalRay.o = originalRay.o + originalRay.d * tmin; 
        auto transmittance = 1.0;
        lm::Float negligibleTransmittance = +0.0;
        auto cached = cachedTransmittanceSample(); //work directly on cached, as the next request will start from this request's result
        travel(originalRay, cached,
        [&] (bool inside,lm::Ray currentRay, lm::Float nextT, CachedSample & info) -> bool {
            bool ret = false;
            if(!inside) {
                if(std::isinf(nextT)) { // we wont hit any tetra again
                } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                    accT += nextT;
                    ret = true;
                }
            } else { 
                auto transmittanceT = glm::min(tmax - accT,  nextT);
                transmittance *= sampleTransmittance(currentRay, 0.0,0.0 + transmittanceT  , info);
                if(nextT  < tmax - accT && transmittance >= negligibleTransmittance)
                    ret = true;
                accT += nextT;//TODO ist das richtig?! wird hier nicht zu wenig addiert, wenns in nächster iteration auf jeden fall vom boundary weitergeht?
            }
            return ret;
        });

        return transmittance;



    }
    


    virtual lm::Float eval_scalar(lm::Vec3 p) const override {
        
        thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals1;
        tmpVals1.clear();
        if(naive)
            gatherValsAtPointNaive(p, tmpVals1.vals);
        else
            LM_ERROR("not naive but used eval_scalar method");
        //if(tmpVals1.vals[TF_VAL_DENS] > 0.0f)
        //   LM_INFO("{}", tmpVals1.vals[TF_VAL_DENS]);
        return  tmpVals1.vals[TF_VAL_DENS];
    }

    virtual bool has_color() const override {
        return true;
    }
    virtual lm::Vec3 eval_color(lm::Vec3 p) const override {
        //thread_local ArepoLoaderInternals::ArepoTempQuantities tmpVals2;


        return lm::Vec3(1.0);
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
    lm::Component::Ptr<lm::Accel> accel;
    lm::Component::Ptr<lm::Scene> scene;
    lm::Component::Ptr<lm::Material> dummyMat;
    bool naive;



    static ArepoLoaderInternals::CachedSample & cachedDistanceSample() {
        thread_local ArepoLoaderInternals::CachedSample c;
        return c;
    }

    
    static ArepoLoaderInternals::CachedSample & cachedTransmittanceSample() {
        thread_local ArepoLoaderInternals::CachedSample cm;
        return cm;
    }

    

    void gatherValsAtPointNaive(lm::Vec3 p, std::vector<float> & toVals) const {
        
        lm::Ray r;
        r.o = p;
        r.d = glm::normalize(lm::Vec3(1));
        auto & cached = cachedDistanceSample();
        if (! findAndCacheTetra(cached,r.o,r.d, meshAdapter.get())) {
           // LM_INFO("return nothing");
            toVals[TF_VAL_DENS] = 0.0; 

        }
        else {
            lm::Float a,b,invNorm;
            sampleCachedCDFCoefficients(r,a,b,cached,invNorm);
            //only need b
           // LM_INFO("return {}", b/invNorm);
            toVals[TF_VAL_DENS] = b/invNorm;
        }

    }


};


LM_COMP_REG_IMPL(Volume_Arepo_Impl, "volume::arepo");


