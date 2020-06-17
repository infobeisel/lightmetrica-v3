
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

    #define INSIDE_TOLERANCE 2.0 * std::numeric_limits<lm::Float>::epsilon()

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

    static ArepoMesh * arepoMeshRef = nullptr;
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



    double det3x3(lm::Vec3 b,lm::Vec3 c,lm::Vec3 d) {
        //glm::determinant(lm::Mat3(b0,b1,b2));
        return b[0]*c[1]*d[2] + c[0]*d[1]*b[2] + d[0]*b[1]*c[2] - d[0]*c[1]*b[2] - c[0]*b[1]*d[2] - b[0]*d[1]*c[2];
    }
    
    void connectP(glm::tmat4x3<lm::Float> & verts, lm::Vec3 & p, glm::tmat4x3<lm::Float> & pToVerts) {

        for(int i = 0; i < 4; i++)
            pToVerts[i] = verts[i] - p;
    }

    void computeDeterminants(glm::tmat4x3<lm::Float> & pToVerts, lm::Vec4 & results) {
        results[0] = det3x3(pToVerts[1],pToVerts[2],pToVerts[3]);
        results[1] = det3x3(pToVerts[0],pToVerts[2],pToVerts[3]);
        results[2] = det3x3(pToVerts[0],pToVerts[1],pToVerts[3]);
        results[3] = det3x3(pToVerts[0],pToVerts[1],pToVerts[2]); 
    }

    //this only works iff the tetrahedron is specified with a specific vertex order, namely, see mesh_arepo.cpp !!!
    bool insideOld(lm::Vec4 & determinants) {
        auto ret0 = determinants[0] > 0.0 && determinants[1] < 0.0 && determinants[2] > 0.0 && determinants[3] < 0.0;
        auto ret1 = determinants[0] < 0.0 && determinants[1] > 0.0 && determinants[2] < 0.0 && determinants[3] > 0.0;
        return ret0 || ret1;
    }

    //this test works irrespective of vertex order of the tetrahedron
    bool inside(lm::Vec4 & determinants, lm::Float tetraDeterminant) {
        return abs(abs(tetraDeterminant) - (abs(determinants[0])+abs(determinants[1])+abs(determinants[2])+abs(determinants[3]))) < INSIDE_TOLERANCE;       
    }


    bool insideCachedTetra(lm::Vec3 p, CachedSample & c) {
        connectP(c.tetraVs,p,c.tmpPVs);
        computeDeterminants(c.tmpPVs,c.tmpDets);
        return inside(c.tmpDets, c.mainDeterminant);
    }

    void sampleCachedCDFCoefficients(lm::Ray ray, lm::Float & a, lm::Float &  b, CachedSample const & cached, lm::Float & out_invNorm) {
        auto lambda012_a = cached.baryInvT * ray.d * 0.5;
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
        b = glm::dot( lm::Vec4(lambda012_b, 1.0 - lambda012_b.x - lambda012_b.y - lambda012_b.z),densities);
        a = glm::dot( lm::Vec4(lambda012_a, 1.0 - lambda012_a.x - lambda012_a.y - lambda012_a.z), densities);

    }

    lm::Float sampleCDF(lm::Ray ray, lm::Float fromT, lm::Float toT, CachedSample const & cached) {
        lm::Float a,b, invNorm;
        sampleCachedCDFCoefficients(ray, a, b, cached,invNorm);
        return ( b * (toT-fromT) / invNorm + a * (toT * toT - fromT * fromT) / invNorm );
    }

    lm::Float sampleTransmittance(lm::Ray ray, lm::Float fromT, lm::Float toT, CachedSample const & cached) {
        lm::Float a,b,invNorm;
        sampleCachedCDFCoefficients(ray, a, b, cached,invNorm);
        return glm::exp(- ( b * (toT-fromT) / invNorm + a * (toT * toT - fromT * fromT) / invNorm ));
    }



    lm::Float sampleCachedICDF_andCDF(lm::Ray ray, lm::Float logxi, lm::Float tmin, lm::Float tmax, lm::Float & out_cdf, CachedSample const & cached) {
        //now sample the distance, using the inverse cdf method.
        //our cdf is e to the power of the integral over the weighted sum of density values of the tetra vertices
        //the weights are the tetrahedral barycentric coordinates
        //so we have to extract the cdf and then invert it. the resulting inverse function can be cached partly 
        //as long we stay in the same tetrahedron!
        //e^tau(t) , tau is integral_0_t( sum_0_3( bary_i(s) * density_i ) ) ds
        //and t moves along ray ray(t) = o + d * t
        //bary_i(s) is the barycentric coordinate of point p on ray at s:
        //and according to wikipedia (doublecheck!) 
        // Tinv * (o-r4) + Tinv * d * t
        // with Tinv being a special matrix which we cache and r4 being the 4th tetra point
        // a order 2 polynomial, set a = Tinv * d * 0.5  and b = Tinv * (o - r4), c= 0
        //then the inverse is x = sqrt(b^2/4a^2 + y/a) - b/a
        // for y put in -ln(rnd) with rnd being sample value
        //if x exceeds distance within tetra (see further up), we need to continue with neighbor tetra. 
        
        lm::Float a,b,invNorm;
        sampleCachedCDFCoefficients(ray, a, b, cached,invNorm);
        
        //use tau*_t1 (t) which is the integral from t1 to t minus the integral from 0 to t1
        auto y = logxi + a * tmin*tmin  / invNorm + b * tmin / invNorm - out_cdf;
        
        lm::Float freeT = glm::sqrt(
            ((b*b)/(4.0*a*a)) + ( y / (a / invNorm )) ) - b/(2.0*a);
        freeT = isnan(freeT) ? std::numeric_limits<lm::Float>::max() : freeT;
        //for evaluating tau, limit free path to tmax
        auto t = glm::min(freeT + tmin, tmax) ;
        out_cdf += (t - tmin) * b / invNorm + a  * (t * t - tmin * tmin) / invNorm; //cdf within tmin and  min of (t , tmax) 
        
        return freeT;//returns sth between tmin - tmin (so 0) and tmax - tmin 
    }

    lm::Float intersectCachedTetra(lm::Ray ray, CachedSample & cached)  {

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
    int getCorrectedHydroInd(int dpInd)
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



    void updateCachedBaryInvT(CachedSample & cachedS)  {
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

    

    //performs point in tetrahedron test, returns result. 
    //stores all relevant test data into cachedS if p is in the tetra
    bool insideTetra(int tetraIndex, tetra tetra, lm::Vec3 p, CachedSample & cachedS)  {
        glm::tmat4x3<lm::Float> pVs;//point to vertex connections
        glm::tmat4x3<lm::Float> verts;
        glm::ivec4 vertInds;
        lm::Vec4 determinants;
        
        
        for(int i = 0; i < 4; i++) {
            vertInds[i] = tetra.p[i];
            auto av = arepoMeshRef->DP[vertInds[i]];
            verts[i] = lm::Vec3(av.x,av.y,av.z);
            pVs[i] = verts[i] - p;
        }

        //also transports sign (choose vertex 3 as "roof"): negative means ccw, positive cw tetrahedron definition
        lm::Float mainDeterminant = det3x3(verts[0] - verts[3],verts[1] - verts[3],verts[2] - verts[3]);
        if(mainDeterminant > 0.0)
            LM_INFO( "cw");
        
        
        //skip points 
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
        }
       
        return insideTet;
        
    }

    //overload
    bool insideTetra(int tetra, lm::Vec3 p, CachedSample & cachedS)  {
        return insideTetra(tetra, arepoMeshRef->DT[tetra], p, cachedS);
    }

    void updateCachedNeighbors(CachedSample & cached) {
        auto ofTetra = cached.tetraI;
        bool foundNeighbor = false;
        bool foundBoundary = false;
        //double check if in original tetra
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
            int * p1s = arepoMeshRef->DT[tet].p;
            int * p2s = arepoMeshRef->DT[of_tet].p;
            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++)
                    if(p1s[i] == p2s[j])
                        return true;
            return false;
        };

        auto isNeighbor2 = [&] (int tet, int of_tet) {
            bool foundN2 = isNeighbor(tet, of_tet);
            for(int i = 0;i < 4; i++) {
                int firstNeighbor = arepoMeshRef->DT[tet].t[i];
                if(firstNeighbor >= 0)
                    foundN2 = foundN2 || isNeighbor(firstNeighbor,tet);
            }
            return foundN2;
            
        };

        std::function<void(int,int)> addNeighbors = [&] (int of_tet, int original_tet) -> void {
            auto tetStrct = arepoMeshRef->DT[of_tet]; 
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
        auto mostprobable = arepoMeshRef->DT[ofTetra].t[(cached.looksAtTriId - 1) % 4];
        //neighbors.push_back(mostprobable);
        //neighborInds.push_back(mostprobable.t[(cached.looksAtTriId - 1) % 3]);
        if(mostprobable >= 0) {
            neighbors.push_back(arepoMeshRef->DT[mostprobable]);
            neighborInds.push_back(mostprobable);
            //addNeighbors(mostprobable,ofTetra);
        }
            
        //else
        //    addNeighbors(ofTetra,ofTetra);
        
    
    }

    //assumes the cached sample has been updated,
    //returns the index to the neighbor tetraeder the queried point is in. returns -1 if it is
    //outside the current and not in any neighbor tetraeder (then it should be outside the arepoMeshRef)
    int checkCachedNeighbors(lm::Vec3 p, CachedSample & cached)  {
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


    void evaluateDensityCached( std::vector<float> & toVals, CachedSample & cachedS)  {
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
            int num = arepoMeshRef->DP[cachedS.tetraInds[minDistIndex]].index;
            int hydroIndex = getCorrectedHydroInd(num);
            for (int i = 0; i < 9; i++)
                    cachedS.values[i] = 0.0f;
            if (hydroIndex > -1 && num  < NumGas ) {
                addValsContribution(cachedS.values,hydroIndex,1.0);//lengths[minDistIndex] / totalD );
                toVals[TF_VAL_DENS] = glm::max(0.0f,glm::min(1.0f,cachedS.values[TF_VAL_DENS]));
            } 
        }                
        cachedS.hydroI = hydroIndex;
        cachedS.minDistI = minDistIndex;
    }

    void cacheCornerValues(CachedSample & cachedS)  {

        for(int i = 0; i < 4; i++) {
            for (int j = 0; j < 9; j++)
                cachedS.cornerVals[i][j] = 0.0f;
            int num = arepoMeshRef->DP[cachedS.tetraInds[i]].index;
            if(num >= 0) {
                int hydroIndex = getCorrectedHydroInd(num);
                
                if(hydroIndex > -1 && num  < NumGas) {
                    addValsContribution(cachedS.cornerVals[i],hydroIndex,100.0);//lengths[minDistIndex] / totalD );
                }  
            }
        }
        
         
    }

    

    bool findAndCacheTetra( CachedSample & cachedS, lm::Vec3 p, lm::Vec3 dir, lm::ArepoLMMesh * toQueryTetraId)  {
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
            lm::stats::add<lm::stats::SampleIdCacheMisses,int,long long>(h,1);
            cachedS.sampleIndex = std::numeric_limits<long long>::max();
        }
        
        lm::stats::add<lm::stats::ResampleAccel,int,long long>(h,1);
        //uncached, need to ray intersect with volume
        lm::Ray r; 
        r.o = p;
        r.d = dir;
        auto hit = accelRef->intersect(r,0.0,std::numeric_limits<lm::Float>::max());

        if (hit != std::nullopt && hit.value().face >= 0 && hit.value().face < toQueryTetraId->num_triangles()) { //check if inside the tetra of hit triangle
            cachedS.lastHit = hit.value(); //save hit value
            int localFaceIndex = (hit.value().face - 1) % 4;
            int tetraIndex =  toQueryTetraId->correspondingTetra(hit.value().face);
            //auto ni = arepoMeshRef->DT[tetraIndex].t[localFaceIndex];
            bool inside = insideTetra(tetraIndex, p, cachedS) ;//|| insideTetra(ni, p, cachedS);
            //can be outside because we have double triangles in mesh, 
            //having exact same positions but belonging to two opposing tetrahedra,
            //need to check for both.
            //check neighbors
           // if(!inside) {
            //    LM_INFO("test shit?");
            //}
           for(int i = 0; i < 4; i++) {
                auto ni = arepoMeshRef->DT[tetraIndex].t[i];
                if(ni >= 0) {
                    inside = inside ||  insideTetra(ni, p, cachedS);
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
        arepoMeshRef = arepoMesh.get();
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
        max_scalar_ = max_scalar();
        LM_INFO("max scalar  {}",max_scalar_);
        LM_INFO("mean scalar  {}",arepo->valBounds[TF_VAL_DENS*3 + 2]);
        LM_INFO("num gas  {}",NumGas);
            LM_INFO("test delaunay mesh");

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

        
        Volume_Arepo_Impl::cachedDistanceSample().sampleIndex = std::numeric_limits<long long>::max();
        
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
        return arepo->valBounds[TF_VAL_DENS*3 + 1];
    }
    virtual bool has_scalar() const override {
        return true;
    }

    virtual lm::Float sample_distance(lm::Ray originalRay,lm::Float tmin, lm::Float tmax, lm::Float xi) const override {
        //return 1000.0;
        bool sampleWasUpdated = false;
        auto cached = cachedDistanceSample(); //work directly on cached, as the next request will start from this request's result
        auto freeT = tmin;

        auto probeRay = originalRay;
        auto probeT = tmin;
        probeRay.o = originalRay.o + probeRay.d * probeT; 
        bool found = false;
        lm::Float inoutcdfValue = 0.0;
        auto logxi = -glm::log(1.0-xi);

        while(!found && freeT < tmax) {
            if(! findAndCacheTetra(cached,probeRay.o,probeRay.d, meshAdapter.get())) {
                //not inside, there has been performed an up to date ray intersection to the volume bound and we can use it
                if(std::isinf(cached.lastHit.t)) { // we wont hit any tetra again
                    //do nothing, return because the sample managed to travel through everything without scattering
                    freeT += tmax;//cached.lastHit.t;
                } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                    //x (our cdf value) stays the same, move probe until the next tetra
                    probeRay.o = probeRay.o + probeRay.d * (cached.lastHit.t + 0.01);
                    freeT += cached.lastHit.t;
                }
            } else {
                //otherwise the last Hit information is outdated, we are in a tetrahedron and 
                //have to accumulate the free path for this tetrahedron
                
            
            //TODO switch back!
                //auto t2_fromt1test = intersectCachedTetra(probeRay,cached) ; // because we are already inside the tetra by 0.01
                auto t2_fromt1 = intersectCachedTetra(probeRay,cached) ; // because we are already inside the tetra by 0.01
                //auto t2_fromt1 = cached.lastHit.t; // because we are already inside the tetra by 0.01
               // LM_INFO("distance vs truth : {} vs {}", t2_fromt1test ,t2_fromt1);
      //returns unlimited free path, but accumulates transmittance only for inside the tetrahedron
                auto freeTCandidate = sampleCachedICDF_andCDF(probeRay, logxi, 0.0, t2_fromt1 , inoutcdfValue, cached);
                if(freeTCandidate > t2_fromt1) {  //we have to continue with the next tetra
                    probeRay.o = probeRay.o + probeRay.d * (t2_fromt1 + 0.01); //be sure to land in next tetra
                    freeT += t2_fromt1;
                } else {
                    freeT += freeTCandidate; //we stop inside 
                    found = true;
                }
            }
        }
        return freeT;
   }

    virtual lm::Float eval_transmittance(lm::Ray originalRay, lm::Float tmin, lm::Float tmax) const override {
        //find the starting tetra if existent
        auto & cached = cachedTransmittanceSample();
        auto t1 = tmin;
        auto t = tmin;
        auto probeRay = originalRay;
        auto probeT = t1;
        probeRay.o = originalRay.o + probeRay.d * probeT; 
        auto transmittance = 1.0;
        lm::Float negligibleTransmittance = +0.0;
        while(t <= tmax && transmittance >= negligibleTransmittance) {
            if(! findAndCacheTetra(cached,probeRay.o,probeRay.d, meshAdapter.get())) {
                //not inside, there has been performed an up to date ray intersection to the volume bound and we can use it
                if(std::isinf(cached.lastHit.t)) { // we wont hit any tetra again
                    //do nothing, return because the sample managed to travel through everything without scattering
                    return transmittance;//cached.lastHit.t;
                    
                } else { //currently we aren't in any tetra, but this will change (and potentially we travelled through tetras before too)
                    //transmittance stays the same, move probe until the next tetra
                    probeRay.o = probeRay.o + probeRay.d * (cached.lastHit.t + 0.01);
                    t += cached.lastHit.t;
                }
            } else {
                //otherwise the last Hit information is outdated, we are in a tetrahedron and 
                //have to accumulate the transmittance

            //TODO switch back!
               auto t2_fromt = intersectCachedTetra(probeRay,cached); 
            //auto t2_fromt = cached.lastHit.t; 

                auto transmittanceT = glm::min(tmax - t, t2_fromt);
                transmittance *= sampleTransmittance(probeRay, 0.0,transmittanceT  , cached);
                t += t2_fromt;
                probeRay.o = probeRay.o + probeRay.d * (t2_fromt + 0.01); //be sure to land in next tetra
            }
        }
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
        return scale_ * tmpVals1.vals[TF_VAL_DENS];
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


    std::unique_ptr<ArepoMesh>  arepoMesh;
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
        if (! findAndCacheTetra(cached,r.o,r.d, meshAdapter.get())) 
            toVals[TF_VAL_DENS] = 0.0; //wtf just add sth to land in tetra
        else {
            //
            lm::Float a,b,invNorm;
            sampleCachedCDFCoefficients(r,a,b,cached,invNorm);
            //only need b
            toVals[TF_VAL_DENS] = b/invNorm;
        }

    }


};


LM_COMP_REG_IMPL(Volume_Arepo_Impl, "volume::arepo");
