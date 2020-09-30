#include <lm/core.h>
#include <lm/lm.h>
#include <lm/volume.h>
#include <lm/stats.h>


//#include "../../../ArepoVTK/arepo/include/mesh/voronoi/voronoi.h"
#include "voronoi_3db.h" //shit didnt want to include any arepo code here


//scaling factors for optical thickness
#define A_B_A_V_T 1.324
#define A_R_A_V_T 0.748


//scalling factors for scattering coefficient
#define A_B_A_V_S 3.0343721 
#define A_R_A_V_S 0.49774692


#define INSIDE_TOLERANCE 10000.0 * std::numeric_limits<lm::Float>::epsilon()

inline lm::Float sampleCachedICDF_andCDF(lm::Float logxi,lm::Float xi, lm::Float tmax, lm::Float & out_cdf,lm::Float a, lm::Float b){
    

    //use tau*_t1 (t) which is the integral from t1 to t minus the integral from 0 to t1
    //auto y = logxi + a *  0.5 * tmin*tmin   + b * tmin  - out_cdf;
    auto y = logxi - out_cdf;
    lm::Float freeT;

    if (abs(a) < INSIDE_TOLERANCE
    && abs(b) < INSIDE_TOLERANCE) {
        freeT = std::numeric_limits<lm::Float>::max();
    } else if(abs(a) < INSIDE_TOLERANCE) {
        //freeT = 2.0 * y / (glm::sqrt( b*b +  2.0 * y * a  ) + b);
        freeT = y / b;
    }
    else if(abs(b) < INSIDE_TOLERANCE) {
        //freeT = (glm::sqrt( b*b +  2.0 * y * a  ) - b) / a;
        freeT = glm::sqrt(2.0 * y * a) / a;
    } else if (abs(a) < abs(b)) {
        freeT = 2.0 * y / (glm::sqrt( b*b +  2.0 * y * a  ) + b);
    } else if (abs(b) < abs(a)) {
        freeT = (glm::sqrt( b*b +  2.0 * y * a  ) - b) / a;
    }
  
    freeT = isnan(freeT) ? std::numeric_limits<lm::Float>::max() : freeT;
    //for evaluating tau, limit free path to tmax
    lm::Float t = glm::min(freeT , tmax);
    lm::Float acc_cdf = t  * b  +  a  * 0.5 * t * t;
    
    out_cdf += glm::max(0.0,
        acc_cdf) ; //cdf within tmin and  min of (t , tmax) 
    return freeT;//returns sth between tmin - tmin (so 0) and tmax - tmin 
}
inline lm::Float sampleCDF(  lm::Float toT,lm::Float a, lm::Float b) {
    return ( b * toT  +  a * 0.5 *toT * toT   );
}

namespace ArepoLoaderInternals {
    struct IArepoMeshMock {
        virtual tetra * getDT() = 0;
        virtual point * getDP() = 0;
        virtual std::vector<lm::Float> & getdensities() = 0;
        virtual int getNdt() = 0;
        virtual int getNdp() = 0;

        virtual BBox WorldBound() = 0;
        virtual lm::Float getDensity(int index) = 0;
        virtual lm::Float max_density() = 0;
     };
}



LM_NAMESPACE_BEGIN(LM_NAMESPACE)



typedef std::vector<std::vector<lm::Float>&> stdvec2d;
    

struct LightToCameraRaySegmentCDF {
    Vec3 weight;
    Vec3 p,d;
    lm::Float cdfSoFar;
    lm::Float localcdf;
    lm::Float t;
    lm::Float tSoFar;
    lm::Float a;
    lm::Float b;
    

};
void to_json(lm::Json& j, const LightToCameraRaySegmentCDF& p); 


void from_json(const lm::Json& j, LightToCameraRaySegmentCDF& p);

namespace stats {
    struct CachedSampleId {};
    struct SampleIdCacheHits {};
    struct SampleIdCacheMisses {};
    struct UsedCachedTetra {};
    struct UsedNeighborTetra {};
    struct InvalidNeighbor {};
    struct ResampleAccel {};
    struct TotalTetraTests {};

    struct MaxTransmittance {};
    struct FreePathTransmittance {};
    struct OpticalThickness {};

    struct RegularTrackingStrategyDistanceSample {};

    struct RegularTrackingStrategyTotalT{};//the distance until the last volume boundary
    struct RegularTrackingStrategyMinT{}; //the distance until the first volume boundary

    struct RegularTrackingStrategyXi {};
    struct RegularTrackingStrategyTotalTau{}; //the optical thickness for the whole current ray
    struct RegularTrackingStrategyTauUntilScatter{}; //the optical thickness accumulated until the regular distance sample scattered
    struct RegularTrackingStrategyNormFac {};//the normalization factor to make all uniform samples scatter 

    struct EquiangularStrategyDistanceSample {};

    struct ScatteringAlbedo {};


    struct DistanceSampleRandomValues {};
    struct EquiDistanceSampleRandomValueVertexIndex {};
    struct RegularDistanceSampleRandomValueVertexIndex {};

    struct EquiContribution;
    struct EquiEquiPDF;

    struct DistanceSamplesPDFs{};
    //2 distance samples, first one from equiangular, second one from regular 
    enum IJ {
        //equiangular sampling, distance sample 0
        _0_0,
        //equiangular sampling, distance sample 1
        _0_1,
        //regular sampling, distance sample 0
        _1_0,
        //regular sampling, distance sample 1
        _1_1
    };

    struct BoundaryVisitor{};
    struct LastBoundarySequence {};

    struct VRL{};
    typedef int TetraIndex;



    




}

struct RaySegmentCDF {
    lm::Float localcdf;
    lm::Float t;
    lm::Float a;
    lm::Float b;
    int tetraI;
};




    
class ArepoLMMesh : public Mesh {

public:
    virtual void foreach_triangle(const ProcessTriangleFunc& process_triangle) const = 0;

    virtual lm::Mesh::Tri triangle_at(int face) const = 0;

    virtual int correspondingTetra(int face) const = 0;

    virtual lm::Mesh::InterpolatedPoint surface_point(int face, lm::Vec2 uv) const = 0;

    virtual int num_triangles() const = 0;
};



class Volume_Arepo : public Volume {
public:
    virtual Bound bound() const = 0;
    virtual bool has_scalar() const = 0;
    virtual Float max_scalar() const = 0;
    virtual Float max_scalar(lm::Ray ray, lm::Float & out_t_forhowlong, lm::Float & out_aCoef, lm::Float & out_bCoef) const = 0;
    virtual Float eval_scalar(Vec3 p) const = 0;
    virtual Float eval_scalar(Vec3 p , Vec3 dir) const = 0;
    virtual bool has_color() const = 0;
    virtual Vec3 eval_color(Vec3 p) const = 0;
    virtual Float sample_distance(Ray ray,lm::Float tmin, lm::Float tmax, lm::Rng& rng,lm::Float & weight) const = 0;
    virtual Float eval_transmittance(lm::Ray ray, Float tmin, Float tmax) const = 0;

};


LM_NAMESPACE_END(LM_NAMESPACE)

/*
LM_NAMESPACE_BEGIN(LM_NAMESPACE::objloader)



class Mesh_Arepo : public OBJLoaderContext {
public:
    public:
    virtual bool load(
        const std::string& path,
        OBJSurfaceGeometry& geo,
        const ProcessMeshFunc& process_mesh,
        const ProcessMaterialFunc& process_material) override = 0;
      
};
LM_NAMESPACE_END(LM_NAMESPACE::objloader)

*/