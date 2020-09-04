#include <lm/core.h>
#include <lm/lm.h>
#include <lm/volume.h>
#include <lm/stats.h>


//#include "../../../ArepoVTK/arepo/include/mesh/voronoi/voronoi.h"
#include "voronoi_3db.h" //shit didnt want to include any arepo code here


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

    struct RegularTrackingStrategyDistanceSample {};
    struct EquiangularStrategyDistanceSample {};

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

}
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