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




namespace stats {
    struct CachedSampleId {};
    struct SampleIdCacheHits {};
    struct SampleIdCacheMisses {};
    struct UsedCachedTetra {};
    struct UsedNeighborTetra {};
    struct InvalidNeighbor {};
    struct ResampleAccel {};
    struct TotalTetraTests {};

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
    virtual Float max_scalar(lm::Ray ray, lm::Float & out_t_forhowlong) const = 0;
    virtual Float eval_scalar(Vec3 p) const = 0;
    virtual bool has_color() const = 0;
    virtual Vec3 eval_color(Vec3 p) const = 0;
    virtual Float sample_distance(Ray ray,lm::Float tmin, lm::Float tmax, lm::Rng& rng) const = 0;
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