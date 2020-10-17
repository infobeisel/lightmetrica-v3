

#include "arepoloader.h"
#include "snapio.h"
#include "fileio.h"
#include "geometry.h"
#include "arepo.h"

#include <lm/lm.h>
#include <lm/core.h>
#include <lm/scene.h>
#include <lm/json.h>
#include <lm/jsontype.h>
#include <lm/mesh.h>

#include <unordered_map>

class ArepoMeshImpl final : public lm::ArepoLMMesh {
private:
    
    std::vector<lm::Mesh::Tri> triangles;
    std::vector<int> triangleIndexToTetraIndex;
    std::unordered_map<int,std::vector<int>> adjacentTetras;
  
public:
    //TODO!
    //LM_SERIALIZE_IMPL(ar) {
    //    ar(dps, ndp, dts, ndt);
   // }

public:
    virtual void construct(const lm::Json& prop) override {


        
        //std::uintptr_t ptrToPoints = prop["ps_addr"];
        //dps = reinterpret_cast<point*>(ptrToPoints);

        const auto * dts = prop["tetras"].get<const ArepoLoaderInternals::Tetra *>();
        int ndt = lm::json::value<int>(prop, "nT", 0);
        
        
        //auto scale = lm::json::value<lm::Float>(prop, "scale", 1.0);


        
        triangles.clear();
        triangleIndexToTetraIndex.clear();
        triangles.reserve(ndt * 4);
        triangleIndexToTetraIndex.reserve(ndt * 4);

        adjacentTetras.clear();
        adjacentTetras.reserve(ndt);


        auto addTriangl = [&] (int tetrai, int a, int b, int c ) {
            auto v0 = lm::Vec3(dts[tetrai].positions[a].x, dts[tetrai].positions[a].y, dts[tetrai].positions[a].z);
            auto v1 = lm::Vec3(dts[tetrai].positions[b].x, dts[tetrai].positions[b].y, dts[tetrai].positions[b].z);
            auto v2 = lm::Vec3(dts[tetrai].positions[c].x, dts[tetrai].positions[c].y, dts[tetrai].positions[c].z);
            auto triangleNormal = glm::cross (glm::normalize(v1 - v0),glm::normalize(v2 - v0));
            triangles.push_back({ 
                {v0,triangleNormal,lm::Vec2(0.0)},
                {v1,triangleNormal, lm::Vec2(0.0)},
                {v2,triangleNormal, lm::Vec2(0.0)},
            });
        };
        
        for(size_t tetrai = 0; tetrai < ndt; tetrai++) {
            
            auto & t = dts[tetrai];
            //LM_INFO("MESH: Tetra {}: {}, {},{},{},{}",tetrai,t.valid, t.p[0],t.p[1],t.p[2],t.p[3]);


            //skip points 
            //wtf see arepo vtk
            DPinfinity = -1;
            //skip outside box
            if (dts[tetrai].t[0] < 0 || dts[tetrai].p[0] <= DPinfinity || dts[tetrai].p[1] <= DPinfinity
                    || dts[tetrai].p[2] <= DPinfinity || dts[tetrai].p[3] <= DPinfinity || !dts[tetrai].valid
                    )
            {
                
            } else {

                
                for(int a = 0; a < 4; a++) {
                    auto & adjacents = adjacentTetras[dts[tetrai].p[a]];
                    adjacents.push_back(tetrai);
                
                }
                

            //if(dts[tetrai].t[0] > 0.0) {
                //assume the following tetra:
                //         2-.
                //        / \ `3
                //       /   \ /
                //      0-----1
                //add the triangles in a ccw way, with normal showing inside ?(TODO)

                //MY point inside tetra works IFF triangles have the following vertex order: 
                //012, 123, 230, 301,
                //which means, they are alternating ccw and cw O.o
 
                //addTriangl(tetrai, 0,2,1);
                //addTriangl(tetrai, 0,1,2);
                addTriangl(tetrai, 0,1,2);
                triangleIndexToTetraIndex.push_back(tetrai);
                //addTriangl(tetrai, 1,2,3);
                //addTriangl(tetrai, 3,2,1);
                addTriangl(tetrai, 1,2,3);
                triangleIndexToTetraIndex.push_back(tetrai);
                //addTriangl(tetrai, 2,0,3);
                //addTriangl(tetrai, 3,0,2);
                addTriangl(tetrai, 2,3,0);
                triangleIndexToTetraIndex.push_back(tetrai);
                //addTriangl(tetrai, 3,0,1);
                //addTriangl(tetrai, 1,0,3);
                addTriangl(tetrai, 3,0,1);
                triangleIndexToTetraIndex.push_back(tetrai);

            } 
        }

        LM_INFO(" created {} triangles ",num_triangles());
    
        
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

    int correspondingTetra(int face) const {
        return triangleIndexToTetraIndex[face];
    }

    const std::vector<int> & adjacentTs(int pointIndex) override {
        
        return adjacentTetras[pointIndex];
    }

    //averages a 3
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
        return triangles.size();
    }
};

LM_COMP_REG_IMPL(ArepoMeshImpl, "mesh::arepo");

