

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

class ArepoMeshImpl final : public lm::ArepoLMMesh {
private:
    
    ArepoMesh * am;

    std::vector<lm::Mesh::Tri> triangles;
    std::vector<int> triangleIndexToTetraIndex;
  
public:
    //TODO!
    //LM_SERIALIZE_IMPL(ar) {
    //    ar(dps, ndp, dts, ndt);
   // }

public:
    virtual void construct(const lm::Json& prop) override {

        
        //std::uintptr_t ptrToPoints = prop["ps_addr"];
        //dps = reinterpret_cast<point*>(ptrToPoints);

        am = reinterpret_cast<ArepoMesh*> ( prop["arepoMesh_addr"].get<uintptr_t>() );
        
        auto dps = am->DP;
        auto ndp = am->Ndp;
        
        //std::uintptr_t ptrToTetras = prop["ts_addr"];
        //dts = reinterpret_cast<tetra*>(ptrToPoints);;
        auto dts = am->DT;
        auto ndt = am->Ndt;

        



        triangles.reserve(ndt * 4);
        triangleIndexToTetraIndex.reserve(ndt * 4);
        for(size_t tetrai = 0; tetrai < ndt; tetrai++) {

            //"oriented" tetrahedron points...?
           // size_t vertexIndex0 = (0 + faceIndex) % 4;
            //ize_t vertexIndex1 = (1 + faceIndex) % 4;
           // size_t vertexIndex2 = (2 + faceIndex) % 4;
           // size_t vertexIndex3 = (2 + faceIndex) % 4;
            
           //auto v3 = lm::Vec3(dps[dts[tetrai].p[2]].x, dps[dts[tetrai].p[2]].y, dps[dts[tetrai].p[2]].z);

           //none of them has density?
            /*for(int i  = 0; i < 4; i++) {

            }
            int hydroIndex = arepoMeshRef->DP[cachedS.tetraInds[i]].index;
            if (hydroIndex > -1 && hydroIndex  < NumGas &&  NumGas > 0) {
*/
            //skip outside box
            if (am->DT[tetrai].t[0] < 0 || am->DT[tetrai].p[0] == DPinfinity || am->DT[tetrai].p[1] == DPinfinity
                    || am->DT[tetrai].p[2] == DPinfinity || am->DT[tetrai].p[3] == DPinfinity
                    )
            {
                // the tetra got deleted during the simulation but we have to keep the triangle count consistent, so that later face id lookups for the densities will be correct 
                /*for(int i = 0; i < 4; i++) {
                    triangles.push_back({ 
                        {lm::Vec3(0),lm::Vec3(0),lm::Vec3(-1)},
                        {lm::Vec3(0),lm::Vec3(0), lm::Vec3(-1)},
                        {lm::Vec3(0),lm::Vec3(0), lm::Vec3(-1)},
                    });
                    
                }*/
            } else {

                auto addTriangl = [&] (int tetrai, int a, int b, int c ) {
                    auto v0 = lm::Vec3(dps[dts[tetrai].p[a]].x, dps[dts[tetrai].p[a]].y, dps[dts[tetrai].p[a]].z);
                    auto v1 = lm::Vec3(dps[dts[tetrai].p[b]].x, dps[dts[tetrai].p[b]].y, dps[dts[tetrai].p[b]].z);
                    auto v2 = lm::Vec3(dps[dts[tetrai].p[c]].x, dps[dts[tetrai].p[c]].y, dps[dts[tetrai].p[c]].z);
                    auto triangleNormal = glm::cross (glm::normalize(v1 - v0),glm::normalize(v2 - v0));
                    triangles.push_back({ 
                        {v0,triangleNormal,lm::Vec2(v0.x)},
                        {v1,triangleNormal, lm::Vec2(v0.y)},
                        {v2,triangleNormal, lm::Vec2(v0.z)},
                    });
                };

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

               /*for(int i = 0; i < 4; i++) {
                    int opposingPoint = (i - 1) % 4;
                    int a = (i + 0) % 4;
                    int b = (i + 1) % 4;
                    int c = (i + 2) % 4;
                
                    

                    


                
                LM_INFO("add tetrahedron {} ? {} . triangle {},{},{} with v0x {}" ,std::to_string(tetrai), 
                std::to_string(dts[tetrai].t[0] == -1),
                std::to_string(a),std::to_string(b),std::to_string(c),
                std::to_string(dps[dts[tetrai].p[a]].x) );
                }*/
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

    int correspondingTetra(int face) const {
        return triangleIndexToTetraIndex[face];
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
        return triangles.size();
    }
};

LM_COMP_REG_IMPL(ArepoMeshImpl, "mesh::arepo");
