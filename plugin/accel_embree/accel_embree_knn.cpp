/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/accel.h>
#include <lm/scene.h>
#include <lm/mesh.h>
#include <lm/exception.h>
#include "embree_params.h"
#include <queue>
LM_NAMESPACE_BEGIN(LM_NAMESPACE)


#define KNN_K 10

struct Neighbour
{
  unsigned int primID;
  float d;

  bool operator<(Neighbour const& n1) const { return d < n1.d; }
};

struct KNNResult
{
  KNNResult() 
  {
    visited.reserve(2 * KNN_K);
  }

  unsigned int k;
  std::priority_queue<Neighbour, std::vector<Neighbour>> knn;
  std::vector<unsigned int> visited; // primIDs of all visited points
};


namespace {  


struct FlattenedPrimitiveNode {
    Transform global_transform;  // Global transform of the primitive
    int primitive;              // Primitive node index
};

static FlattenedPrimitiveNode * ref_flattened_nodes_;

void pointBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const FlattenedPrimitiveNode* nodes = (const FlattenedPrimitiveNode*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    const FlattenedPrimitiveNode& node = nodes[args->primID];
    bounds_o->lower_x = node.global_transform.M[3][0];
    bounds_o->lower_y = node.global_transform.M[3][1];
    bounds_o->lower_z = node.global_transform.M[3][2];
    bounds_o->upper_x = node.global_transform.M[3][0];
    bounds_o->upper_y = node.global_transform.M[3][1];
    bounds_o->upper_z = node.global_transform.M[3][2];
}


bool pointQueryFunc(struct RTCPointQueryFunctionArguments* args)
{
  RTCPointQuery* query = (RTCPointQuery*)args->query;
  const unsigned int primID = args->primID;
  const lm::Vec3 q(query->x, query->y, query->z);
  const FlattenedPrimitiveNode& node = ref_flattened_nodes_[primID];
  const lm::Vec3 p(node.global_transform.M[3][0],
                    node.global_transform.M[3][1],
                    node.global_transform.M[3][2]);
  const float d = glm::distance(p, q);
    
  assert(args->query);
  KNNResult* result = (KNNResult*)args->userPtr;
  result->visited.push_back(primID);

  if (d < query->radius && (result->knn.size() < result->k || d < result->knn.top().d))
  {
    Neighbour neighbour;
    neighbour.primID = primID;
    neighbour.d = d;
  
    if (result->knn.size() == result->k)
      result->knn.pop();
      
    result->knn.push(neighbour);
  
    if (result->knn.size() == result->k)
    {
      const float R = result->knn.top().d;
      query->radius = R;
      return true;
    }
  }
  return false;
}

    
}

// ------------------------------------------------------------------------------------------------

/*
\rst
.. function:: accel::embree

   Acceleration structure with Embree library.
\endrst
*/
class Accel_Embree_Knn final : public AccelKnn {
private:
    RTCDevice device_ = nullptr;
    RTCScene scene_ = nullptr;
    RTCBuildArguments settings_;
    RTCSceneFlags sf_;
    
    
    std::vector<FlattenedPrimitiveNode> flattened_nodes_;

public:

     virtual void construct(const Json& prop) override {
        settings_ = prop;
        sf_ = prop;

        //check actual values used for building
        Json j;
        j = settings_;
        LM_INFO(j.dump());
        j.clear();
        j = sf_;
        LM_INFO( j.dump() );
        }

    Accel_Embree_Knn() {
        device_ = rtcNewDevice("");
        handle_embree_error(nullptr, rtcGetDeviceError(device_));
        rtcSetDeviceErrorFunction(device_, handle_embree_error, nullptr);
    }

    ~Accel_Embree_Knn() {
        if (scene_) {
            rtcReleaseScene(scene_);
        }
        if (device_) {
            rtcReleaseDevice(device_);
        }

    }

private:
    void reset() {
        if (scene_) {
            rtcReleaseScene(scene_);
            scene_ = nullptr;
        }
    }



public:
    virtual void build(const Scene& scene) override {
        exception::ScopedDisableFPEx guard_;

        /* create scene */
        size_t N = scene.num_lights();
        flattened_nodes_.resize(N);
        ref_flattened_nodes_ = &flattened_nodes_[0];

        reset();
        scene_ = rtcNewScene(device_);

        RTCGeometry geom = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_USER);
        unsigned int geomID = rtcAttachGeometry(scene_, geom);

        unsigned int lighId = 0;
        scene.traverse_primitive_nodes([&](const SceneNode& node, Mat4 global_transform) {
            if (node.type != SceneNodeType::Primitive) {
                return;
            }
            if (!node.primitive.light) {
                return;
            }
            // Record flattened primitive
            flattened_nodes_[lighId] = { Transform(global_transform), node.index };
            lighId++;
            
        });
        

        rtcSetGeometryUserPrimitiveCount(geom, N);
        rtcSetGeometryUserData(geom, &flattened_nodes_[0]);
        rtcSetGeometryBoundsFunction(geom, pointBoundsFunc, nullptr);
        rtcCommitGeometry(geom);
        rtcReleaseGeometry(geom);

        rtcCommitScene(scene_);
    }

    virtual void queryKnn(lm::Float x,lm::Float y,lm::Float z, lm::Float radius, KNNResult & out_result) const override {
        //exception::ScopedDisableFPEx guard_;
        RTCPointQuery query;
        query.x = x;
        query.y = y;
        query.z = z;
        query.radius = radius;
        query.time = 0.f;
        RTCPointQueryContext context;
        rtcInitPointQueryContext(&context);
        rtcPointQuery(scene_, &query, &context, pointQueryFunc, (void*)&out_result);
        
    }
    virtual std::optional<Hit> intersect(Ray ray, Float tmin, Float tmax) const {
        return {};
    }
};

LM_COMP_REG_IMPL(Accel_Embree_Knn, "accel::embree_knn");

LM_NAMESPACE_END(LM_NAMESPACE)
