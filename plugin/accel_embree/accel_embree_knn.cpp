/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/accel.h>
#include <lm/scene.h>
#include <lm/mesh.h>
#include <lm/stats.h>
#include <lm/exception.h>
#include "embree_params.h"
#include <queue>
#include <lm/light.h>
LM_NAMESPACE_BEGIN(LM_NAMESPACE)


struct KNNNode {
    Vec3 global_pos;  // Global transform of the primitive
    int nodeIndex; //lm scene node index
    Float radius; //impact radius of the node
};

namespace {  

//DATA:POINT QUERY:POINT COMPARISON:LINE--------------------------------------------------------------------------------------------------

    bool pointQueryFuncLineComp(struct RTCPointQueryFunctionArguments* args)
    {
        

        KNNNode * nodes = ((KNNResult *) args->userPtr)->nodes;

        
        RTCPointQuery* query = (RTCPointQuery*)args->query;
        const unsigned int primID = args->primID;
        const lm::Vec3 q(query->x, query->y, query->z);
        KNNNode& node = nodes[primID];

        KNNResult& result = *(KNNResult*)args->userPtr;

        if(result.results.size() < result.numResults) {
            result.results.push_back(node.nodeIndex);
        } else {
            result.results[result.numResults] = node.nodeIndex;
        }
        result.numResults++;

        /*const lm::Vec3 p(node.global_pos[0],
                            node.global_pos[1],
                            node.global_pos[2]);
        //const float d = glm::distance(p, q);

        int key = 0;
        auto start_and_end = stats::get<stats::KNNLineComp,int,std::pair<Vec3,Vec3>>(key);
        auto start_to_end = start_and_end.second - start_and_end.first;
        auto minRad = glm::length(start_to_end);
        auto start_to_point = p - start_and_end.first;

        auto th = glm::dot(start_to_point, start_to_end); 
        auto shortest = start_and_end.first + start_to_end * th - p;
        auto d = glm::length(shortest);

            
        assert(args->query);
        KNNResult* result = (KNNResult*)args->userPtr;
        //result->visited.push_back(primID);

        if (d < query->radius && (result->knn.size() < result->k || d < result->knn.top().d))
        {

            Neighbour neighbour;
            neighbour.nodeIndex = node.nodeIndex;
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
        }*/
        return false;
    }
//--------------------------------------------------------------------------------------------------
    


//DATA:POINT QUERY:POINT COMPARISON:POINT--------------------------------------------------------------------------------------------------
    void pointBoundsFunc(const struct RTCBoundsFunctionArguments* args)
    {
        const std::vector<KNNNode> & nodes = * (const std::vector<KNNNode> *) args->geometryUserPtr;

        

        RTCBounds* bounds_o = args->bounds_o;
        auto primID = args->primID;
        const KNNNode& node = nodes[primID];
        

        bounds_o->lower_x = node.global_pos[0];
        bounds_o->lower_y = node.global_pos[1];
        bounds_o->lower_z = node.global_pos[2];
        bounds_o->upper_x = node.global_pos[0];
        bounds_o->upper_y = node.global_pos[1];
        bounds_o->upper_z = node.global_pos[2];

        bounds_o->lower_x -= node.radius;
        bounds_o->lower_y -=  node.radius;
        bounds_o->lower_z -=  node.radius;
        bounds_o->upper_x +=  node.radius;
        bounds_o->upper_y +=  node.radius;
        bounds_o->upper_z +=  node.radius;
    }


    bool pointQueryFunc(struct RTCPointQueryFunctionArguments* args)
    {
        KNNNode * nodes = ((KNNResult *) args->userPtr)->nodes;

        
        RTCPointQuery* query = (RTCPointQuery*)args->query;
        const unsigned int primID = args->primID;
        //const lm::Vec3 q(query->x, query->y, query->z);
        KNNNode& node = nodes[primID];
        //const lm::Vec3 p(node.global_pos[0],
        //                    node.global_pos[1],
        //                    node.global_pos[2]);
        //const float d = glm::distance(node.global_pos, q);
            
        assert(args->query);
        //Neighbour neighbour;
        //neighbour.nodeIndex = node.nodeIndex;
        //neighbour.d = d;
        KNNResult& result = *(KNNResult*)args->userPtr;
        //result->visited.push_back(primID);
        if(result.results.size() <= result.numResults) {
            result.results.push_back(node.nodeIndex);
        } else {
            result.results[result.numResults] = node.nodeIndex;
        }
        result.numResults++;
        
        /*if (d < query->radius && (result->knn.size() < result->k || d < result->knn.top().d))
        {

            
        
            if (result->knn.size() == result->k)
            result->knn.pop();
            
            result->knn.push(neighbour);
        
            if (result->knn.size() == result->k)
            {
            const float R = result->knn.top().d;
            query->radius = R;
            return true;
            }
        }*/
        return false;
    }
//--------------------------------------------------------------------------------------------------
    
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
    std::string strategy_;
    
    std::vector<KNNNode> flattened_nodes_;

public:

     virtual void construct(const Json& prop) override {
        LM_INFO("FlattenedPrimitive node size in bytes: {}", sizeof(KNNNode));
        //assert(sizeof(KNNNode) % 16 == 0);

        strategy_ = json::value<std::string>(prop, "strategy","point");


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
        //ref_flattened_nodes_ = &flattened_nodes_;

        reset();
        scene_ = rtcNewScene(device_);

        RTCGeometry geom = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_USER);
        unsigned int geomID = rtcAttachGeometry(scene_, geom);

        unsigned int lighId = 0;
        scene.traverse_primitive_nodes([&](const SceneNode& node, Mat4 global_transform) {

            if (node.type == SceneNodeType::Primitive && node.primitive.light) {


                auto lightind = scene.light_index_at(node.index);
                auto lightprim = scene.light_primitive_index_at(lightind);
                auto possmpl = node.primitive.light->sample_position({{0.0,0.0},0.0},lightprim.global_transform);
                
                assert(possmpl);
                Vec4 v = Vec4(possmpl->geom.p,1);
                Mat4 t = Mat4(1);
                t[3] = v;
                t = global_transform * t;
                
                /*LM_INFO("new node at {},{},{},{}, {},{},{},{},  {},{},{},{},  {},{},{},{}, pos {},{},{}",
                global_transform[0][0],
                global_transform[0][1],
                global_transform[0][2],
                global_transform[0][3],
                global_transform[1][0],
                global_transform[1][1],
                global_transform[1][2],
                global_transform[1][3],
                global_transform[2][0],
                global_transform[2][1],
                global_transform[2][2],
                global_transform[2][3],
                global_transform[3][0],
                global_transform[3][1],
                global_transform[3][2],
                global_transform[3][3],
                t[3][0],
                t[3][1],
                t[3][2]

                );*/
                // Record flattened primitive
                flattened_nodes_[lighId].global_pos = t[3];
                flattened_nodes_[lighId].nodeIndex = node.index;
                flattened_nodes_[lighId].radius = 0.0; //for the moment
                
                lighId++;
            }
            
        });

    
        rtcSetGeometryUserPrimitiveCount(geom, N);
        rtcSetGeometryUserData(geom, &flattened_nodes_);
        rtcSetGeometryBoundsFunction(geom, pointBoundsFunc, nullptr);

        if(strategy_ == "point")
            rtcSetGeometryPointQueryFunction(geom,pointQueryFunc);
        else if(strategy_ == "line")
            rtcSetGeometryPointQueryFunction(geom,pointQueryFuncLineComp);

        rtcCommitGeometry(geom);
        rtcReleaseGeometry(geom);

        rtcCommitScene(scene_);
    }


    virtual void build(size_t numObjects, std::function<bool(Vec3&,int&,Float&)> nextObject) override {
        exception::ScopedDisableFPEx guard_;

        /* create scene */
        size_t N = numObjects;
        flattened_nodes_.resize(N);
        LM_INFO("build accel {}",numObjects);
        //ref_flattened_nodes_ = &flattened_nodes_;

        reset();
        scene_ = rtcNewScene(device_);

        RTCGeometry geom = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_USER);
        unsigned int geomID = rtcAttachGeometry(scene_, geom);

        Vec3 global_transform;
        int someIndex;
        Float someradius;
        for ( int currentIndex = 0; currentIndex < numObjects && nextObject(global_transform,someIndex,someradius); currentIndex++){
            flattened_nodes_[currentIndex].global_pos = global_transform;
            flattened_nodes_[currentIndex].nodeIndex = someIndex;
            flattened_nodes_[currentIndex].radius = someradius;
            
        }


        

        rtcSetGeometryUserPrimitiveCount(geom, N);
        rtcSetGeometryUserData(geom, &flattened_nodes_);
        rtcSetGeometryBoundsFunction(geom, pointBoundsFunc, nullptr);
        if(strategy_ == "point")
            rtcSetGeometryPointQueryFunction(geom,pointQueryFunc);
        else if(strategy_ == "line")
            rtcSetGeometryPointQueryFunction(geom,pointQueryFuncLineComp);

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
        out_result.nodes = (KNNNode*)&flattened_nodes_[0];
        out_result.numResults = 0;
        //rtcPointQuery(scene_, &query, &context, pointQueryFunc, (void*)&out_result);
        rtcPointQuery(scene_, &query, &context,nullptr, (void*)&out_result); //out_result is the user ptr, therefore it has to contain flattened nodes omg
        
    }
    virtual std::optional<Hit> intersect(Ray ray, Float tmin, Float tmax) const {
        return {};
    }
};

LM_COMP_REG_IMPL(Accel_Embree_Knn, "accel::embree_knn");

LM_NAMESPACE_END(LM_NAMESPACE)
