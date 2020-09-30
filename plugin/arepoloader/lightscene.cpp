/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/core.h>
#include <lm/scene.h>
#include <lm/accel.h>
#include <lm/mesh.h>
#include <lm/camera.h>
#include <lm/material.h>
#include <lm/light.h>
#include <lm/model.h>
#include <lm/medium.h>
#include <lm/phase.h>
#include <lm/stats.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)
namespace stats {
    struct SelectedLightPDF{};
}
namespace {
    enum LightSamplingMode {
        UNIFORM,
        KNN
    };


    const Float bTable_0 = 1.4*glm::pow(10.0,-10.0);
    const Float bTable_1 = 0.9*glm::pow(10.0,-10.0);
    const Float bTable_2 = 1.2*glm::pow(10.0,-10.0);
    const Float bTable_3 = 1.8*glm::pow(10.0,-10.0);
    const Float bTable_4 = 7.4*glm::pow(10.0,-10.0);
    
    void ugriz_to_jy(Float & u, Float & g, Float & r, Float & i, Float & z) {
        u = 3631.0 * 2.0 * bTable_0 * glm::sinh(-u*glm::log(10)/2.5 - glm::log(bTable_0));
        g = 3631.0 * 2.0 * bTable_1 * glm::sinh(-g*glm::log(10)/2.5 - glm::log(bTable_1));
        r = 3631.0 * 2.0 * bTable_2 * glm::sinh(-r*glm::log(10)/2.5 - glm::log(bTable_2));
        i = 3631.0 * 2.0 * bTable_3 * glm::sinh(-i*glm::log(10)/2.5 - glm::log(bTable_3));
        z = 3631.0 * 2.0 * bTable_4 * glm::sinh(-z*glm::log(10)/2.5 - glm::log(bTable_4));
    }

    Float poly10th(std::vector<Float> & coefs,Float x) {
        auto ret = 0.0;
        Float x_ = 1.0;
        for(int i = 0; i < 10; i++) {
            ret += x_ * coefs[i];
            x_ *= x;
        }
        return ret;
            
    }


    Float const h = 6.62607004e-34;
    Float const stb = 1.380649e-23;
    Float const c = 299792458.0;
    Float const hc = h * c;
    inline Float blackbody(Float wavelengthNM,Float temperatureKelvin) {
        Float lambd = wavelengthNM *  1e-9;//m
        Float t = 2.0 * hc*c/(lambd*lambd*lambd*lambd*lambd);
        Float t2 = 1.0 / (glm::exp(hc/(lambd * stb * temperatureKelvin))-1.0);
        return t * t2;
    }

    inline Float xFit_1931( Float wave )
    {
        Float t1 = (wave-442.0f)*((wave<442.0f)?0.0624f:0.0374f);
        Float t2 = (wave-599.8f)*((wave<599.8f)?0.0264f:0.0323f);
        Float t3 = (wave-501.1f)*((wave<501.1f)?0.0490f:0.0382f);
        return 0.362f*glm::exp(-0.5f*t1*t1) + 1.056f*glm::exp(-0.5f*t2*t2)
        - 0.065f*glm::exp(-0.5f*t3*t3);
    }
    inline Float yFit_1931( Float wave )
    {
        Float t1 = (wave-568.8f)*((wave<568.8f)?0.0213f:0.0247f);
        Float t2 = (wave-530.9f)*((wave<530.9f)?0.0613f:0.0322f);
        return 0.821f*glm::exp(-0.5f*t1*t1) + 0.286f*glm::exp(-0.5f*t2*t2);
    }
    inline Float zFit_1931( Float wave )
    {
        Float t1 = (wave-437.0f)*((wave<437.0f)?0.0845f:0.0278f);
        Float t2 = (wave-459.0f)*((wave<459.0f)?0.0385f:0.0725f);
        return 1.217f*glm::exp(-0.5f*t1*t1) + 0.681f*glm::exp(-0.5f*t2*t2);
    }

    inline void tempToXYZ(Float temp, Float & X,Float & Y,Float & Z) {
        X = Y = Z = 0.0;
        for(int nmi = 340; nmi < 900; nmi++) {
            Float nm = (Float)nmi;
            X += xFit_1931(nm) * blackbody(nm,temp);
            Y += yFit_1931(nm) * blackbody(nm,temp);
            Z += zFit_1931(nm) * blackbody(nm,temp);
        }
    }

        
            
    
}

class LightScene_ final : public Scene {
private:
    AccelKnn* accel_;                                   // Acceleration structure for closest light sampling
    std::vector<SceneNode> nodes_;                   // Scene nodes (index 0: root node)
    std::optional<int> camera_;                      // Camera index
    std::vector<LightPrimitiveIndex> lights_;        // Primitive node indices of lights and global transforms
    std::unordered_map<int, int> light_indices_map_; // Map from node indices to light indices.
    std::optional<int> env_light_;                   // Environment light index
    std::optional<int> medium_;                      // Medium index
    LightSamplingMode lightsamplingmode;
    std::vector<Float> star_coords_;
    std::vector<Float> star_ugriz_;
    std::vector<Float> f_coefs_g_to_temp_;
    std::vector<Float> f_coefs_temp_to_x_;
    std::vector<Float> f_coefs_temp_to_y_;
    std::vector<Float> f_coefs_temp_to_z_;
    std::vector<Ptr<Light>> createdLights_;
public:
    LM_SERIALIZE_IMPL(ar) {
        ar(accel_, nodes_, camera_, lights_, light_indices_map_, env_light_,createdLights_);
    }


    virtual Component* underlying(const std::string& name) const override {
        for(auto & l : createdLights_) {
            if(l->name() == name)
                return l.get();
        }
        return nullptr;
                // Underlying component must be accessible with the same name specified in create function
       // return name == comp->name() ? comp.get() : nullptr;
    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        comp::visit(visit, accel_);
        for (auto& node : nodes_) {
            if (node.type == SceneNodeType::Primitive) {
                comp::visit(visit, node.primitive.mesh);
                comp::visit(visit, node.primitive.material);
                comp::visit(visit, node.primitive.light);
                comp::visit(visit, node.primitive.camera);
            }
        }
        for(auto & light : createdLights_) {
            auto * p =  light.get();
            comp::visit(visit,p);
        }
    }

public:
    virtual void construct(const Json& prop) override {
        accel_ = json::comp_ref_or_nullptr<AccelKnn>(prop, "accel");
        lightsamplingmode = lm::json::value<std::string>(prop,"lightsampling","knn") == "knn" ? 
        LightSamplingMode::KNN : LightSamplingMode::UNIFORM;
        star_coords_ = lm::json::value<std::vector<Float>>(prop,"star_coords");
        star_ugriz_ = lm::json::value<std::vector<Float>>(prop,"star_ugriz");
        f_coefs_g_to_temp_ = lm::json::value<std::vector<Float>>(prop,"f_coefs_g_to_temp");
        //f_coefs_temp_to_x_ = lm::json::value<std::vector<Float>>(prop,"f_coefs_temp_to_x");
        //f_coefs_temp_to_y_ = lm::json::value<std::vector<Float>>(prop,"f_coefs_temp_to_y");
        //f_coefs_temp_to_z_ = lm::json::value<std::vector<Float>>(prop,"f_coefs_temp_to_z");

        

        //to sRGB
        auto m0 = Vec3(3.2404542,-0.9692660,0.0556434);     
        auto m1 = Vec3( -1.5371385, 1.8760108,-0.2040259);  
        auto m2 = Vec3(  -0.4985314, 0.0415560, 1.0572252);
        auto M = Mat3(m0,m1,m2);

        LM_INFO("star coords: {}, star photos {}",star_coords_.size()/3,star_ugriz_.size()/8);
        createdLights_.reserve(star_coords_.size() / 3);

        Json lightprop;
        for(int ind = 0; ind < star_coords_.size() / 3; ind++) {
            auto ugriz_i = ind * 8;
            auto coord_i = ind * 3;
            auto u = star_ugriz_[ugriz_i + 0];
            auto g = star_ugriz_[ugriz_i + 4];
            auto r = star_ugriz_[ugriz_i + 5];
            auto i = star_ugriz_[ugriz_i + 6];
            auto z = star_ugriz_[ugriz_i + 7];
            ugriz_to_jy(u,g,r,i,z);
            auto sum = g + r + i + z; //do not use u as it disturbs the fit
            auto temperature = poly10th(f_coefs_g_to_temp_,g / sum);



            //auto x_ = poly10th(f_coefs_temp_to_x_,temperature);
            //auto y_ = poly10th(f_coefs_temp_to_y_,temperature);
            //auto z_ = poly10th(f_coefs_temp_to_z_,temperature);
            Float x_,y_,z_;
            tempToXYZ(temperature,x_,y_,z_);
            auto rgb = M * Vec3(x_,y_,z_);

            lightprop["Le"] = Vec3(rgb.r,rgb.g,rgb.b);
            lightprop["position"] = Vec3(star_coords_[coord_i],star_coords_[coord_i+1],star_coords_[coord_i+2]);
            std::string name = "light" + std::to_string(ind);
            createdLights_.push_back(
                std::move(
                    comp::create<Light>(
                        "light::point", make_loc(name), lightprop)
                )
            );
            /*LM_INFO("star with jansky {},{},{},{},  with temp {}, rgb {},{},{}, at pos {},{},{}",
            g,r,i,z,
            temperature,
            rgb.r,
            rgb.g,
            rgb.b,
            star_coords_[coord_i],
            star_coords_[coord_i+1],star_coords_[coord_i+2]
            );*/

        }
        
        /*
        for (int i = 0; i < star_coords_.size(); i+=3) {
            LM_INFO("coord {},{},{}",star_coords_[i],star_coords_[i+1],star_coords_[i+2]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }
        for (int i = 0; i < star_ugriz_.size(); i+=5) {
            LM_INFO("ugriz {},{},{},{},{}",
            star_ugriz_[i],star_ugriz_[i+1],star_ugriz_[i+2],star_ugriz_[i+3],star_ugriz_[i+4]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }
        for (int i = 0; i < f_coefs_g_to_temp_.size(); i+=1) {
            LM_INFO("f_coefs_g_to_temp_ {}",f_coefs_g_to_temp_[i]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }


        for (int i = 0; i < f_coefs_temp_to_x_.size(); i+=1) {
            LM_INFO("f_coefs_temp_to_x_ {}",f_coefs_temp_to_x_[i]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }

        for (int i = 0; i < f_coefs_temp_to_y_.size(); i+=1) {
            LM_INFO("f_coefs_temp_to_y_ {}",f_coefs_temp_to_y_[i]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }
        for (int i = 0; i < f_coefs_temp_to_z_.size(); i+=1) {
            LM_INFO("f_coefs_temp_to_z_ {}",f_coefs_temp_to_z_[i]);
            //ps_.push_back(Vec3(ps[i],ps[i+1],ps[i+2]));
        }*/
        



        reset();
    }

public:
    virtual void reset() override {
        nodes_.clear();
        camera_ = {};
        lights_.clear();
        light_indices_map_.clear();
        env_light_ = {};
        medium_ = {};
        nodes_.push_back(SceneNode::make_group(0, false, {}));
    }

    // --------------------------------------------------------------------------------------------

    #pragma region Scene graph manipulation and access

    virtual int root_node() override {
        return 0;
    }

    virtual int create_primitive_node(const Json& prop) override {
        // Find an asset by property name
        const auto get_asset_ref_by = [&](const std::string& propName) -> Component * {
            const auto it = prop.find(propName);
            if (it == prop.end()) {
                return nullptr;
            }
            return comp::get<Component>(it.value().get<std::string>());
        };

        // Node index
        const int index = int(nodes_.size());

        // Get asset references
        auto* mesh = dynamic_cast<Mesh*>(get_asset_ref_by("mesh"));
        auto* material = dynamic_cast<Material*>(get_asset_ref_by("material"));
        auto* light = dynamic_cast<Light*>(get_asset_ref_by("light"));
        auto* camera = dynamic_cast<Camera*>(get_asset_ref_by("camera"));
        auto* medium = dynamic_cast<Medium*>(get_asset_ref_by("medium"));

        // Check validity
        if (!mesh && !material && !light && !camera && !medium) {
            LM_THROW_EXCEPTION(Error::InvalidArgument, "Invalid primitive node. Given assets are invalid.");
            return -1;
        }
        if (camera && light) {
            LM_THROW_EXCEPTION(Error::InvalidArgument, "Primitive cannot be both camera and light.");
            return -1;
        }

        // If you specify the mesh, you also need to specify the material
        if ((!mesh && material) || (mesh && !material)) {
            LM_THROW_EXCEPTION(Error::InvalidArgument, "You must specify both mesh and material.");
            return -1;
        }

        // Camera
        if (camera) {
            camera_ = index;
        }

        // Envlight
        if (light && light->is_env()) {
            if (env_light_) {
                LM_THROW_EXCEPTION(Error::InvalidArgument, "Environment light is already registered. "
                    "You can register only one environment light in the scene.");
                return -1;
            }
            env_light_ = index;
        }

        // Medium
        if (medium) {
            // For now, consider the medium as global asset.
            medium_ = index;
        }

        // Create primitive node
        nodes_.push_back(SceneNode::make_primitive(index, mesh, material, light, camera, medium));

        return index;
    }

    virtual int create_group_node(Mat4 transform) override {
        const int index = int(nodes_.size());
        nodes_.push_back(SceneNode::make_group(index, false, transform));
        return index;
    }

    virtual int create_instance_group_node() override {
        const int index = int(nodes_.size());
        nodes_.push_back(SceneNode::make_group(index, true, {}));
        return index;
    }

    virtual void add_child(int parent, int child) override {
        if (parent < 0 || parent >= int(nodes_.size())) {
            LM_ERROR("Missing parent index [index='{}'", parent);
            return;
        }

        auto& node = nodes_.at(parent);
        if (node.type != SceneNodeType::Group) {
            LM_ERROR("Adding child to non-group node [parent='{}', child='{}']", parent, child);
            return;
        }

        node.group.children.push_back(child);
    }

    virtual void add_child_from_model(int parent, const std::string& modelLoc) override {
        if (parent < 0 || parent >= int(nodes_.size())) {
            LM_ERROR("Missing parent index [index='{}'", parent);
            return;
        }

        auto* model = comp::get<Model>(modelLoc);
        if (!model) {
            return;
        }

        model->create_primitives([&](Component* mesh, Component* material, Component* light) {
            const int index = int(nodes_.size());
            nodes_.push_back(SceneNode::make_primitive(
                index,
                dynamic_cast<Mesh*>(mesh),
                dynamic_cast<Material*>(material),
                dynamic_cast<Light*>(light),
                nullptr,
                nullptr));
            add_child(parent, index);
        });
    }

    virtual int create_group_from_model(const std::string& modelLoc) override {
        auto* model = comp::get<Model>(modelLoc);
        if (!model) {
            LM_ERROR("Invalid model [loc={}]", modelLoc);
            return -1;
        }

        // Iterate underlyng scene nodes in the model
        // To copy the underlying scene graph as a child of scene graph of the scene,
        // we first want to copy the nodes and modify the references to the other nodes.
        // Actually, we merely want to offset the indices by nodes_.size().
        const int offset = int(nodes_.size());
        model->foreach_node([&](const SceneNode& modelNode) {
            // Copy the node
            nodes_.push_back(modelNode);
            auto& node = nodes_.back();

            // Update the references to other nodes
            node.index += offset;
            if (node.type == SceneNodeType::Group) {
                for (int& child : node.group.children) {
                    child += offset;
                }
            }
            else if (node.type == SceneNodeType::Primitive) {
                if (node.primitive.camera) {
                    camera_ = node.index;
                }
            }
        });

        return offset;
    }

    virtual void traverse_primitive_nodes(const NodeTraverseFunc& traverseFunc) const override {
        std::function<void(int, Mat4)> visit = [&](int index, Mat4 global_transform) {
            const auto& node = nodes_.at(index);
            traverseFunc(node, global_transform);
            if (node.type == SceneNodeType::Group) {
                const auto M = node.group.local_transform
                    ? global_transform * *node.group.local_transform
                    : global_transform;
                for (int child : node.group.children) {
                    visit(child, M);
                }
            }
        };
        visit(0, Mat4(1_f));
    }

    virtual void visit_node(int node_index, const VisitNodeFunc& visit) const override {
        visit(nodes_.at(node_index));
    }

    virtual const SceneNode& node_at(int node_index) const override {
        return nodes_.at(node_index);
    }

    virtual int num_nodes() const override {
        return (int)(nodes_.size());
    }

    virtual int env_light_node() const override {
        return env_light_ ? *env_light_ : -1;
    }

    virtual int camera_node() const override {
        return camera_ ? *camera_ : -1;
    }

    virtual int medium_node() const override {
        return medium_ ? *medium_ : -1;
    }

    virtual int num_lights() const override {
        return (int)(lights_.size());
    }

    #pragma endregion

    // --------------------------------------------------------------------------------------------

    #pragma region Ray-scene intersection

    virtual Accel* accel() const override {
        return accel_;
    }

    virtual void set_accel(const std::string& accel_loc) override {
        accel_ = comp::get<AccelKnn>(accel_loc);
    }

    virtual void build() override {
        // Update light indices
        // We keep the global transformation of the light primitive as well as the references.
        // We need to recompute the indices when an update of the scene happens,
        // because the global tranformation can only be obtained by traversing the nodes.
        light_indices_map_.clear();
        lights_.clear();

        for(auto & l : createdLights_)
            add_primitive({
                {"light",l->loc()}
            });

        traverse_primitive_nodes([&](const SceneNode& node, Mat4 global_transform) {
            if (node.type == SceneNodeType::Primitive && node.primitive.light) {
                light_indices_map_[node.index] = int(lights_.size());
                lights_.push_back({ Transform(global_transform), node.index });
            }
        });

        

        // Compute scene bound
        Bound bound;
        traverse_primitive_nodes([&](const SceneNode& node, Mat4 global_transform) {
            if (node.type != SceneNodeType::Primitive) {
                return;
            }
            if (!node.primitive.mesh) {
                return;
            }
            node.primitive.mesh->foreach_triangle([&](int, const Mesh::Tri& tri) {
                const auto p1 = global_transform * Vec4(tri.p1.p, 1_f);
                const auto p2 = global_transform * Vec4(tri.p2.p, 1_f);
                const auto p3 = global_transform * Vec4(tri.p3.p, 1_f);
                bound = merge(bound, p1);
                bound = merge(bound, p2);
                bound = merge(bound, p3);
            });
        });
        
        // Set scene bound to the lights
        for (auto& l : lights_) {
            auto* light = nodes_.at(l.index).primitive.light;
            light->set_scene_bound(bound);
        }

        // Build acceleration structure
        LM_INFO("Building acceleration structure [name='{}']", accel_->name());
        LM_INDENT();
        accel_->build(*this);
    }

    virtual std::optional<SceneInteraction> intersect(Ray ray, Float tmin, Float tmax) const override {
        return {};
    }

    #pragma endregion

    // --------------------------------------------------------------------------------------------

    #pragma region Primitive type checking

    virtual bool is_light(const SceneInteraction& sp) const override {
        const auto& primitive = nodes_.at(sp.primitive).primitive;
        if (sp.is_type(SceneInteraction::MediumInteraction)) {
            return primitive.medium->is_emitter();
        }
        else {
            // Note: SurfaceInterfaction can contain light component.
            return primitive.light != nullptr;
        }
    }

    virtual bool is_camera(const SceneInteraction& sp) const override {
        return sp.primitive == *camera_;
    }

    #pragma endregion

    // --------------------------------------------------------------------------------------------

    #pragma region Light sampling


    virtual LightSelectionSample sample_light_selection_from_pos(Float u, Vec3 const pos) const override {
        if(lightsamplingmode == LightSamplingMode::KNN) {
            const int n = int(lights_.size());
            int onepercent = int(0.01 * Float(n));  
            onepercent = onepercent > 0 ? onepercent : 1;
            KNNResult knnres;

/*
            unsigned int k;
  std::priority_queue<Neighbour, std::vector<Neighbour>> knn;
  std::vector<unsigned int> visited; // primIDs of all visited points


unsigned int primID;
  float d;
*/

            int key = 0;
            auto & visitor = *stats::get<stats::LightKnnVisitor,int,std::function<void(int,Float)>*>(key);

            knnres.k = onepercent;
            accel_->queryKnn(pos.x,pos.y,pos.z, std::numeric_limits<Float>::max(), knnres);
            
            std::vector<Neighbour> sorted;
            sorted.reserve(knnres.knn.size());
            //LM_INFO("found {} nearest lights",knnres.knn.size());
            Float d_total = 0.0;
            while(!knnres.knn.empty()) {
                auto & neighbor = knnres.knn.top();
                d_total += neighbor.d;
                sorted.push_back(neighbor);
                knnres.knn.pop();
            }

            


            Float cdf = 0.0;
            int finalSelectedNodeIndex = 0;
            Float pL = onepercent   ;
            bool stopUpdating = false;
            for(auto & neighbor : sorted) {
                Float d = neighbor.d / d_total;
                visitor(neighbor.nodeIndex,d);
                cdf += d;
                if(u < cdf && !stopUpdating) {
                    finalSelectedNodeIndex = neighbor.nodeIndex;
                    pL = d / d_total;
                    stopUpdating = true;
                }
            }
            //LM_INFO("chose light index {} with pdf {}",finalSelectedNodeIndex,pL);
            stats::set<stats::SelectedLightPDF,int,Vec3>(finalSelectedNodeIndex,pos);
            stats::set<stats::SelectedLightPDF,int,Float>(finalSelectedNodeIndex,pL);

            const int i = glm::clamp(int(u * static_cast<Float>(sorted.size())), 0, static_cast<int>(sorted.size()) - 1);
            pL = 1_f / n;
            return LightSelectionSample{
                i,
                pL
            };


            //return LightSelectionSample{
            //    finalSelectedNodeIndex,
            //    pL
            //};
        }

        if(lightsamplingmode == LightSamplingMode::UNIFORM) {
            const int n = int(lights_.size());
            const int i = glm::clamp(int(u * n), 0, n - 1);
            const auto pL = 1_f / n;
            return LightSelectionSample{
                i,
                pL
            };
        }
    }


    virtual LightSelectionSample sample_light_selection(Float u) const override {
        LM_ERROR("dont use this");
        const int n = int(lights_.size());
        const int i = glm::clamp(int(u * n), 0, n - 1);
        const auto pL = 1_f / n;
        return LightSelectionSample{
            i,
            pL
        };
    
    }


    virtual Float pdf_light_selection(int finalSelectedNodeIndex) const override {
        if(lightsamplingmode == LightSamplingMode::UNIFORM) {
            const int n = int(lights_.size());
            return 1_f / n;
        } 
        if(lightsamplingmode == LightSamplingMode::KNN) {
            LM_ERROR("could not implement");
           
            
            return 0.0;
        } 
        
        
    };



    virtual Float pdf_light_selection_from_pos(int light_index, Vec3 const pos) const override {
        if(lightsamplingmode == LightSamplingMode::UNIFORM) {
            return pdf_light_selection(light_index);
        } 
        if(lightsamplingmode == LightSamplingMode::KNN) {
            
            LM_ERROR("could not implement, but need to");

            return 0.0;
        } 
        
        
    };


    virtual LightPrimitiveIndex light_primitive_index_at(int light_index) const override {
        return lights_.at(light_index);
    }

    virtual int light_index_at(int node_index) const override {
        return light_indices_map_.at(node_index);
    }

    #pragma endregion
};

LM_COMP_REG_IMPL(LightScene_, "scene::arepo");

LM_NAMESPACE_END(LM_NAMESPACE)
