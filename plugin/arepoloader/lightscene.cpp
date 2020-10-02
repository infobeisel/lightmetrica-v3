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

    static const Float g_band_wavelengths_nm[89] = 
    {
        3630.0,3655.0,3680.0,3705.0,3730.0,3755.0,3780.0,3805.0,3830.0,3855.0,3880.0,3905.0
        ,3930.0,3955.0,3980.0,4005.0,4030.0,4055.0,4080.0,4105.0,4130.0,4155.0,4180.0,4205.0
        ,4230.0,4255.0,4280.0,4305.0,4330.0,4355.0,4380.0,4405.0,4430.0,4455.0,4480.0,4505.0
        ,4530.0,4555.0,4580.0,4605.0,4630.0,4655.0,4680.0,4705.0,4730.0,4755.0,4780.0,4805.0
        ,4830.0,4855.0,4880.0,4905.0,4930.0,4955.0,4980.0,5005.0,5030.0,5055.0,5080.0,5105.0
        ,5130.0,5155.0,5180.0,5205.0,5230.0,5255.0,5280.0,5305.0,5330.0,5355.0,5380.0,5405.0
        ,5430.0,5455.0,5480.0,5505.0,5530.0,5555.0,5580.0,5605.0,5630.0,5655.0,5680.0,5705.0
        ,5730.0,5755.0,5780.0,5805.0,5830.0
    };

    static const Float g_band_response[89] = 
    {
        0.000e+00,3.000e-04,8.000e-04,1.300e-03,1.900e-03,2.400e-03,3.400e-03
        ,5.500e-03,1.030e-02,1.940e-02,3.260e-02,4.920e-02,6.860e-02,9.000e-02
        ,1.123e-01,1.342e-01,1.545e-01,1.722e-01,1.873e-01,2.003e-01,2.116e-01
        ,2.214e-01,2.301e-01,2.378e-01,2.448e-01,2.513e-01,2.574e-01,2.633e-01
        ,2.691e-01,2.747e-01,2.801e-01,2.852e-01,2.899e-01,2.940e-01,2.979e-01
        ,3.016e-01,3.055e-01,3.097e-01,3.141e-01,3.184e-01,3.224e-01,3.257e-01
        ,3.284e-01,3.307e-01,3.327e-01,3.346e-01,3.364e-01,3.383e-01,3.403e-01
        ,3.425e-01,3.448e-01,3.472e-01,3.495e-01,3.519e-01,3.541e-01,3.562e-01
        ,3.581e-01,3.597e-01,3.609e-01,3.613e-01,3.609e-01,3.595e-01,3.581e-01
        ,3.558e-01,3.452e-01,3.194e-01,2.807e-01,2.339e-01,1.839e-01,1.352e-01
        ,9.110e-02,5.480e-02,2.950e-02,1.660e-02,1.120e-02,7.700e-03,5.000e-03
        ,3.200e-03,2.100e-03,1.500e-03,1.200e-03,1.000e-03,9.000e-04,8.000e-04
        ,6.000e-04,5.000e-04,3.000e-04,1.000e-04,0.000e+00
    };

    static const Float g_band_weights_sum = 16.740800000000004;

    static const Vec3 cie_color_match[81] = {
        {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
        {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
        {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
        {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
        {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
        {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
        {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
        {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
        {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
        {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
        {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
        {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
        {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
        {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
        {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
        {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
        {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
        {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
        {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
        {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
        {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
        {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
        {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
        {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
        {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
        {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
        {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
    };

    static const Float sun_solid_angle = 6.793970367229902e-05;


    static const Float g_to_temp_params[10] = {
         9.17721455e+02,  2.00722937e+05, -7.69752165e+06,  1.57132181e+08,
        -1.76560740e+09,  1.17234417e+10, -4.71587498e+10,  1.12931504e+11,
        -1.48156059e+11,  8.20461651e+10
    };

    void mag_to_flux(Float & u, Float & g, Float & r, Float & i, Float & z) {

        u = glm::pow(10.0,-u/2.5) * (3631.0*1e-26);
        g = glm::pow(10.0,-g/2.5) * (3631.0*1e-26);
        r = glm::pow(10.0,-r/2.5) * (3631.0*1e-26);
        i = glm::pow(10.0,-i/2.5) * (3631.0*1e-26);
        z = glm::pow(10.0,-z/2.5) * (3631.0*1e-26);

    }


    void flux_to_mag(Float & u, Float & g, Float & r, Float & i, Float & z) {

        u = -2.5*log10(u/(3631.0*1e-26));
        g = -2.5*log10(g/(3631.0*1e-26));
        r = -2.5*log10(r/(3631.0*1e-26));
        i = -2.5*log10(i/(3631.0*1e-26));
        z = -2.5*log10(z/(3631.0*1e-26));

    }

    void flux_to_mag_g(Float & g) {

        g = -2.5*log10(g/(3631.0*1e-26));
       
    }



    Float poly10th(Float x) {
        auto ret = 0.0;
        Float x_ = 1.0;
        for(int i = 0; i < 10; i++) {
            ret += x_ * g_to_temp_params[i];
            x_ *= x;
        }
        return ret;
            
    }


    Float const h = 6.62607004e-34;
    Float const kb = 1.380649e-23;
    Float const c = 299792458.0;
    Float const hc = h * c;
    inline Float blackbody(Float wavelengthNM,Float temperatureKelvin) {
        auto fr = c / (wavelengthNM *  1e-9);
        auto t1 = 2.0 * h *fr *fr*fr;
        auto t2 = 1.0 / (c*c* (glm::exp(h*fr/(kb * temperatureKelvin))-1.0));
        return  t1 * t2;
    }



    Float applySdssGBand(Float blackBodyTemperature) {
        Float ret = 0.0;
        for(int i = 0; i < 89; i++) {
            ret += blackbody(g_band_wavelengths_nm[i], blackBodyTemperature) * g_band_response[i];
        }
        ret /= g_band_weights_sum;
        return ret;
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
        //f_coefs_g_to_temp_ = lm::json::value<std::vector<Float>>(prop,"f_coefs_g_to_temp");

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
            auto u_mag = star_ugriz_[ugriz_i + 0];
            auto g_mag = star_ugriz_[ugriz_i + 4];
            auto r_mag = star_ugriz_[ugriz_i + 5];
            auto i_mag = star_ugriz_[ugriz_i + 6];
            auto z_mag = star_ugriz_[ugriz_i + 7];

            auto u_flux = u_mag;
            auto g_flux = g_mag;
            auto r_flux = r_mag;
            auto i_flux = i_mag;
            auto z_flux = z_mag;

            mag_to_flux(u_flux,g_flux,r_flux,i_flux,z_flux);

            auto sum = g_flux +  r_flux +  i_flux +  z_flux;
            auto g_ratio = g_flux / sum;
            auto temperature = poly10th(g_ratio);

            auto g_response_flux = sun_solid_angle*applySdssGBand(temperature);
            auto scaling = g_flux / g_response_flux;

            auto reproduced_mag = scaling * g_response_flux;
            flux_to_mag_g(reproduced_mag);

            //LM_INFO("incoming mag : {},{},{},{}, reproduced g mag: {}",
            //g_mag,r_mag,i_mag,z_mag,
           // reproduced_mag);

            Float x_,y_,z_;
            tempToXYZ(temperature,x_,y_,z_);
            auto rgb = M * Vec3(x_,y_,z_);

            lightprop["Le"] = scaling * Vec3(rgb.r,rgb.g,rgb.b);
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
