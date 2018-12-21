/*
    Lightmetrica - Copyright (c) 2018 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/scene.h>
#include <lm/assets.h>
#include <lm/accel.h>
#include <lm/mesh.h>
#include <lm/camera.h>
#include <lm/material.h>
#include <lm/light.h>
#include <lm/model.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

struct Primitive {
    int index;
    Mat4 transform; 
    const Mesh* mesh;
    const Material* material;
    const Light* light;
    const Camera* camera;
};

class Scene_ final : public Scene {
public:
    virtual bool loadPrimitive(const Component& assetGroup, Mat4 transform, const Json& prop) override {
        const auto getAssetRefBy = [&](const std::string& propName) -> Component* {
            if (propName.empty()) {
                return nullptr;
            }
            // Find reference to an asset
            const auto it = prop.find(propName);
            if (it == prop.end()) {
                return nullptr;
            }
            // Obtain the referenced asset
            return assetGroup.underlying(it.value().get<std::string>().c_str());
        };
        primitives_.push_back(Primitive{
            int(primitives_.size()),
            transform,
            getAssetRefBy("mesh")->cast<Mesh>(),
            getAssetRefBy("material")->cast<Material>(),
            getAssetRefBy("light")->cast<Light>(),
            getAssetRefBy("camera")->cast<Camera>()
        });
        if (primitives_.back().camera) {
            camera_ = primitives_.back().index;
        }
        return true;
    }

    virtual bool loadPrimitives(const Component& assetGroup, Mat4 transform, const std::string& modelName) override {
        auto* model = assetGroup.underlying<Model>(modelName);
        if (!model) {
            return false;
        }
        model->createPrimitives([&](Component* mesh, Component* material, Component* light) {
            primitives_.push_back(Primitive{
                int(primitives_.size()),
                transform,
                dynamic_cast<Mesh*>(mesh),
                dynamic_cast<Material*>(material),
                dynamic_cast<Light*>(light),
                nullptr
            });
        });
        return true;
    }

    virtual void foreachTriangle(const ProcessTriangleFunc& processTriangle) const override {
        for (const auto& primitive : primitives_) {
            if (!primitive.mesh) {
                continue;
            }
            primitive.mesh->foreachTriangle([&](int face, Vec3 p1, Vec3 p2, Vec3 p3) {
                processTriangle(
                    primitive.index,
                    face,
                    Vec3(primitive.transform * Vec4(p1, 1_f)),
                    Vec3(primitive.transform * Vec4(p2, 1_f)),
                    Vec3(primitive.transform * Vec4(p3, 1_f))
                );
            });
        }
    }

    // ------------------------------------------------------------------------

    virtual void build(const std::string& name) override {
        accel_ = comp::create<Accel>(name, this);
        if (!accel_) {
            return;
        }
        accel_->build(*this);
    }

    virtual std::optional<SurfacePoint> intersect(Ray ray, Float tmin, Float tmax) const override {
        const auto hit = accel_->intersect(ray, tmin, tmax);
        if (!hit) {
            return {};
        }
        const auto [t, uv, primitive, face] = *hit;
        const auto* mesh = primitives_.at(primitive).mesh;
        const auto p = mesh->surfacePoint(face, uv);
        return SurfacePoint(primitive, -1, p.p, p.n, p.t);
    }

    // ------------------------------------------------------------------------

    virtual bool isLight(const SurfacePoint& sp) const override {
        return primitives_.at(sp.primitive).light;
    }

    virtual bool isSpecular(const SurfacePoint& sp) const override {
        const auto& primitive = primitives_.at(sp.primitive);
        if (sp.endpoint) {
            if (primitive.light) {
                return primitive.light->isSpecular(sp);
            }
            else if (primitive.camera) {
                return primitive.camera->isSpecular(sp);
            }
            LM_UNREACHABLE_RETURN();
        }
        return primitive.material->isSpecular(sp);
    }

    // ------------------------------------------------------------------------

    virtual Ray primaryRay(Vec2 rp) const {
        return primitives_.at(camera_).camera->primaryRay(rp);
    }

    virtual std::optional<RaySample> sampleRay(Rng& rng, const SurfacePoint& sp, Vec3 wi) const override {
        return primitives_.at(sp.primitive).material->sampleRay(rng, sp, wi);
    }

    virtual std::optional<RaySample> samplePrimaryRay(Rng& rng, Vec4 window) const override {
        auto rs = primitives_.at(camera_).camera->samplePrimaryRay(rng, window);
        if (!rs) {
            return {};
        }
        return rs->asPrimitive(camera_).asEndpoint(true);
    }

	virtual std::optional<LightSample> sampleLight(Rng& rng, const SurfacePoint& sp) const override {
		
	}

	// ------------------------------------------------------------------------

	virtual Vec3 evalFs(const SurfacePoint& sp, Vec3 wi, Vec3 wo) const override {
		
	}

    virtual Vec3 evalContrbEndpoint(const SurfacePoint& sp, Vec3 wo) const override {
        const auto& primitive = primitives_.at(sp.primitive);
        if (!primitive.light) {
            return {};
        }
        return primitive.light->eval(sp, wo);
    }

    virtual std::optional<Vec3> reflectance(const SurfacePoint& sp) const override {
        const auto& primitive = primitives_.at(sp.primitive);
        if (!primitive.material) {
            return {};
        }
        return primitive.material->reflectance(sp);
    }

private:
    int camera_;    // Camera primitive index
    std::vector<Primitive> primitives_;
    Component::Ptr<Accel> accel_;
};

LM_COMP_REG_IMPL(Scene_, "scene::default");

LM_NAMESPACE_END(LM_NAMESPACE)