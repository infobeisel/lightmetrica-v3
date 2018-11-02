/*
    Lightmetrica - Copyright (c) 2018 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/api.h>
#include <lm/assets.h>
#include <lm/scene.h>
#include <lm/renderer.h>
#include <lm/logger.h>
#include <lm/film.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

// ----------------------------------------------------------------------------

namespace {

/*
    \brief User API context.
    Manages all global states manipulated by user apis.
*/
class Context : public Component {
public:
    static Context& instance() {
        static Context instance;
        return instance;
    }
    
public:
    virtual Component* underlying(const std::string& name) const {
        const auto [s, r] = comp::splitFirst(name);
        if (s == "assets") {
            return comp::getCurrentOrUnderlying(r, assets_.get());
        }
        else if (s == "scene") {
            return comp::getCurrentOrUnderlying(r, scene_.get());
        }
        else if (s == "renderer") {
            return renderer_.get();
        }
        return nullptr;
    }

public:
    void init() {
        log::init();
        LM_INFO("Initializing Lightmetrica [version='{}']");
        if (init_) {
            LM_WARN("Lightmetrica is already initialized");
            return;
        }
        assets_ = comp::create<Assets>("assets::default", this);
        scene_  = comp::create<Scene>("scene::default", this);
        init_ = true;
    }

    void shutdown() {
        if (!checkInitialized()) { return; }
        renderer_.reset();
        scene_.reset();
        assets_.reset();
        log::shutdown();
    }

    void asset(const std::string& name, const std::string& implKey, const Json& prop) {
        if (!checkInitialized()) { return; }
        assets_->loadAsset(name, implKey, prop);
    }

    void primitive(Mat4 transform, const Json& prop) {
        if (!checkInitialized()) { return; }
        scene_->loadPrimitive(*assets_.get(), transform, prop);
    }

    void primitives(Mat4 transform, const std::string& modelName) {
        if (!checkInitialized()) { return; }
        scene_->loadPrimitives(*assets_.get(), transform, modelName);
    }

    void render(const std::string& rendererName, const std::string& accelName, const Json& prop) {
        if (!checkInitialized()) { return; }
        scene_->build(accelName);
        renderer_ = lm::comp::create<Renderer>(rendererName, this, prop);
        if (!renderer_) {
            LM_ERROR("Failed to render [renderer='{}',accel='{}']", rendererName, accelName);
            return;
        }
        LM_INFO("Starting render [name='{}']", rendererName);
        renderer_->render(*scene_.get());
    }

    void save(const std::string& filmName, const std::string& outpath) {
        if (!checkInitialized()) { return; }
        const auto* film = assets_->underlying<Film>(filmName);
        if (!film) {
            return;
        }
        if (!film->save(outpath)) {
            return;
        }
    }

private:
    bool checkInitialized() {
        if (!init_) {
            LM_ERROR("Lightmetrica is not initialized");
            return false;
        }
        return true;
    }

private:
    bool init_ = false;
    Ptr<Assets> assets_;
    Ptr<Scene> scene_;
    Ptr<Renderer> renderer_;
};

}

// ----------------------------------------------------------------------------

LM_PUBLIC_API void init() {
    Context::instance().init();
}

LM_PUBLIC_API void shutdown() {
    Context::instance().shutdown();
}

LM_PUBLIC_API void asset(const std::string& name, const std::string& implKey, const Json& prop) {
    Context::instance().asset(name, implKey, prop);
}

LM_PUBLIC_API void primitive(Mat4 transform, const Json& prop) {
    Context::instance().primitive(transform, prop);
}

LM_PUBLIC_API void primitives(Mat4 transform, const std::string& modelName) {
    Context::instance().primitives(transform, modelName);
}

LM_PUBLIC_API void render(const std::string& rendererName, const std::string& accelName, const Json& prop) {
    Context::instance().render(rendererName, accelName, prop);
}

LM_PUBLIC_API void save(const std::string& filmName, const std::string& outpath) {
    Context::instance().save(filmName, outpath);
}

LM_NAMESPACE_END(LM_NAMESPACE)
