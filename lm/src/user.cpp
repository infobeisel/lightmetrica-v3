/*
    Lightmetrica - Copyright (c) 2018 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/user.h>
#include <lm/assets.h>
#include <lm/scene.h>
#include <lm/renderer.h>
#include <lm/logger.h>
#include <lm/film.h>
#include <lm/parallel.h>
#include <lm/exception.h>

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
    void init(const Json& prop) {
        exception::init();
        log::init();
        LM_INFO("Initializing Lightmetrica [version='{}']");
        parallel::init("parallel::openmp", prop);
        assets_ = comp::create<Assets>("assets::default", this);
        scene_  = comp::create<Scene>("scene::default", this);
    }

    void shutdown() {
        renderer_.reset();
        scene_.reset();
        assets_.reset();
        parallel::shutdown();
        log::shutdown();
        exception::shutdown();
    }

    void asset(const std::string& name, const std::string& implKey, const Json& prop) {
        assets_->loadAsset(name, implKey, prop);
    }

    void primitive(Mat4 transform, const Json& prop) {
        scene_->loadPrimitive(*assets_.get(), transform, prop);
    }

    void primitives(Mat4 transform, const std::string& modelName) {
        scene_->loadPrimitives(*assets_.get(), transform, modelName);
    }

    void render(const std::string& rendererName, const std::string& accelName, const Json& prop) {
        scene_->build(accelName);
        renderer_ = lm::comp::create<Renderer>(rendererName, this, prop);
        if (!renderer_) {
            LM_ERROR("Failed to render [renderer='{}',accel='{}']", rendererName, accelName);
            return;
        }
        LM_INFO("Starting render [name='{}']", rendererName);
        LM_INDENT();
        renderer_->render(*scene_.get());
    }

    void save(const std::string& filmName, const std::string& outpath) {
        const auto* film = assets_->underlying<Film>(filmName);
        if (!film) {
            return;
        }
        if (!film->save(outpath)) {
            return;
        }
    }

    FilmBuffer buffer(const std::string& filmName) {
        auto* film = assets_->underlying<Film>(filmName);
        if (!film) {
            return {};
        }
        return film->buffer();
    }

private:
    Ptr<Assets> assets_;
    Ptr<Scene> scene_;
    Ptr<Renderer> renderer_;
};

}

// ----------------------------------------------------------------------------

LM_PUBLIC_API void init(const Json& prop) {
    Context::instance().init(prop);
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

LM_PUBLIC_API FilmBuffer buffer(const std::string& filmName) {
    return Context::instance().buffer(filmName);
}

LM_NAMESPACE_END(LM_NAMESPACE)
