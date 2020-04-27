/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/core.h>
#include <lm/texture.h>
#pragma warning(push)
#pragma warning(disable:4244) // possible loss of data
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#pragma warning(pop)

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

std::string sanitize_directory_separator(std::string p) {
    std::replace(p.begin(), p.end(), '\\', '/');
    return p;
}

/*
\rst
.. function:: texture::bitmap

    Bitmap texture.

    :param str path: Path to texture.
    :param bool flip: Flip loaded texture if true.
\endrst
*/
class Texture_Bitmap final : public Texture {
private:
    int w_;     // Width of the image
    int h_;     // Height of the image
    int c_;     // Number of components
    std::vector<float> data_;

public:
    LM_SERIALIZE_IMPL(ar) {
        ar(w_, h_, c_, data_);
    }

public:
    virtual TextureSize size() const override {
        return { w_, h_ };
    }

    virtual void construct(const Json& prop) override {
        // Image path
        const std::string path = sanitize_directory_separator(json::value<std::string>(prop, "path"));
        LM_INFO("Loading texture [path='{}']", fs::path(path).filename().string());

        // Load as HDR image
        // LDR image is internally converted to HDR
        const bool flip = json::value<bool>(prop, "flip", true);
        stbi_set_flip_vertically_on_load(flip);
        float* data = stbi_loadf(path.c_str(), &w_, &h_, &c_, 0);
        if (data == nullptr) {
            LM_ERROR("Failed to load image: {} [path='{}']", stbi_failure_reason(), path);
            LM_THROW_EXCEPTION_DEFAULT(Error::IOError);
        }
        // Allocate and copy the data
        data_.assign(data, data + (w_*h_*c_));
        stbi_image_free(data);
    }

    virtual Vec3 eval(Vec2 t) const override {
        const auto u = t.x - floor(t.x);
        const auto v = t.y - floor(t.y);
        const int x = std::clamp(int(u * w_), 0, w_ - 1);
        const int y = std::clamp(int(v * h_), 0, h_ - 1);
        const int i = w_ * y + x;
        return Vec3(data_[c_*i], data_[c_*i+1], data_[c_*i+2]);
    }

    virtual Vec3 eval_by_pixel_coords(int x, int y) const override {
        const int i = w_*y + x;
        return Vec3(data_[c_*i], data_[c_*i+1], data_[c_*i+2]);
    }

    virtual Float eval_alpha(Vec2 t) const override {
        const auto u = t.x - floor(t.x);
        const auto v = t.y - floor(t.y);
        const int x = std::clamp(int(u * w_), 0, w_ - 1);
        const int y = std::clamp(int(v * h_), 0, h_ - 1);
        const int i = w_ * y + x;
        return data_[c_*i+3];
    }

    virtual bool has_alpha() const override {
        return c_ == 4;
    }

    virtual TextureBuffer buffer() override {
        return { w_, h_, c_, data_.data() };
    }
};

LM_COMP_REG_IMPL(Texture_Bitmap, "texture::bitmap");

LM_NAMESPACE_END(LM_NAMESPACE)
