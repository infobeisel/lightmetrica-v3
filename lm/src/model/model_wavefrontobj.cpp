/*
    Lightmetrica - Copyright (c) 2018 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/model.h>
#include <lm/mesh.h>
#include <lm/material.h>
#include <lm/texture.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

// ----------------------------------------------------------------------------

// Surface geometry shared among meshes
struct OBJSurfaceGeometry {
    std::vector<Vec3> ps;   // Positions
    std::vector<Vec3> ns;   // Normals
    std::vector<Vec2> ts;   // Texture coordinates
};

// Face indices
struct OBJMeshFaceIndex {
    int p = -1;     // Index of position
    int t = -1;     // Index of texture coordinates
    int n = -1;     // Index of normal
};

// Texture parameters
struct MTLTextureParams {
    std::string name;   // Name
    std::string path;   // Texture path
};

// MLT material parameters
struct MTLMatParams {
    std::string name;   // Name
    int illum;          // Type
    Vec3 Kd;            // Diffuse reflectance
    Vec3 Ks;            // Specular reflectance
    Vec3 Ke;            // Luminance
    int mapKd;          // Texture index for Kd
    Float Ni;           // Index of refraction
    Float Ns;           // Specular exponent for phong shading
    Float an;           // Anisotropy
};

// Wavefront OBJ/MTL file parser
class WavefrontOBJParser {
private:
    // Material parameters
    std::vector<MTLMatParams> ms_;
    std::unordered_map<std::string, int> msmap_;
    // Texture parameters
    std::vector<MTLTextureParams> ts_;
    std::unordered_map<std::string, int> tsmap_;

public:
    // Callback functions
    using ProcessMeshFunc = std::function<
        std::optional<int>(const std::vector<OBJMeshFaceIndex>& fs, int materialIndex)>;
    using ProcessMaterialFunc = std::function<
        bool(const MTLMatParams& m)>;
    using ProcessTextureFunc = std::function<
        bool(const MTLTextureParams& path)>;

    // Parses .obj file
    bool parse(const std::string& p, OBJSurfaceGeometry& geo,
        const ProcessMeshFunc& processMesh,
        const ProcessMaterialFunc& processMaterial,
        const ProcessTextureFunc& processTexture)
    {
        LM_INFO("Loading OBJ file [path='{}']", p);
        char l[4096], name[256];
        std::ifstream f(p);
        if (!f) {
            LM_ERROR("Missing OBJ file [path='{}']", p);
            return false;
        }

        // Active face indices and material index
        int currMaterialIdx = -1;
        std::vector<OBJMeshFaceIndex> currfs;

        // Parse .obj file line by line
        while (f.getline(l, 4096)) {
            char *t = l;
            skipSpaces(t);
            if (command(t, "v", 1)) {
                geo.ps.emplace_back(nextVec3(t += 2));
            } else if (command(t, "vn", 2)) {
                geo.ns.emplace_back(nextVec3(t += 3));
            } else if (command(t, "vt", 2)) {
                geo.ts.emplace_back(nextVec3(t += 3));
            } else if (command(t, "f", 1)) {
                t += 2;
                if (ms_.empty()) {
                    // Process the case where MTL file is missing
                    ms_.push_back({ "default", -1, Vec3(1) });
                    if (!processMaterial(ms_.back())) {
                        return false;
                    }
                    currMaterialIdx = 0;
                }
                OBJMeshFaceIndex is[4];
                for (auto& i : is) {
                    i = parseIndices(geo, t);
                }
                currfs.insert(currfs.end(), {is[0], is[1], is[2]});
                if (is[3].p != -1) {
                    // Triangulate quad
                    currfs.insert(currfs.end(), {is[0], is[2], is[3]});
                }
            } else if (command(t, "usemtl", 6)) {
                t += 7;
                nextString(t, name);
                if (!currfs.empty()) {
                    // 'usemtl' indicates end of mesh groups
                    processMesh(currfs, currMaterialIdx);
                    currfs.clear();
                }
                currMaterialIdx = msmap_.at(name);
            } else if (command(t, "mtllib", 6)) {
                nextString(t += 7, name);
                if (!loadmtl((std::filesystem::path(p).remove_filename() / name).string(),
                    processMaterial, processTexture)) {
                    return false;
                }
            } else {
                continue;
            }
        }
        if (!currfs.empty()) {
            processMesh(currfs, currMaterialIdx);
        }

        return true;
    }

private:
    // Checks a character is space-like
    bool whitespace(char c) { return c == ' ' || c == '\t'; };

    // Checks the token is a command 
    bool command(char *&t, const char *c, int n) { return !strncmp(t, c, n) && whitespace(t[n]); }

    // Skips spaces
    void skipSpaces(char *&t) { t += strspn(t, " \t"); }

    // Skips spaces or /
    void skipSpacesOrComments(char *&t) { t += strcspn(t, "/ \t"); }

    // Parses floating point value
    Float nf(char *&t) {
        skipSpaces(t);
        Float v = Float(atof(t));
        skipSpacesOrComments(t);
        return v;
    }

    // Parses int value
    int nextInt(char *&t) {
        skipSpaces(t);
        int v = atoi(t);
        skipSpacesOrComments(t);
        return v;
    }

    // Parses 3d vector
    Vec3 nextVec3(char *&t) {
        Vec3 v;
        v.x = nf(t);
        v.y = nf(t);
        v.z = nf(t);
        return v;
    }

    // Parses vertex index. See specification of obj file for detail.
    int parseIndex(int i, int vn) { return i < 0 ? vn + i : i > 0 ? i - 1 : -1; }
    OBJMeshFaceIndex parseIndices(OBJSurfaceGeometry& geo, char *&t) {
        OBJMeshFaceIndex i;
        skipSpaces(t);
        i.p = parseIndex(atoi(t), int(geo.ps.size()));
        skipSpacesOrComments(t);
        if (t++[0] != '/') { return i; }
        i.t = parseIndex(atoi(t), int(geo.ts.size()));
        skipSpacesOrComments(t);
        if (t++[0] != '/') { return i; }
        i.n = parseIndex(atoi(t), int(geo.ns.size()));
        skipSpacesOrComments(t);
        return i;
    }

    // Parses a string
    template <int N>
    void nextString(char *&t, char (&name)[N]) { sscanf_s(t, "%s", name, N); };

    // Parses .mtl file
    bool loadmtl(std::string p, const ProcessMaterialFunc& processMaterial, const ProcessTextureFunc& processTexture)
    {
        LM_INFO("Loading MTL file [path='{}']", p);
        std::ifstream f(p);
        if (!f) {
            LM_ERROR("Missing MLT file [path='{}']", p);
            return false;
        }
        char l[4096], name[256];
        while (f.getline(l, 4096)) {
            auto *t = l;
            skipSpaces(t);
            if (command(t, "newmtl", 6)) {
                nextString(t += 7, name);
                msmap_[name] = int(ms_.size());
                ms_.emplace_back();
                ms_.back().name = name;
                continue;
            }
            if (ms_.empty()) {
                continue;
            }
            auto& m = ms_.back();
            if      (command(t, "Kd", 2))     { m.Kd = nextVec3(t += 3); }
            else if (command(t, "Ks", 2))     { m.Ks = nextVec3(t += 3); }
            else if (command(t, "Ni", 2))     { m.Ni = nf(t += 3); }
            else if (command(t, "Ns", 2))     { m.Ns = nf(t += 3); }
            else if (command(t, "aniso", 5))  { m.an = nf(t += 5); }
            else if (command(t, "Ke", 2))     { m.Ke = nextVec3(t += 3); }
            else if (command(t, "illum", 5))  { m.illum = nextInt(t += 6); }
            else if (command(t, "map_Kd", 6)) {
                nextString(t += 7, name);
                // Check if the texture is alread loaded
                auto it = tsmap_.find(name);
                if (it != tsmap_.end()) {
                    m.mapKd = it->second;
                }
                else {
                    // Register a new texture
                    m.mapKd = tsmap_[name] = int(ts_.size());
                    ts_.push_back({name, (std::filesystem::path(p).remove_filename()/name).string()});
                }
            }
        }
        // Let the user to process texture and materials
        for (const auto& t : ts_) {
            if (!processTexture(t)) {
                return false;
            }
        }
        for (const auto& m : ms_) {
            if (!processMaterial(m)) {
                return false;
            }
        }
        return true;
    }
};

// ----------------------------------------------------------------------------

class Mesh_WavefrontObj : public Mesh {
private:
    OBJSurfaceGeometry& geo_;
    std::vector<OBJMeshFaceIndex> fs_;

public:
    Mesh_WavefrontObj(OBJSurfaceGeometry& geo, const std::vector<OBJMeshFaceIndex>& fs)
        : geo_(geo), fs_(fs) {}
    
    virtual void foreachTriangle(const ProcessTriangleFunc& processTriangle) const override {
        for (int fi = 0; fi < int(fs_.size()); fi += 3) {
            const auto p1 = geo_.ps[fs_[fi].p];
            const auto p2 = geo_.ps[fs_[fi + 1].p];
            const auto p3 = geo_.ps[fs_[fi + 2].p];
            processTriangle(fi, p1, p2, p3);
        }
    }

    virtual Point surfacePoint(int face, Vec2 uv) const override {
        const auto i1 = fs_[face];
        const auto i2 = fs_[face + 1];
        const auto i3 = fs_[face + 2];
        return {
            math::mixBarycentric(
                geo_.ps[i1.p], geo_.ps[i2.p], geo_.ps[i3.p], uv),
            glm::normalize(math::mixBarycentric(
                geo_.ns[i1.n], geo_.ns[i2.n], geo_.ns[i3.n], uv)),
            math::mixBarycentric(
                geo_.ts[i1.t], geo_.ts[i2.t], geo_.ts[i3.t], uv)
        };
    }
};

// ----------------------------------------------------------------------------

class Material_WavefrontObj : public Material {
private:
    MTLMatParams objmat_;

public:
    Material_WavefrontObj(const MTLMatParams& m) : objmat_(m) {}
    
public:


};

// ----------------------------------------------------------------------------

class Model_WavefrontObj : public Model {
private:
    // Surface geometry
    OBJSurfaceGeometry geo_;
    // Underlying meshes. The instances are directly created instead of using comp::create.
    std::vector<Ptr<Mesh_WavefrontObj>> meshes_;
    // Underlying materials
    std::vector<Ptr<Material_WavefrontObj>> materials_;
    // Underlying textures
    std::vector<Ptr<Texture>> textures_;
    // Mesh group which assocites a mesh and a material
    std::vector<std::tuple<int, int>> groups_;

public:
    virtual bool construct(const Json& prop) override {
        WavefrontOBJParser parser;
        parser.parse(prop["path"], geo_,
            // Process mesh
            [&](const std::vector<OBJMeshFaceIndex>& fs, int materialIndex) -> std::optional<int> {
                auto mesh = comp::detail::createDirect<Mesh_WavefrontObj>(this, geo_, fs);
                groups_.push_back({int(meshes_.size()), materialIndex});
                meshes_.push_back(std::move(mesh));
                return true;
            },
            // Process material
            [&](const MTLMatParams& m) -> bool {
                auto mat = comp::detail::createDirect<Material_WavefrontObj>(this, m);
                materials_.push_back(std::move(mat));
                return true;
            },
            // Process texture
            [&](const MTLTextureParams& tex) -> bool {
                auto texture = comp::create<Texture>("texture::bitmap", this, {{"path", tex.path}});
                if (!texture) {
                    return false;
                }
                textures_.push_back(std::move(texture));
                return true;
            });
        return true;
    }
    
    virtual void createPrimitives(const CreatePrimitiveFunc& createPrimitive) const override {
        for (auto [mesh, material] : groups_) {
            createPrimitive(meshes_.at(mesh).get(), materials_.at(material).get());
        }
    }
};

LM_COMP_REG_IMPL(Model_WavefrontObj, "model::wavefrontobj");

LM_NAMESPACE_END(LM_NAMESPACE)
