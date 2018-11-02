/*
    Lightmetrica - Copyright (c) 2018 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/lm.h>

int main(int argc, char** argv) {
    // Command line arguments
    // ----------------------
    // ./003_obj
    //   scene_path output_path image_width image_height
    //   camera_pos_{x,y,z} camera_lookat_{x,y,z} vertical_fov
    char* a = argv[1];
    const std::string objPath(a++);
    const std::string outputPath(a++);
    const int w = atoi(a++);
    const int h = atoi(a++);
    const glm::vec3 cameraPosition(atof(a++), atof(a++), atof(a++));
    const glm::vec3 cameraLookat(atof(a++), atof(a++), atof(a++));
    const auto vfov = lm::Float(atof(a++));

    // Initialize the framework
    // ------------------------
    lm::init();

    // Define assets
    // -------------
    // Film for the rendered image
    lm::asset("film", "film::bitmap", {
        {"w", w},
        {"h", h}
    });

    // Pinhole camera
    lm::asset("camera1", "camera::pinhole", {
        {"position", lm::castToJson(cameraPosition)},
        {"center", lm::castToJson(cameraLookat)},
        {"up", {0,1,0}},
        {"vfov", vfov},
        {"aspect", (lm::Float)(w) / h}
    });

    // OBJ model
    lm::asset("obj1", "model::wavefrontobj", {
        {"path", objPath}
    });
    
    // Define scene primitives
    // -----------------------
    // Camera
    lm::primitive(lm::Mat4(1), {
        {"camera", "camera1"}
    });

    // Create primitives from model asset
    lm::primitives(lm::Mat4(1), "obj1");

    // Render an image
    // ---------------
    lm::render("renderer::raycast", "accel::sahbvh", {
        {"output", "film"},
        {"color", lm::castToJson(lm::Vec3(0))}
    });

    // Save rendered image
    lm::save("film", outputPath);

    // Finalize the framework
    // ----------------------
    lm::shutdown();

    return 0;
}