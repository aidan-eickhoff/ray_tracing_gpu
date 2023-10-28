#include "bvh.h"
#include "config.h"
#include "sampler.h"
#include "modern_screen.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
#include "compute_shader.h"
#include "glm/gtc/type_ptr.hpp"
DISABLE_WARNINGS_PUSH()
#include <fmt/chrono.h>
#include <fmt/core.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <glm/gtx/quaternion.hpp>
#include <imgui/imgui.h>
#include <nativefiledialog/nfd.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <framework/imguizmo.h>
#include <framework/trackball.h>
#include <framework/variant_helper.h>
#include <framework/window.h>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <variant>
#include <iostream>

void updateSettingsBuffer(GLuint settings_buffer, Features features);


// This is the main application. The code in here does not need to be modified.
enum class ViewMode {
    Rasterization = 0,
    RayTracing = 1
};

int debugBVHLeafId = 0;
uint32_t debugRaySeed = 4; // Chosen by fair dice roll

int main(int argc, char** argv)
{
    Config config = {};
    if (argc > 1) {
        config = readConfigFile(argv[1]);
    } else {
        // Add a default camera if no config file is given.
        config.cameras.emplace_back(CameraConfig {});
    }
	config.features.enableAccelStructure = true;


    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
                << std::endl;

	//create a window
    Window window { "Final Project", config.windowSize, OpenGLVersion::GL2, true };

	//create a "screen" -- handles rendering of the calculated view stored in a 2D texture
	ModernScreen modern_screen("vertex.glsl", "fragment.glsl", config.windowSize);
	modern_screen.set_fov(config.cameras[0].fieldOfView);

	//create an object which handles our ray-tracing compute shader
	ComputeShader cs;
	cs.set_shader_file("ray-trace.glsl");
	cs.compile_shader();
	ModernScreen::gen_texture();


	//creates a trackball object which allows panning and zooming -- adapted from TU Delft professors :)
    Trackball camera { &window, glm::radians(config.cameras[0].fieldOfView), config.cameras[0].distanceFromLookAt };
	config.features.extra.dof.fov = glm::radians(config.cameras[0].fieldOfView);
    camera.setCamera(config.cameras[0].lookAt, glm::radians(config.cameras[0].rotation), config.cameras[0].distanceFromLookAt + 6);

	//select which scene to load
    SceneType sceneType { SceneType::Dragon };

	//load the scene
	Scene scene = loadScenePrebuilt(sceneType, config.dataPath);
    BVH bvh(scene, config.features);

	//allow window resize
	modern_screen.bind_resize_callback(window.get_window());

	//push the BVH Nodes to the GPU
	GLuint node_buffer;
	glGenBuffers(1, &node_buffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, node_buffer);
	glNamedBufferStorage(node_buffer, bvh.nodes().size_bytes(), bvh.nodes().data(), GL_MAP_READ_BIT);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, node_buffer);

	//push the primitives to the GPU
	GLuint primitive_buffer;
	glGenBuffers(1, &primitive_buffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, primitive_buffer);
	glNamedBufferStorage(primitive_buffer, bvh.primitives().size_bytes(), bvh.primitives().data(), GL_MAP_READ_BIT);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, primitive_buffer);


	//uniform variables which *shouldn't* change frequently
	unsigned int uniform_data[] = {static_cast<unsigned int>(bvh.primitives().size()), static_cast<unsigned int>(bvh.nodes().size()), static_cast<unsigned int>(scene.meshes.size())};
	GLuint uniform_buffer;
	glGenBuffers(1, &uniform_buffer);
	glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer);
	glNamedBufferData(uniform_buffer, sizeof(uniform_data), &(uniform_data[0]), GL_STATIC_READ);
	glBindBufferBase(GL_UNIFORM_BUFFER, 4, uniform_buffer);


	//push the materials to the GPU
	std::vector<Material> mesh_materials;
	for(const Mesh &m : scene.meshes) {
		mesh_materials.push_back(m.material);
	}

	GLuint mesh_material_buffer;
	glGenBuffers(1, &mesh_material_buffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, mesh_material_buffer);
	glNamedBufferStorage(mesh_material_buffer, mesh_materials.size() * sizeof(Material), mesh_materials.data(), GL_MAP_READ_BIT);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, mesh_material_buffer);


	//push the settings to the GPU
	config.features.enableShading = false;
	config.features.shadingModel = ShadingModel::Lambertian;
	config.features.enableNormalInterp = false;
	config.features.enableAccelStructure = true;

	unsigned int settings[] = {
		config.features.enableAccelStructure,
		static_cast<unsigned int>(config.features.enableShading) + static_cast<unsigned int>(config.features.shadingModel),
		config.features.enableShadows,
		config.features.enableTransparency,
		config.features.enableReflections,
		config.features.enableNormalInterp};

	GLuint settings_buffer;
	glGenBuffers(1, &settings_buffer);
	glBindBuffer(GL_UNIFORM_BUFFER, settings_buffer);
	glNamedBufferData(settings_buffer, sizeof(settings), &(settings[0]), GL_STATIC_READ);
	glBindBufferBase(GL_UNIFORM_BUFFER, 5, settings_buffer);


	//main loop
	while(!window.shouldClose()) {
		window.updateInput();

		using clock = std::chrono::high_resolution_clock;
		const auto start = clock::now();
		cs.use();

		//uniform variables which must be updated on a frame by frame basis
		glUniformMatrix3fv(1, 1, GL_FALSE, glm::value_ptr(glm::toMat3(camera.get_rotation_quat())));
		glUniform3fv(2, 1, glm::value_ptr(camera.position()));
		glUniform2fv(3, 1, glm::value_ptr(ModernScreen::get_screen_space()));
		try {
			glUniform3fv(6, 1, glm::value_ptr(std::get<0>(scene.lights[0]).position));
		} catch (const std::bad_variant_access& ex) {
			glm::vec3 location(1.0f, 1.0f, 1.0f);
			glUniform3fv(6, 1, glm::value_ptr(location));
		}

		//call the ray-tracing function
		cs.compute(ModernScreen::get_screen_size().x, ModernScreen::get_screen_size().y, 1);
		//wait until the ray-tracing function is complete
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		//draw the computed screen through rasterized texture mapping
		modern_screen.draw();

		//swap window buffers & time the frame
		glfwSwapBuffers(window.get_window());
		const auto end = clock::now();
		const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		fmt::print("Rendering took {} ms.\n", duration);

	}
    return 0;
}

void updateSettingsBuffer(GLuint settings_buffer, Features features) {
	unsigned int settings[] = {
		features.enableAccelStructure,
		features.enableShading,
		features.enableShadows,
		features.enableTransparency,
		features.enableReflections,
		features.enableNormalInterp};

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, settings_buffer);
	glNamedBufferData(settings_buffer, sizeof(settings), &(settings[0]), GL_STATIC_READ);
}