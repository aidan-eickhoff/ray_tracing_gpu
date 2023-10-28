#pragma once
#include <glad/glad.h>
#include <filesystem>

#include "GLFW/glfw3.h"
#include "glm/vec2.hpp"




class ModernScreen {

public:
	ModernScreen();
	ModernScreen(const std::filesystem::path& path_vertex, const std::filesystem::path& path_fragment, glm::ivec2 size);
	~ModernScreen();
	void draw();
	void set_fov(float new_fovy);
	void bind_resize_callback(GLFWwindow* window);


	static float get_aspect_ratio();
	static GLuint gen_texture();
	static glm::vec2 get_screen_space();
	static glm::ivec2 get_screen_size();
	static void set_window(GLFWwindow* w);

private:
	inline static GLuint m_texture;
	inline static GLuint m_program;
	inline static GLFWwindow* m_window;
	inline static GLuint quadVAO = 0, quadVBO;

	inline static float fovy, half_ssw, half_ssh;
	inline static glm::ivec2 screen_size;

	static void update_screen_space();
	static void window_resize_callback(GLFWwindow* window, int width, int height);

};
