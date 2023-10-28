#pragma once
#include <filesystem>
#include <glad/glad.h>

std::string readFile(const std::filesystem::path path);

class ComputeShader {



public:

	ComputeShader(GLuint program);
	ComputeShader(ComputeShader&& shader);
	ComputeShader();
	~ComputeShader();

	void set_shader_file(const std::filesystem::path &path);
	void add_shader_header_file(const std::filesystem::path &path);
	void ad_shader_definition_file(const std::filesystem::path &path);

	void use();
	void compute(GLuint size_x, GLuint size_y, GLuint size_z = 1u);
	void compile_shader();
	GLuint m_program;

private:

	std::filesystem::path shader_path;
	std::vector<std::filesystem::path> shader_headers;
	std::vector<std::filesystem::path> shader_definitions;
};

