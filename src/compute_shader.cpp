#include <glad/glad.h>
#include <iostream>
#include <fstream>
#include <fmt/format.h>
#include "compute_shader.h"

static constexpr GLuint invalid = 0xFFFFFFFF;
static void check_error(GLuint program, bool shader);

ComputeShader::ComputeShader(GLuint program)
	: m_program(program) {}

ComputeShader::ComputeShader()
	: m_program(invalid) {}

ComputeShader::ComputeShader(ComputeShader&& other) {
	m_program = other.m_program;
	other.m_program = invalid;
}

ComputeShader::~ComputeShader() {
	if (m_program != invalid)
		glDeleteProgram(m_program);
}

void ComputeShader::set_shader_file(const std::filesystem::path &path) {
	shader_path = std::filesystem::path(SHADER_DIR).append(path.c_str());
}

void ComputeShader::compute(GLuint size_x, GLuint size_y, GLuint size_z) {
	if(m_program == invalid) {
		std::cout << "Invalid Shader" << std::endl;
	}
	glDispatchCompute(size_x, size_y, size_z);
}

void ComputeShader::use() {
	glUseProgram(m_program);
}

static std::string read_file(const std::filesystem::path& path) {
	std::ifstream file(path, std::ios::binary);
	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

void ComputeShader::compile_shader() {
	if (!std::filesystem::exists(shader_path)) {
		throw std::runtime_error::runtime_error(fmt::format("File {} does not exist", shader_path.string().c_str()));
	}

	//delete old program in case of update
	if(m_program != invalid) {
		glDeleteProgram(m_program);
	}

	std::string computeCode;
	std::ifstream cShaderFile;
	// ensure ifstream objects can throw exceptions:
	cShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		// open files
		cShaderFile.open(shader_path);

		std::stringstream cShaderStream;
		// read file's buffer contents into streams
		cShaderStream << cShaderFile.rdbuf();
		// close file handlers
		cShaderFile.close();
		// convert stream into string
		computeCode = cShaderStream.str();
	}
	catch (std::ifstream::failure& e)
	{
		std::cout << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << std::endl;
	}
	const char* cShaderCode = computeCode.c_str();
	// 2. compile shaders
	GLuint compute;
	// compute shader
	compute = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(compute, 1, &cShaderCode, NULL);
	glCompileShader(compute);
	check_error(compute, true);

	// shader Program
	m_program = glCreateProgram();
	glAttachShader(m_program, compute);
	glLinkProgram(m_program);
	check_error(m_program, false);
	// delete the shaders as they're linked into our program now and no longer necessary
	glDeleteShader(compute);

}

static void check_error(GLuint program, bool shader) {
	GLint success;
	if(shader) {
		glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

	} else {
		glGetProgramiv(program, GL_LINK_STATUS, &success);
	}
	if(!success) {
		GLchar infoLog[1024];
		if(shader) {
			glGetShaderInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		} else {
			glGetProgramInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		}
	}
}