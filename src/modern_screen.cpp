#include "modern_screen.h"
#include <filesystem>
#include <fstream>
#include <iostream>

#include "GLFW/glfw3.h"


static constexpr GLuint invalid = 0xFFFFFFFF;

static void check_error(GLuint program, bool shader) {
	GLint success;
	if(shader) {
		glGetShaderiv(program, GL_COMPILE_STATUS, &success);

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

ModernScreen::ModernScreen() {
	m_program = invalid;
	m_texture = invalid;
}

//loads the vertex shader and fragment shader from a file path, as well as sets the size of the screen
ModernScreen::ModernScreen(const std::filesystem::path& path_vertex,const std::filesystem::path& path_fragment, glm::ivec2 size) {
	// 1. retrieve the vertex/fragment source code from filePath
	std::string vertexCode;
	std::string fragmentCode;
	std::ifstream vShaderFile;
	std::ifstream fShaderFile;
	// ensure ifstream objects can throw exceptions:
	vShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
	fShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
	try 
	{
		// open files
		vShaderFile.open( std::filesystem::path(SHADER_DIR).append("rendering").append(path_vertex.c_str()));
		fShaderFile.open( std::filesystem::path(SHADER_DIR).append("rendering").append(path_fragment.c_str()));
		std::stringstream vShaderStream, fShaderStream;
		// read file's buffer contents into streams
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();		
		// close file handlers
		vShaderFile.close();
		fShaderFile.close();
		// convert stream into string
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();			
	}
	catch (std::ifstream::failure& e)
	{
		std::cout << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << std::endl;
	}
	const char* vShaderCode = vertexCode.c_str();
	const char * fShaderCode = fragmentCode.c_str();
	// 2. compile shaders
	GLuint vertex, fragment;
	// vertex shader
	vertex = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex, 1, &vShaderCode, NULL);
	glCompileShader(vertex);
	check_error(vertex, true);

	// fragment Shader
	fragment = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment, 1, &fShaderCode, NULL);
	glCompileShader(fragment);
	check_error(fragment, true);
	// shader Program
	m_program = glCreateProgram();
	glAttachShader(m_program, vertex);
	glAttachShader(m_program, fragment);
	glLinkProgram(m_program);
	check_error(m_program, false);
	// delete the shaders as they're linked into our program now and no longer necessary
	glDeleteShader(vertex);
	glDeleteShader(fragment);

	screen_size = size;
	m_texture = invalid;
}

GLuint ModernScreen::gen_texture() {
	if(m_texture != invalid) {
		glDeleteTextures(1, &m_texture);
	}

	glGenTextures(1, &m_texture);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, ModernScreen::get_screen_size().x, ModernScreen::get_screen_size().y, 0, GL_RGBA, 
				 GL_FLOAT, NULL);

	glBindImageTexture(0, m_texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

	glUniform1i(glGetUniformLocation(m_program, "tex"), 0); 

	return m_texture;
}

ModernScreen::~ModernScreen() {
	glDeleteProgram(m_program);
	glDeleteTextures(1, &m_texture);
}

//draws a 2x2 cube 1 unit away from the camera, displays the texture saved in GL_TEXTURE0
void ModernScreen::draw() {
	glClearColor(1.0f, 1.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(m_program);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_texture);
	if (quadVAO == 0)
	{
		float quadVertices[] = {
			// positions        // texture Coords
			-1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
			-1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
			1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
			1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
		};
		// setup plane VAO
		glGenVertexArrays(1, &quadVAO);
		glGenBuffers(1, &quadVBO);
		glBindVertexArray(quadVAO);
		glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
	}
	glBindVertexArray(quadVAO);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	glBindVertexArray(0);
}

//returns the aspect ratio of the screen
float ModernScreen::get_aspect_ratio() {
	return static_cast<float>(screen_size.x) / static_cast<float>(screen_size.y);
}

//returns the screen size
glm::ivec2 ModernScreen::get_screen_size() {
	return screen_size;
}

//called upon window resize
void ModernScreen::window_resize_callback(GLFWwindow* window, int width, int height) {
	screen_size = glm::ivec2(width, height);
	update_screen_space();
	glViewport(0, 0, width, height);
	gen_texture();
}

//set the FOV of the screen
void ModernScreen::set_fov(float new_fovy) {
	fovy = new_fovy;
	update_screen_space();
}

//update the screen space when either the window size or FOV has been updated
void ModernScreen::update_screen_space() {
	half_ssh = std::tan(fovy / 2.0f);
	half_ssw = static_cast<float>(screen_size.x) * half_ssh / static_cast<float>(screen_size.y);
}

void ModernScreen::bind_resize_callback(GLFWwindow* window) {
	glfwSetWindowSizeCallback(window, window_resize_callback);
}

glm::vec2 ModernScreen::get_screen_space() {
	return glm::vec2(half_ssw, half_ssh);
}


void ModernScreen::set_window(GLFWwindow* w) {
	m_window = w;
}
