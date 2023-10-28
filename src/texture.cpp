#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
	size_t i = (size_t)(image.width * texCoord.x);
	size_t j = image.width * (size_t)(image.height * texCoord.y);
    return image.pixels[i + j];
}


// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
	float x = ((float)image.width * texCoord.x);
	float y = ((float)image.height * texCoord.y);

	int i = (int)(image.width * texCoord.x);
	int j = (int)(image.height * texCoord.y);
	int i_1 = glm::min(i + 1, image.width);
	int j_1 = glm::min(j + 1, image.height);

	float alpha = x - (float)i;
	float beta = y - (float)j;

	glm::vec3 top, bottom;

	top = (1.0f - alpha) * image.pixels[i + image.width * j] + (alpha) *  image.pixels[i_1 + image.width * j];
	bottom = (1.0f - alpha) * image.pixels[i + image.width * j_1] + (alpha) *  image.pixels[i_1 + image.width * j_1];

	return beta * bottom + (1.0f - beta) * top;
}