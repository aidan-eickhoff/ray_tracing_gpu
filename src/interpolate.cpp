#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
	glm::vec3 x = v1 - v0;
	glm::vec3 y = v2 - v0;
	glm::vec3 z = p - v0;

	float area = 0.5f * glm::length(glm::cross(x, y));	

	if(area == 0.0f) {
		return glm::vec3(0);
	}

	glm::vec3 area_side1 = glm::cross(x, z);
	glm::vec3 area_side2 = glm::cross(z, y);


	float gamma = (glm::length(area_side1) / 2.0f) / area;
	float beta = (glm::length(area_side2) / 2.0f) / area;

    return glm::vec3(1.0f - beta - gamma, beta, gamma);
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    return n0 * bc.x + n1 * bc.y + n2 * bc.z;
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    return t0 * bc.x + t1 * bc.y + t2 * bc.z;
}
