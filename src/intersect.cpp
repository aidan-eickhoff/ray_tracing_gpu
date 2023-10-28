#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
	return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
	return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
	Plane plane;
	return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo & hitInfo)
{
	glm::vec3 edge1 = v1 - v0, edge2 = v2 - v0;
	glm::vec3 h = glm::cross(ray.direction, edge2 );
	float a = glm::dot(edge1, h);
	if (a > -1e-5f && a < 1e-5f) {
		return false; // ray parallel to triangle
	}
	float f = 1 / a;
	glm::vec3 s = ray.origin - v0;
	float u = f * glm::dot(s, h);
	if (u < 0 || u > 1) {
		return false;
	}
	glm::vec3 q = glm::cross(s, edge1);
	float v = f * glm::dot(ray.direction, q);
	if (v < 0.0f || u + v > 1.0f) {
		return false;
	}
	float t = f * glm::dot(edge2, q);
	if (t > 0.0f && t < ray.t) {
		ray.t = t;
		return true;
	}
	return false;
}


/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo & hitInfo)
{
	float A = glm::dot(ray.direction, ray.direction);
	float B = 2 * glm::dot(ray.direction, ray.origin - sphere.center);
	float C = glm::dot(ray.origin - sphere.center, ray.origin - sphere.center) - sphere.radius * sphere.radius;
	float det = B * B - 4.0f * A * C;
	if(det < 0) {
		return false;
	}
	float sol1 = (-B + std::sqrt(det)) / (2 * A);
	float sol2 = (-B - std::sqrt(det)) / (2 * A);

	if(sol1 <= 0 && sol2 <= 0) {
		return false;
	} else if(sol1 > 0 && sol2 > 0) {
		if(ray.t > (sol1 <= sol2 ? sol1 : sol2)) {
			ray.t = sol1 <= sol2 ? sol1 : sol2;
		}
	} else {
		if(ray.t > (sol1 <= 0 ? sol2 : sol1)) {
			ray.t = sol1 <= 0 ? sol2 : sol1;
		}
	}
	return true;

}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
	float x_inv = 1.0f / ray.direction.x;
	float y_inv = 1.0f / ray.direction.y;
	float z_inv = 1.0f / ray.direction.z;


	float t_x_lower = (box.lower.x - ray.origin.x) * x_inv;
	float t_x_upper = (box.upper.x - ray.origin.x) * x_inv;
	float t_y_lower = (box.lower.y - ray.origin.y) * y_inv;
	float t_y_upper = (box.upper.y - ray.origin.y) * y_inv;
	float t_z_lower = (box.lower.z - ray.origin.z) * z_inv;
	float t_z_upper = (box.upper.z - ray.origin.z) * z_inv;

	float t_in_all = glm::max(glm::min(t_x_lower, t_x_upper), glm::max(glm::min(t_y_lower, t_y_upper), glm::min(t_z_lower, t_z_upper)));
	float t_out_one = glm::min(glm::max(t_x_lower, t_x_upper), glm::max(glm::max(t_y_lower, t_y_upper), glm::max(t_z_lower, t_z_upper)));

	if(t_in_all > t_out_one || t_out_one <= 0) {
		return false;
	} 
	return true;
}