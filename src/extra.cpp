#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "draw.h"
#include "bvh_interface.h"
#include "intersect.h"
#include <framework/trackball.h>
#include <iostream>
#include <omp.h>


//techniques of how to implement DOF taken from https://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/drt.pdf
//couldn't find any explanation in the book

using std::endl;
using std::cout;
using Primitive = BVH::Primitive;


// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }
	// Either directly render the image, or pass through to extra.h methods

	const float screen_space_height = std::tan(features.extra.dof.fov / 2.0f);
	const float screen_space_width = screen_space_height * static_cast<float>(screen.resolution().x) / static_cast<float>(screen.resolution().y);
	glm::vec2 screenSpaceSize(screen_space_width, screen_space_height);
	const float invAperture = 0.25f / features.extra.dof.aperture;
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
	
	for (int y = 0; y < screen.resolution().y; y++) {
		for (int x = 0; x != screen.resolution().x; x++) {
			// Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
			// Note; we seed the sampler for consistent behavior across frames
			RenderState state = {
				.scene = scene,
				.features = features,
				.bvh = bvh,
				.sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
			};
		
			auto rays = generateDOFRays(state, camera, { x, y }, screen.resolution(), screenSpaceSize, invAperture);
			auto L = renderRays(state, rays);
			screen.setPixel(x, y, L);
		}
	}

	// Pass through to extra.h for post processing
	if (features.extra.enableBloomEffect) {
		postprocessImageWithBloom(scene, features, camera, screen);
	}
}


std::vector<Ray> generateDOFRays(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution, glm::vec2 screenSpaceSize, float invAperture) {
	glm::vec2 position = (glm::vec2(pixel) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f;
	//ripped directly from Trackball.cpp
	const glm::vec3 camera_direction = glm::quat(camera.rotationEulerAngles()) * glm::normalize(glm::vec3(-position.x * screenSpaceSize.x, position.y * screenSpaceSize.y, 1.0f));
	const float offset_ratio = 1.0f / state.features.extra.dof.distance;


	std::vector<Ray> rval;
	for(uint32_t i = 0; i < state.features.extra.dof.numSamples; i++) {
		//aperture
		glm::vec3 offset((state.sampler.next_2d() * 2.0f - 1.0f) * invAperture , 0.0f);
		glm::vec3 coord_offset = offset.x * camera.up() - offset.y * camera.left();
		Ray ray;
		//offset the rays starting position by a random amount on the 2d plane perpendicular to the viewing angle
		ray.origin = camera.position() + coord_offset;
		//then you must alter the rays direction in the oppposite direction based on how far you want the focused area to be
		//if focused area is N units away, alter ray direction by -1 * offset / N
		ray.direction = camera_direction - coord_offset * offset_ratio;
		ray.t = std::numeric_limits<float>::max();
		rval.push_back(ray);
	}
	return rval;
}


// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
	for (int y = 0; y < screen.resolution().y; y++) {
		for (int x = 0; x != screen.resolution().x; x++) {
			// Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
			// Note; we seed the sampler for consistenct behavior across frames
			RenderState state = {
				.scene = scene,
				.features = features,
				.bvh = bvh,
				.sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
			};
			auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
			auto L = renderRays(state, rays);
			screen.setPixel(x, y, L);
		}
	}

}


AxisAlignedBox transformAABB(const AxisAlignedBox& aabb, const glm::mat4& transformation) {
	AxisAlignedBox rval;

	rval.upper = glm::vec3(glm::vec4(aabb.upper, 1.0f) * transformation);
	rval.lower = glm::vec3(glm::vec4(aabb.lower, 1.0f) * transformation);

	return rval;
}

BVHInterface::Primitive transformPrimitive(const BVHInterface::Primitive& primitive, const glm::mat4& transformation) {
	BVHInterface::Primitive rval{};

	glm::mat3 rotations = glm::mat3(transformation);

	rval.v0 = {.position = glm::vec3(glm::vec4(primitive.v0.position, 1.0f) * transformation), .normal = primitive.v0.normal * rotations, .texCoord= primitive.v0.texCoord};
	rval.v1 = {.position = glm::vec3(glm::vec4(primitive.v1.position, 1.0f) * transformation), .normal = primitive.v0.normal * rotations, .texCoord= primitive.v1.texCoord};
	rval.v2 = {.position = glm::vec3(glm::vec4(primitive.v2.position, 1.0f) * transformation), .normal = primitive.v0.normal * rotations, .texCoord= primitive.v2.texCoord};
	rval.meshID = primitive.meshID;

	return rval;
}

bool intersectRayWithBVHTime(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo) {

	// Relevant data in the constructed BVH
	std::span<const BVHInterface::Node> nodes = bvh.nodes();
	std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

	const float time = state.time;
	const glm::mat4 transformation = glm::rotate(glm::identity<glm::mat4>(), glm::radians(10.0f * time), glm::vec3(0.0f, 1.0f, 0.0f));

	bool is_hit = false;



	if (state.features.enableAccelStructure) {
		std::vector<size_t> stack;

		stack.push_back(0);
		BVHInterface::Node curr;

		while(!stack.empty()) {
			curr = nodes[stack.back()];
			stack.pop_back();

			if(!curr.isLeaf()) {
				if(intersect_ray_shape(transformAABB(nodes[curr.leftChild()].aabb, transformation), ray)) {
					stack.push_back(curr.leftChild());
				}
				if(intersect_ray_shape(transformAABB(nodes[curr.rightChild()].aabb, transformation), ray)) {
					stack.push_back(curr.rightChild());
				}
			} else {
				for(size_t i = curr.primitiveOffset(); i < curr.primitiveCount() + curr.primitiveOffset(); i++) {
					const auto& t = [transformation](glm::vec3 a) {
						return glm::vec3(glm::vec4(a, 1.0f) * transformation);
					};
					if (intersectRayWithTriangle(t(primitives[i].v0.position), t(primitives[i].v1.position), t(primitives[i].v2.position), ray, hitInfo)) {
						updateHitInfo(state, transformPrimitive(primitives[i], transformation), ray, hitInfo);
						is_hit = true;
					}
				}
			}
		}
	} else {
		const auto& t = [transformation](glm::vec3 a) {
			return glm::vec3(glm::vec4(a, 1.0f) * transformation);
		};
		// Naive implementation; simply iterates over all primitives
		for (const auto& prim : primitives) {
			const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
			if (intersectRayWithTriangle(t(v0.position), t(v1.position), t(v2.position), ray, hitInfo)) {
				updateHitInfo(state, transformPrimitive(prim, transformation), ray, hitInfo);
				is_hit = true;
			}
		}
	}


	for (const auto& sphere : state.scene.spheres) {
		glm::mat4 rotation = glm::rotate(glm::identity<glm::mat4>(), glm::radians(45.0f * time), glm::vec3(1.0f, 0.0f, 0.0f));
		glm::vec4 s_center = glm::vec4(sphere.center, 1.0f) * rotation;
		Sphere new_sphere {.center=glm::vec3(s_center), .radius=sphere.radius, .material=sphere.material};
		if(intersectRayWithShape(new_sphere, ray, hitInfo)) {
			hitInfo.material = sphere.material;
			is_hit = true;
			hitInfo.normal = glm::normalize(ray.origin + ray.t * ray.direction  - new_sphere.center);
		}
	}

	return is_hit;

}



// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

	glm::vec2 size = image.resolution();
	size_t w = size.x, h = size.y;


	const size_t filter_size = 15u;
	const size_t offset = (filter_size + 1u) / 2u;
	float sum = static_cast<float>(1 << filter_size);
	std::vector<float> filter;

	for(uint32_t i = 0u; i <= filter_size; i++) {
		filter.push_back(calculate_nCr(filter_size, i) / sum);
	}

	long pos;
	Screen h_bloom(size);

	for(size_t i = 0u; i < image.pixels().size(); i++) {
		for(long j = 0; j <= filter_size; j++) {
			pos = i + j;
			if((pos) / w != i / w) {
				continue;
			}
			if(glm::dot(image.pixels()[pos], image.pixels()[pos]) < 0.81f) {
				continue;
			}
			h_bloom.pixels()[i + offset] += image.pixels()[pos] * filter[j];
		}
	}

	Screen v_bloom(size);
	const size_t v_offset = offset * size.y;

	for(long i = 0; i < h_bloom.pixels().size(); i++) {
		for(long j = 0; j <= filter_size; j++) {
			pos = i + j * w;
			if (pos < 0 || pos >= h_bloom.pixels().size() || i + v_offset >= v_bloom.pixels().size()) {
				continue;
			}
			if(glm::dot(h_bloom.pixels()[pos], h_bloom.pixels()[pos]) == 0.0f) {
				continue;
			}
			v_bloom.pixels()[i + v_offset] += h_bloom.pixels()[pos] * filter[j];
		}
	}

	for(uint32_t i = 0; i < v_bloom.pixels().size(); i++) {
		image.pixels()[i] +=  v_bloom.pixels()[i];
	}
}

float calculate_nCr(uint32_t n, uint32_t r) {
	float rval = 1.0f;
	for(uint32_t i = 0; i < r; i++) {
		rval *= (float)(n - i) / (float)(r - i);
	}
	return rval;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    uint32_t numSamples = state.features.extra.numGlossySamples;
	const float main_ray_weight = 1.0f / ((float)numSamples + 0.2f);
	const float glossy_ray_weight = (1.0f - main_ray_weight) / (float)numSamples;
	const float shinynessAmount = hitInfo.material.shininess / 128.0f;

	Ray r = generateReflectionRay(ray, hitInfo);

	hitColor += main_ray_weight * renderRay(state, r, rayDepth + 1);


	glm::vec3 perpendicular;
	if (ray.direction.z != 0) {
		perpendicular = glm::normalize(glm::vec3(1.0f, 1.0f, -(ray.direction.x + ray.direction.y) / ray.direction.z));
	} else {
		perpendicular = { 0, 0, 1 };
	}

	glm::vec3 second_perpendicular = glm::normalize(glm::cross(ray.direction, perpendicular));

	for(uint32_t i = 0u; i < numSamples; i++) {
		float theta = state.sampler.next_1d() * 2.0f * glm::pi<float>();
		glm::vec2 basis_translations(std::cos(theta), std::sin(theta));
		basis_translations *= shinynessAmount;
		glm::vec3 translation = basis_translations.x * perpendicular + basis_translations.y * second_perpendicular;

		Ray glossy_ray = {r.origin, r.direction + translation};
		hitColor += renderRay(state, glossy_ray, rayDepth + 1) * glossy_ray_weight;
	}

}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
	std::sort(primitives.begin(), primitives.end(), 
			  [axis](Primitive a, Primitive b) {
		return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
	});

	size_t split_num = 0;
	float min_value = std::numeric_limits<float>::max();

	size_t inc_amount = std::max(primitives.size() / 40u, 1ull);

	for(size_t i = 0; i < primitives.size(); i+=inc_amount) {
		float val = surface_area_of_primitives(primitives.subspan(0, i)) * i + surface_area_of_primitives(primitives.subspan(i, primitives.size() - i)) * (primitives.size() - i);

		if(val < min_value) {
			min_value = val;
			split_num = i;
		}
	}
	return split_num;
}

float surface_area_of_primitives(std::span<BVHInterface::Primitive> primitives) {
	AxisAlignedBox box = computeSpanAABB(primitives);
	glm::vec3 size = box.upper - box.lower;
	return 2.0f * (size.x * size.y + size.x * size.z + size.y * size.z);
}
