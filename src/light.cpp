#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()
#include <framework/opengl_includes.h>

#include <iostream>

using std::cout;
using std::endl;


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    position = sample * light.endpoint1 + (1.0f - sample) * light.endpoint0;
    color = sample * light.color1 + (1.0f - sample) * light.color0;
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
	glm::vec3 col_edge01 =  sample.x * light.color1 + (1.0f - sample.x) * light.color0;
	glm::vec3 col_edge23 =  sample.x * light.color3 + (1.0f - sample.x) * light.color2;

	glm::vec3 pos_edge01 = sample.x * (light.v0 + light.edge01) + (1.0f - sample.x) * light.v0;
	glm::vec3 pos_edge23 = sample.x * (light.v0 + light.edge01 + light.edge02) + (1.0f - sample.x) * (light.v0 + light.edge02);

    position = sample.y * pos_edge23 + (1.0f - sample.y) * pos_edge01;
    color = sample.y * col_edge23 + (1.0f - sample.y) * col_edge01;
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
		HitInfo h;
		glm::vec3 origin = ray.direction * (ray.t - 1e-4f) + ray.origin;
		Ray outRay = {.origin=origin, .direction=glm::normalize(lightPosition - origin)};
        if(state.bvh.intersect(state, outRay, h)) {
			return glm::length(lightPosition - origin) <= glm::length(outRay.direction * outRay.t);
		}
        return true;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 origin = ray.origin + ray.direction * (ray.t - 1e-4f);
	Ray r = {.origin = origin, .direction = glm::normalize(lightPosition - origin)};
	HitInfo h;
	return helper_transparency(state, lightPosition, lightColor, r, h);
}

glm::vec3 helper_transparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, Ray& ray, HitInfo& hitInfo) {
	if(state.bvh.intersect(state, ray, hitInfo)) {
		if(glm::length(ray.t * ray.direction) > glm::length(lightPosition - ray.origin)) {
			return lightColor;
		} else if(hitInfo.material.transparency == 1) {
			return glm::vec3(0, 0, 0);

		} else {
			glm::vec3 origin = ray.origin + ray.direction * (ray.t + 1e-4f);
			glm::vec3 direction = glm::normalize(lightPosition - origin);
			Ray r = {.origin = origin, .direction =direction, .inv_d = 1.0f / direction};
			HitInfo h;
			return sampleMaterialKd(state, hitInfo) * (1 - hitInfo.material.transparency) * helper_transparency(state, lightPosition, lightColor, r, h);
		}
	}
	return lightColor;

}



// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;
	glm::vec3 l_color = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);
    return computeShading(state, v, l, l_color, hitInfo);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model

	glm::vec3 p = ray.origin + ray.t * ray.direction;
	glm::vec3 pos;
	glm::vec3 color;
	glm::vec3 v = -ray.direction;
	glm::vec3 final_color(0.0f);

	for(uint32_t i = 0u; i < numSamples; i++) {
		sampleSegmentLight(state.sampler.next_1d(), light, pos, color);
		color = visibilityOfLightSample(state, pos, color, ray, hitInfo);
		final_color += computeShading(state, v, glm::normalize(pos - p), color, hitInfo);
	}
    return final_color / (float)numSamples;
}


void debugSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
	// TODO: implement this function; repeat numSamples times:
	// - sample the segment light
	// - test the sample's visibility
	// - then evaluate the phong model
	Ray r = ray;
	HitInfo h = hitInfo;
	state.bvh.intersect(state, r, h);
	glm::vec3 p = ray.origin + r.t * ray.direction;
	glm::vec3 pos;
	glm::vec3 color;
	glm::vec3 v = -ray.direction;

	glm::vec3 final_color(0.0f);

	for(uint32_t i = 0u; i < numSamples; i++) {
		sampleSegmentLight(state.sampler.next_1d(), light, pos, color);
		color = visibilityOfLightSample(state, p, color, ray, hitInfo);
		final_color += computeShading(state, v, glm::normalize(pos - p), color, hitInfo);
		glColor3f(1.0f, 1.0f, 1.0f);
		glBegin(GL_LINES);
			glVertex3fv(glm::value_ptr(pos));
			glVertex3fv(glm::value_ptr(p));
		glEnd();
	}
}



// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
	glm::vec3 p = ray.origin + ray.t * ray.direction;
	glm::vec3 pos;
	glm::vec3 color;
	glm::vec3 v = -ray.direction;
	glm::vec3 final_color(0.0f);

	for(uint32_t i = 0u; i < numSamples; i++) {
		sampleParallelogramLight(state.sampler.next_2d(), light, pos, color);
		color = visibilityOfLightSample(state, pos, color, ray, hitInfo);
		final_color += computeShading(state, v, glm::normalize(pos - p), color, hitInfo);
	}
	return final_color / (float)numSamples;
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}