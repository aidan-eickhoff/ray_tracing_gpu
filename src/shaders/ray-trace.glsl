#version 430 core

precision highp float;

//include structs, uniforms, and buffers
#incl_def "common"
#incl_def "buffers"

//include forward function definitions
#include "intersect"
#include "shading"

const uint leafBit = 1 << 31;

//function definitions
void generateRay(out Ray rval);



void main() {
    vec3 color = vec3(0.0, 0.0, 0.0);
    ivec2 texelCoord = ivec2(gl_GlobalInvocationID.xy);
	
    Ray ray;
    generateRay(ray);
    HitInfo h;

    if(intersectRayBVH(ray, h)) {
        color = evaluateColor(ray, h);
    }

    vec4 value = vec4(color, 1.0);
    imageStore(imgOutput, texelCoord, value);
}

void generateRay(out Ray rval) {
    rval.xo = camera_position.x;
    rval.yo = camera_position.y;
    rval.zo = camera_position.z;

	float pix_x = ((float(gl_GlobalInvocationID.x) + .5)/(gl_NumWorkGroups.x)) * 2.0 - 1.0;
    float pix_y = ((float(gl_GlobalInvocationID.y) + .5)/(gl_NumWorkGroups.y)) * 2.0 - 1.0;

    vec3 camera_space_dir = normalize(vec3(pix_x * half_screen_space.x, -pix_y * half_screen_space.y, 1.0));
    vec3 real_dir = camera_rotation * camera_space_dir;

    rval.xd = real_dir.x;
    rval.yd = real_dir.y;    
    rval.zd = real_dir.z;    

    rval.t = 1e30;
}

