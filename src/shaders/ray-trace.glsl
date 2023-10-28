#version 430 core

struct Ray {
	float xo, yo, zo;
    float xd, yd, zd;
    float t;
    float padding;
};

struct Material {
    vec3 kd; // Diffuse color.
	vec3 ks;
	float shininess;
	float transparency;
};

struct HitInfo {
    vec3 barycentric_coords;
    vec3 normal;
    Material m;
};

struct AxisAlignedBox {
    float x_min, y_min, z_min;
    float x_max, y_max, z_max;
};

struct BVHNode {
    AxisAlignedBox aabb;
    uint data1, data2;
};

struct Vertex {
    float pos_x, pos_y, pos_z;
    float nx, ny, nz;
    float tx, ty;
};

struct Primitive {
    uint meshID;
    Vertex v0, v1, v2;
};

layout (local_size_x = 1, local_size_y = 1, local_size_z = 1) in;
layout(location = 0) uniform mat3 camera_rotation;
layout(location = 1) uniform vec3 camera_position;
layout(location = 2) uniform vec2 half_screen_space;
layout(location = 3) uniform vec3 light_position;


//texture to store data in
layout(rgba32f, binding = 0) uniform image2D imgOutput;

//BVH Nodes
layout(std430, binding = 1) buffer node_buffer
{
    BVHNode nodes[];
};

//Scene primitives
layout(std430, binding = 2) buffer primitive_buffer
{
    Primitive primitives[];
};

//Mesh colors/reflectivity/transparency
layout(std430, binding = 3) buffer meshColorBuffer
{
    Material mesh_materials[];
};

//BVH Node sizes
layout(packed, binding = 4) uniform UBuffer
{
    uint node_buffer_size;
    uint primitive_buffer_size;
    uint mesh_materials_size;
};

//rendering settings
layout(binding = 5) uniform options
{
    uint acceleration;
    uint shading; //0 - off, 1 - lambertian, 2 - phong, 3 - blinn-phong
    uint shadows;
    uint transparency;
    uint reflection;
    uint interpolation;
    // bool multisampling;
};


const uint leafBit = 1 << 31;

//function definitions
void generateRay(out Ray rval);
bool intersectRayBVH(inout Ray ray, inout HitInfo h);
float intersectRayAABB(in Ray ray, in AxisAlignedBox aabb);
bool intersectRayTriangle(in vec3 v0, in vec3 v1, in vec3 v2, inout Ray ray, inout HitInfo hitInfo);
vec3 evaluateColor(in Ray ray, in HitInfo hitInfo);
void fillOutHitInfo(in Ray ray, in Primitive p, inout HitInfo h);

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

	float pix_x = float(gl_GlobalInvocationID.x)/(gl_NumWorkGroups.x) * 2.0 - 1.0;
    float pix_y = float(gl_GlobalInvocationID.y)/(gl_NumWorkGroups.y) * 2.0 - 1.0;

    vec3 camera_space_dir = normalize(vec3(pix_x * half_screen_space.x, -pix_y * half_screen_space.y, 1.0));
    vec3 real_dir = camera_rotation * camera_space_dir;

    // vec3 real_dir = normalize(-pix_x * camera_information[0] + pix_y + camera_information[1] + camera_information[2]);

    rval.xd = real_dir.x;
    rval.yd = real_dir.y;    
    rval.zd = real_dir.z;    

    rval.t = 1e30;
}


vec3 evaluateColor(in Ray ray, in HitInfo h) {
    if(shading == 0) {
        return h.m.kd;
    } else if(shading == 1) {
        vec3 position = vec3(ray.xd * ray.t + ray.xo, ray.yd * ray.t + ray.yo, ray.zd * ray.t + ray.zo);
        vec3 light_direction = normalize(light_position - position);
        if(dot(light_direction, h.normal) < 0) {
            return vec3(0.0, 0.0, 0.0);
        }
        return h.m.kd * dot(light_direction, h.normal);
    }
    return vec3(0.0, 1.0, 0.0);
}

void fillOutHitInfo(in Ray ray, in Primitive p, inout HitInfo h) {
    h.m = mesh_materials[p.meshID];
    if(interpolation == 1) {
        h.normal.x = p.v0.nx * h.barycentric_coords.x + p.v1.nx * h.barycentric_coords.y + p.v2.nx * h.barycentric_coords.z;
        h.normal.y = p.v0.ny * h.barycentric_coords.x + p.v1.ny * h.barycentric_coords.y + p.v2.ny * h.barycentric_coords.z;
        h.normal.z = p.v0.nz * h.barycentric_coords.x + p.v1.nz * h.barycentric_coords.y + p.v2.nz * h.barycentric_coords.z;
    } else {
        h.normal.x = (p.v0.nx + p.v1.nx + p.v2.nx) / 3.0;
        h.normal.y = (p.v0.ny + p.v1.ny + p.v2.ny) / 3.0;
        h.normal.z = (p.v0.nz + p.v1.nz + p.v2.nz) / 3.0;
    }
    h.normal = normalize(h.normal);
    if (dot(vec3(ray.xd, ray.yd, ray.zd), h.normal) > 0.0) {
		h.normal = -h.normal;
	}
}

bool intersectRayBVH(inout Ray ray, inout HitInfo h) {
    bool rval = false;
    vec3 v0, v1, v2;
    
    if(acceleration == 1) {
        uint curr = 0;
        uint stack[64];
        uint stack_ptr = 0;

        while(true) {
            if((nodes[curr].data1 & leafBit) != 0) {
                uint offset = nodes[curr].data1 & (~leafBit);
                for(uint i = offset; i < offset + nodes[curr].data2; i++) {
                    v0 = vec3(primitives[i].v0.pos_x, primitives[i].v0.pos_y, primitives[i].v0.pos_z);
                    v1 = vec3(primitives[i].v1.pos_x, primitives[i].v1.pos_y, primitives[i].v1.pos_z);
                    v2 = vec3(primitives[i].v2.pos_x, primitives[i].v2.pos_y, primitives[i].v2.pos_z);

                    if(intersectRayTriangle(v0, v1, v2, ray, h)) {
                        rval = true;
                        fillOutHitInfo(ray, primitives[i], h);
                    }
                }
                
                if(stack_ptr == 0) {
                    break;
                } else {
                    curr = stack[--stack_ptr];
                }
                continue;
            }

            uint child_1 = nodes[curr].data1;
            uint child_2 = nodes[curr].data2;
            float t_1 = intersectRayAABB(ray, nodes[child_1].aabb);
            float t_2 = intersectRayAABB(ray, nodes[child_2].aabb);

            if(t_1 > t_2) {
                uint temp = child_2;
                child_2 = child_1;
                child_1 = temp;
                float temp_2 = t_2;
                t_1 = t_2;
                t_2 = temp_2;
            }

            if(t_1 == 1e30) {
                if(stack_ptr == 0) {
                    break;
                } else {
                    curr = stack[--stack_ptr];
                }
            } else {
                curr = child_1;
                if(t_2 != 1e30) {
                    stack[stack_ptr++] = child_2;
                }
            }
        }
    } else {
        for(int i = 0; i < primitive_buffer_size; i++) {
            v0 = vec3(primitives[i].v0.pos_x, primitives[i].v0.pos_y, primitives[i].v0.pos_z);
            v1 = vec3(primitives[i].v1.pos_x, primitives[i].v1.pos_y, primitives[i].v1.pos_z);
            v2 = vec3(primitives[i].v2.pos_x, primitives[i].v2.pos_y, primitives[i].v2.pos_z);

            if(intersectRayTriangle(v0, v1, v2, ray, h)) {
                rval = true;
                fillOutHitInfo(ray, primitives[i], h);
            }
        }
    }
    return rval;
}

 
float intersectRayAABB(in Ray ray, in AxisAlignedBox aabb) {
    float x_inv = 1.0 / ray.xd;
	float y_inv = 1.0 / ray.yd;
	float z_inv = 1.0 / ray.zd;

	float t_x_lower = (aabb.x_min - ray.xo) * x_inv;
	float t_x_upper = (aabb.x_max - ray.xo) * x_inv;

	float t_y_lower = (aabb.y_min - ray.yo) * y_inv;
	float t_y_upper = (aabb.y_max - ray.yo) * y_inv;

	float t_z_lower = (aabb.z_min - ray.zo) * z_inv;
	float t_z_upper = (aabb.z_max - ray.zo) * z_inv;

	float t_in_all = max(min(t_x_lower, t_x_upper), max(min(t_y_lower, t_y_upper), min(t_z_lower, t_z_upper)));
	float t_out_one = min(max(t_x_lower, t_x_upper), min(max(t_y_lower, t_y_upper), max(t_z_lower, t_z_upper)));

	if(t_in_all > t_out_one || t_out_one <= 0 || t_in_all > ray.t) {
		return 1e30;
	} 
	return t_in_all;
}

bool intersectRayTriangle(in vec3 v0, in vec3 v1, in vec3 v2, inout Ray ray, inout HitInfo hitInfo) {
    vec3 ray_d = vec3(ray.xd, ray.yd, ray.zd);
    vec3 ray_o = vec3(ray.xo, ray.yo, ray.zo);
	vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;

	vec3 h = cross(ray_d, edge2);
	float a = dot(edge1, h);
	if (a > -1e-5 && a < 1e-5) {
		return false; // ray parallel to triangle
	}
	float f = 1 / a;
	vec3 s = ray_o - v0;
	float u = f * dot(s, h);
	if (u < 0.0 || u > 1.0) {
		return false;
	}
	vec3 q = cross(s, edge1);
	float v = f * dot(ray_d, q);
	if (v < 0.0 || u + v > 1.0) {
		return false;
	}
	float t = f * dot(edge2, q);
	if (t > 0.0 && t < ray.t) {
		ray.t = t;
        hitInfo.barycentric_coords.x = 1.0 - u - v;
        hitInfo.barycentric_coords.y = u;
        hitInfo.barycentric_coords.z = v;
        return true;
	}
    return false;
}