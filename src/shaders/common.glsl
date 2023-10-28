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
