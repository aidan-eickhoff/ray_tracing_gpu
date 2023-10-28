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
