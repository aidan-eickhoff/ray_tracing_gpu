#include "bvh.h"
#include "intersect.h"
#include "scene.h"
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <queue>
#include <utility>

using Primitive = BVHInterface::Primitive;
using std::cout;
using std::endl;


// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
	// Store start of bvh build for timing
	using clock = std::chrono::high_resolution_clock;
	const auto start = clock::now();
#endif

	// Count the total nr. of triangles in the scene
	size_t numTriangles = 0;
	for (const auto& mesh : scene.meshes)
		numTriangles += mesh.triangles.size();

	// Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
	std::vector<Primitive> primitives;
	primitives.reserve(numTriangles);
	for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
		const auto& mesh = scene.meshes[meshID];
		for (const auto& triangle : mesh.triangles) {
			primitives.push_back(Primitive {
				.meshID = meshID,
				.v0 = mesh.vertices[triangle.x],
				.v1 = mesh.vertices[triangle.y],
				.v2 = mesh.vertices[triangle.z] });
		}
	}

	// Tell underlying vectors how large they should approximately be
	m_primitives.reserve(numTriangles);
	m_nodes.reserve(numTriangles + 1);

	// Recursively build BVH structure; this is where your implementation comes in
	m_nodes.emplace_back(); // Create root node
	m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
	buildRecursive(scene, features, primitives, RootIndex);

	// Fill in boilerplate data
	buildNumLevels();
	buildNumLeaves();

#ifndef NDEBUG
	// Output end of bvh build for timing
	const auto end = clock::now();
	std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
	const auto idx = static_cast<uint32_t>(m_nodes.size());
	m_nodes.emplace_back();
	return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
	glm::vec3 l = glm::vec3(std::min({primitive.v0.position.x, primitive.v1.position.x, primitive.v2.position.x}), 
							std::min({primitive.v0.position.y, primitive.v1.position.y, primitive.v2.position.y}),
							std::min({primitive.v0.position.z, primitive.v1.position.z, primitive.v2.position.z}));
	glm::vec3 u = glm::vec3(std::max({primitive.v0.position.x, primitive.v1.position.x, primitive.v2.position.x}), 
							std::max({primitive.v0.position.y, primitive.v1.position.y, primitive.v2.position.y}),
							std::max({primitive.v0.position.z, primitive.v1.position.z, primitive.v2.position.z}));
	return { .lower = l, .upper = u-l };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
	if(primitives.size() == 0) {
		return {};
	}
	float x_min = primitives[0].v0.position.x;
	float x_max = primitives[0].v0.position.x;
	float y_min = primitives[0].v0.position.y;
	float y_max = primitives[0].v0.position.y;
	float z_min = primitives[0].v0.position.z;
	float z_max = primitives[0].v0.position.z;

	for(const auto &p : primitives) {
		for(const auto &v : {p.v0, p.v1, p.v2}) {
			x_min = x_min < v.position.x ? x_min : v.position.x;
			x_max = x_max > v.position.x ? x_max : v.position.x;
			y_min = y_min < v.position.y ? y_min : v.position.y;
			y_max = y_max > v.position.y ? y_max : v.position.y;
			z_min = z_min < v.position.z ? z_min : v.position.z;
			z_max = z_max > v.position.z ? z_max : v.position.z;
		}

	}

	return { .lower = glm::vec3(x_min, y_min, z_min), .upper = glm::vec3(x_max, y_max , z_max) };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
	return (primitive.v0.position + primitive.v1.position + primitive.v2.position) / 3.0f;
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
	float x_len = aabb.upper.x - aabb.lower.x;
	float y_len = aabb.upper.y - aabb.lower.y;
	float z_len = aabb.upper.z - aabb.lower.z;
	if(x_len >= y_len && x_len >= z_len) {
		return 0;
	}
	if(y_len >= x_len && y_len >= z_len) {
		return 1;
	}
	if(z_len > x_len && z_len > y_len) {
		return 2;
	}
	return 0;
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element 
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
	using Primitive = BVHInterface::Primitive;

	std::sort(primitives.begin(), primitives.end(), 
			  [axis](Primitive a, Primitive b) {
		return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
	});

	return (primitives.size() + 1) / 2; 
}

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


// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
	Node node;
	// TODO fill in the leaf's data; refer to `bvh_interface.h` for details

	node.data[0] = m_primitives.size() | node.LeafBit;
	node.data[1] = primitives.size();
	node.aabb = aabb;

	// Copy the current set of primitives to the back of the primitives vector
	std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

	return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
	Node node;
	node.data = {leftChildIndex, rightChildIndex};
	node.aabb = aabb;
	return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
	// WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
	// because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.
	// Compute the AABB of the current node.
	AxisAlignedBox aabb = computeSpanAABB(primitives);

	if(primitives.size() <= LeafSize) {
		m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
	} else {
		uint32_t next_node_1 = nextNodeIdx();
		uint32_t median;
		if(features.extra.enableBvhSahBinning) {
			median = splitPrimitivesBySAHBin(aabb, computeAABBLongestAxis(aabb), primitives);
		} else {
			median = splitPrimitivesByMedian(aabb, computeAABBLongestAxis(aabb), primitives);
		}
		buildRecursive(scene, features, primitives.subspan(0, median), next_node_1);
		uint32_t next_node_2 = nextNodeIdx();
		buildRecursive(scene, features, primitives.subspan(median, primitives.size() - median), next_node_2);
		m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, next_node_1, next_node_2);
	}


	// Just configure the current node as a giant leaf for now
	//m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
}

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
	std::queue<Node> queue;
	std::queue<uint32_t> lvl_q;

	queue.push(m_nodes[0]);
	lvl_q.push(1);

	Node curr;
	uint32_t curr_lvl;

	m_numLevels = 0;

	while(!queue.empty()) {
		curr = queue.front();
		curr_lvl = lvl_q.front();
		queue.pop();
		lvl_q.pop();

		if(curr_lvl >= m_numLevels) {
			m_numLevels = curr_lvl;
		}

		if(!curr.isLeaf()) {
			queue.push(m_nodes[curr.leftChild()]);
			queue.push(m_nodes[curr.rightChild()]);
			lvl_q.push(curr_lvl + 1);
			lvl_q.push(curr_lvl + 1);
		}
	}
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
	std::queue<Node> queue;

	queue.push(m_nodes[0]);

	Node curr;

	m_numLeaves = 0;

	while(!queue.empty()) {
		curr = queue.front();
		queue.pop();

		if(!curr.isLeaf()) {
			queue.push(m_nodes[curr.leftChild()]);
			queue.push(m_nodes[curr.rightChild()]);
		} else {
			m_numLeaves++;
		}
	}
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
	// Example showing how to draw an AABB as a (white) wireframe box.
	// Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
	// colors, transparencies, etc.
	std::queue<Node> queue;
	std::queue<int> lvl_q;

	queue.push(m_nodes[0]);
	lvl_q.push(0);

	Node curr;
	int curr_lvl = 0;

	while(!queue.empty()) {
		curr = queue.front();
		curr_lvl = lvl_q.front();
		queue.pop();
		lvl_q.pop();

		if(curr_lvl == level) {
			drawAABB(curr.aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 1.0f);
		} else if(!curr.isLeaf()) {
			queue.push(m_nodes[curr.leftChild()]);
			queue.push(m_nodes[curr.rightChild()]);
			lvl_q.push(curr_lvl + 1);
			lvl_q.push(curr_lvl + 1);
		}
	}
}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
	// Example showing how to draw an AABB as a (white) wireframe box.
	// Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
	// colors, transparencies, etc.
	std::queue<Node> queue;

	queue.push(m_nodes[0]);

	Node curr;
	int currLeaf = 1;

	while(!queue.empty()) {
		curr = queue.front();
		queue.pop();

		if(!curr.isLeaf()) {
			queue.push(m_nodes[curr.leftChild()]);
			queue.push(m_nodes[curr.rightChild()]);
		} else if(currLeaf < leafIndex) {
			currLeaf++;
		} else {
			HitInfo h;
			drawAABB(curr.aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 1.0f);
			break;
		}
	}
}