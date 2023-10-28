
bool intersectRayBVH(inout Ray ray, inout HitInfo h);

float intersectRayAABB(in Ray ray, in AxisAlignedBox aabb);

bool intersectRayTriangle(in uint i, inout Ray ray, inout HitInfo hitInfo);

void fillOutHitInfo(in Ray ray, in Primitive p, inout HitInfo h);
