vec3 evaluateColor(in Ray ray, in HitInfo h) {
    if(shading == 0) {
        return h.m.kd;
    } else if(shading == 1) {
        vec3 position = vec3(ray.xd * ray.t + ray.xo, ray.yd * ray.t + ray.yo, ray.zd * ray.t + ray.zo);
        vec3 light_direction = normalize(light_position - position);
        return evalutateLambertian(h, vec3(1.0, 1.0, 1.0), light_direction);
    }
    return vec3(0.0, 1.0, 0.0);
}

vec3 evalutateLambertian(in HitInfo h, in vec3 light_color, in vec3 light_direction) {
    if(dot(h.normal, light_direction) < 0) {
        return vec3(0.0, 0.0, 0.0);
    }
    return h.m.kd * light_color * dot(h.normal, light_direction);
}