#include "sphere.h"
#include "geometry.h"

#include <optional>

std::optional<Intersection> Sphere::intersect(const Ray& ray) const {
    double temp = ray.get_direction().dot(ray.get_origin() - center);
    double delta = temp * temp - ((ray.get_origin() - center).norm2() - radius * radius);
    if (delta >= 0) {
        double t1 = -temp - sqrt(delta);
        double t2 = -temp + sqrt(delta);
        if (t2 < 0)
            return std::nullopt;

        double t = t1 >= 0 ? t1 : t2;
        Vector3 intersection_point = ray.get_origin() + t * ray.get_direction();
        Vector3 normal = (intersection_point - center).normalized();
        return Intersection(this, ray, intersection_point, normal, t);
    }

    return std::nullopt;
}