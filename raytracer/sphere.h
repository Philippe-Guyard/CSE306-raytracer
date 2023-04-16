#ifndef SPHERE_H
#define SPHERE_H

#include <optional>

#include "vector.h"
#include "geometry.h"

class Sphere : public Geometry {
private:
    const Vector3 center;
    double radius;
public:
    Sphere() = default;
    Sphere(const Vector3& center, 
           double radius, 
           const Vector3& color,
           bool is_mirror = false,
           bool is_transparent = false): center(center), 
                                         radius(radius),
                                         Geometry(color, is_mirror, is_transparent) {};

    double get_r() const { return radius; }
    const Vector3& get_center() { return center; }
    std::optional<Intersection> intersect(const Ray& ray) const override;
};

#endif