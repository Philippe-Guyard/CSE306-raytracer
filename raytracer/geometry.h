#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <optional>

#include "vector.h"
#include "ray.h"

class Geometry;

struct Intersection {
    const Geometry *object;
    const Ray& source;
    Vector3 point, normal;
    double t;

    // Useless to use temporaries here since Vector3 has no heap allocation.
    // See this: https://www.youtube.com/watch?v=ehMg6zvXuMY
    // And this: https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
    Intersection(const Geometry *object, 
                 const Ray& source, 
                 const Vector3& point, 
                 const Vector3& normal, 
                 double t): object(object), source(source), point(point), normal(normal), t(t) {};

    Intersection(Intersection&& other) noexcept: source(other.source) {
        object = other.object;
        point = std::move(other.point);
        normal = std::move(other.normal);
        t = other.t;
        other.object = nullptr;
    };
};

class Geometry {
protected:
    Vector3 color;
    bool m_is_mirror, m_is_transparent;
public:
    Geometry(const Vector3& color, 
            bool is_mirror = false, 
            bool is_transparent = false): color(color), m_is_mirror(is_mirror), m_is_transparent(is_transparent) {};
    virtual ~Geometry() {};

    virtual std::optional<Intersection> intersect(const Ray& ray) const = 0;
    const Vector3& get_color() const { return color; }

    bool is_mirror() const { return m_is_mirror; }
    bool is_transparent() const { return m_is_transparent; }
};

#endif 