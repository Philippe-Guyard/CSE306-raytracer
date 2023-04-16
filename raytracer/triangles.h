#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vector.h"
#include "geometry.h"

#include <vector>
#include <optional>
#include <iostream>

class Triangle : public Geometry {
private:
    Vector3 a, e1, e2;
    Vector3 normal;
public:
    Triangle(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& color): a(a), Geometry(color) {
        e1 = b - a;
        e2 = c - a;
        normal = e1.cross(e2).normalized();
    };

    std::optional<Intersection> intersect(const Ray& ray) const {
        Vector3 p = (a - ray.get_origin()).cross(ray.get_direction());
        // What do we do about 0 determinant?
        double det = ray.get_direction().dot(normal);
        double beta = p.dot(e2) / det;
        if (beta < 0 || beta > 1)
            return std::nullopt;
        double gamma = p.dot(e1) / det;
        if (gamma < 0 || gamma > 1)
            return std::nullopt;

        // std::cout << "beta: " << beta << " gamma: " << gamma << std::endl;
        double t = (a - ray.get_origin()).dot(normal) / det;
        return Intersection(this, ray, ray.get_origin() + t * ray.get_direction(), normal, t);
    }    
};

class TriangleMesh : public Geometry {
private:
    std::vector<Triangle> triangles;
public:
    TriangleMesh(const std::vector<Triangle>& triangles) {
        this->triangles = triangles;
    }

    TriangleMesh(std::vector<Triangle>&& triangles) {
        this->triangles = std::move(triangles);
    }

    std::optional<Intersection> intersect(const Ray& r) const {
        std::optional<Intersection> intersection;
        for (const Triangle& triangle : triangles) {
            std::optional<Intersection> temp = triangle.intersect(r);
            if (temp.has_value()) {
                if (!intersection.has_value() || temp.value().t < intersection.value().t) {
                    intersection.emplace(std::move(temp.value));
                }
            }
        }

        return intersection;
    }
};


#endif 