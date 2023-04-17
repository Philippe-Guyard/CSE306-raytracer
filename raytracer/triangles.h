#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vector.h"
#include "geometry.h"

#include <vector>
#include <optional>
#include <iostream>

class Triangle {
private:
    Vector3 A, e1, e2;
    Vector3 n1, n2, n3;
    Vector3 N;
public:
    struct RayInter {
        double t;
        Vector3 normal;
    };

    Triangle(const Vector3& a, const Vector3& b, const Vector3& c,
             const Vector3& n1, const Vector3& n2, const Vector3& n3): A(a), n1(n1), n2(n2), n3(n3) {
        e1 = b - a;
        e2 = c - a;
        N = e1.cross(e2);
    };

    std::optional<RayInter> intersect(const Ray& ray) const {
        const Vector3& O = ray.get_origin();
        const Vector3& u = ray.get_direction();
        Vector3 p = (A - O).cross(u);
        // What do we do about 0 determinant?
        double det = u.dot(N);
        double beta = p.dot(e2) / det;
        if (beta < 0 || beta > 1)
            return std::nullopt;

        double gamma = -p.dot(e1) / det;
        if (gamma < 0 || gamma > 1)
            return std::nullopt;

        double alpha = 1 - beta - gamma;
        if (alpha < 0 || alpha > 1)
            return std::nullopt;

        // std::cout << "beta: " << beta << " gamma: " << gamma << std::endl;
        double t = (A - O).dot(N) / det;
        if (t < 0)
            return std::nullopt;
            
        return RayInter{t, (1 - beta - gamma) * n1 + beta * n2 + gamma * n3};
    }    
};

class TriangleMesh : public Geometry {
private:
    std::vector<Triangle> triangles;
public:
    TriangleMesh(const std::vector<Triangle>& triangles, const Vector3& color, 
                bool is_mirror = false, 
                bool is_transparent = false): Geometry(color, is_mirror, is_transparent) {
        this->triangles = triangles;
    }

    TriangleMesh(std::vector<Triangle>&& triangles, const Vector3& color, 
                bool is_mirror = false, 
                bool is_transparent = false): Geometry(color, is_mirror, is_transparent) {
        this->triangles = std::move(triangles);
    }

    std::optional<Intersection> intersect(const Ray& r) const {
        std::optional<Triangle::RayInter> pre_inter;
        for (size_t i = 0; i < triangles.size(); i++) {
            auto temp = triangles[i].intersect(r);
            if (temp.has_value()) {
                if (!pre_inter.has_value() || temp.value().t < pre_inter.value().t) {
                    pre_inter.emplace(std::move(temp.value()));
                }
            }
        }
        if (!pre_inter.has_value())
            return std::nullopt;
        
        double t = pre_inter.value().t;
        Vector3 normal = pre_inter.value().normal;
        return Intersection(this, r, r.get_origin() + t * r.get_direction(), normal, t);
    }
};


#endif 