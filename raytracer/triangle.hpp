#pragma once 

#include "vector.h"
#include "geometry.h"

#include <vector>
#include <optional>
#include <iostream>

class Triangle {
private:
    // TODO: Huge memory loss duplicating normals and vertices
    Vector3 A, B, C;
    Vector3 n1, n2, n3;
    Vector3 N;
public:
    Triangle(const Vector3& a, const Vector3& b, const Vector3& c,
             const Vector3& n1, const Vector3& n2, const Vector3& n3): A(a), B(b), C(c), n1(n1), n2(n2), n3(n3) {
        N = (b - a).cross(c - a);
    };

    std::optional<Intersection> intersect(const Ray& ray, const Geometry *src) const {
        const Vector3& O = ray.get_origin();
        const Vector3& u = ray.get_direction();
        Vector3 p = (A - O).cross(u);
        // What do we do about 0 determinant?
        double det = u.dot(N);
        double beta = p.dot(C - A) / det;
        if (beta < 0 || beta > 1)
            return std::nullopt;

        double gamma = -p.dot(B - A) / det;
        if (gamma < 0 || gamma > 1)
            return std::nullopt;

        double alpha = 1 - beta - gamma;
        if (alpha < 0 || alpha > 1)
            return std::nullopt;

        double t = (A - O).dot(N) / det;
        if (t < 0)
            return std::nullopt;
            
        return Intersection(src, ray, O + t * u, alpha * n1 + beta * n2 + gamma * n3, t);
    }

    const Vector3& get_a() const {
        return A;
    }

    const Vector3& get_b() const {
        return B;
    }

    const Vector3& get_c() const {
        return C;
    }

    Triangle& operator=(Triangle other) {
        // NOTE: This makes me want to cry, but whatever...
        A = other.A;
        B = other.B;
        C = other.C;
        n1 = other.n1;
        n2 = other.n2;
        n3 = other.n3;
        N = other.N;
        return *this;
    }
};