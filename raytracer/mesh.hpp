#pragma once 

#include <vector>
#include <optional>
#include <iostream>

#include "vector.h"
#include "geometry.h"
#include "triangle.hpp"
#include "bvh.hpp"

class TriangleMesh : public Geometry {
private:
    std::vector<Vector3> vertices, normals;
    std::vector<Triangle> triangles;

    Texture *texture;
    BVH::BVH *bvh;
public:
    // TODO: This 
    // TriangleMesh(const std::vector<Triangle>& triangles, const Vector3& color, 
    //             bool is_mirror = false, 
    //             bool is_transparent = false): Geometry(color, is_mirror, is_transparent) {
    //     this->triangles = std::vector<Triangle>(triangles);
    // }

    TriangleMesh(std::vector<Vector3>&& vertices, std::vector<Vector3>&& normals, std::vector<Triangle>&& triangles, const Vector3& color,
                Texture *texture, bool is_mirror = false, bool is_transparent = false): Geometry(color, is_mirror, is_transparent) {
        this->texture = texture;
        this->vertices = std::move(vertices);
        this->triangles = std::move(triangles);
        this->bvh = new BVH::BVH(this, this->triangles.begin(), this->triangles.end());
    }

    std::optional<Intersection> intersect(const Ray& r) const {
        // std::optional<Intersection> closest = std::nullopt;
        // for (const Triangle& t : this->triangles) {
        //     std::optional<Intersection> i = t.intersect(r, this);
        //     if (i && (!closest || i.value().t < closest.value().t))
        //         closest.emplace(std::move(i.value()));
        // }
        // return closest;
        return this->bvh->intersects(r);
    }

    ~TriangleMesh() {
        delete this->bvh;
        delete this->texture;
    }
};
