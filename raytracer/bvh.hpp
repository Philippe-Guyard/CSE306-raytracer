#pragma once 

#include <vector>
#include <limits>
#include <algorithm>
#include <list>

#include "vector.h"
#include "ray.h"
#include "triangle.hpp"
#include "geometry.h"

namespace BVH {
    using vector_index = std::vector<Vector3>::iterator;
    using triangle_index = std::vector<Triangle>::iterator;

    struct Bbox final {
        Vector3 min;
        Vector3 max;

        Bbox() = default;
        Bbox(const Vector3& min, const Vector3& max): min(min), max(max) {};
        Bbox(triangle_index begin, triangle_index end) {
            double xmin = std::numeric_limits<double>::infinity();
            double ymin = std::numeric_limits<double>::infinity();
            double zmin = std::numeric_limits<double>::infinity();

            double xmax = -std::numeric_limits<double>::infinity();
            double ymax = -std::numeric_limits<double>::infinity();
            double zmax = -std::numeric_limits<double>::infinity();

            for(auto it = begin; it != end; it++) {
                auto vertices = {it->get_a(), it->get_b(), it->get_c()};
                for(const auto& vertex: vertices) {
                    xmin = std::min({xmin, vertex.x});
                    ymin = std::min({ymin, vertex.y});
                    zmin = std::min({zmin, vertex.z});

                    xmax = std::max({xmax, vertex.x});
                    ymax = std::max({ymax, vertex.y});
                    zmax = std::max({zmax, vertex.z});
                }
            }

            min = Vector3(xmin, ymin, zmin);
            max = Vector3(xmax, ymax, zmax);
        }

        bool intersects(const Ray& ray, double &tmin, double &tmax) const {
            const Vector3& O = ray.get_origin();
            const Vector3& u = ray.get_direction();

            double txmin = (min.x - O.x) / u.x;
            double txmax = (max.x - O.x) / u.x;
            double tymin = (min.y - O.y) / u.y;
            double tymax = (max.y - O.y) / u.y;
            double tzmin = (min.z - O.z) / u.z;
            double tzmax = (max.z - O.z) / u.z;
            tmin = std::max({std::min(txmin, txmax), std::min(tymin, tymax), std::min(tzmin, tzmax)});
            tmax = std::min({std::max(txmin, txmax), std::max(tymin, tymax), std::max(tzmin, tzmax)});

            if (tmax < 0 || tmin > tmax)
                return false;

            return true;
        }        
    };

    struct Node final {
        const size_t MIN_TRIANGLES = 5;

        Bbox bbox;
        triangle_index begin_triangle, end_triangle;

        Node *left, *right;

        Node(triangle_index begin, triangle_index end): begin_triangle(begin), end_triangle(end) { 
            left = right = nullptr;
            bbox = Bbox(begin, end);
            Vector3 diag = bbox.max - bbox.min;
            Vector3 center = bbox.min + 0.5 * diag;
            int longest_axis = diag.x > diag.y ? (diag.x > diag.z ? 0 : 2) : (diag.y > diag.z ? 1 : 2);
            triangle_index pivot_index = begin;
            for(auto it = begin; it != end; it++) {
                Vector3 barycenter = (it->get_a() + it->get_b() + it->get_c()) / 3;
                if (barycenter[longest_axis] < center[longest_axis]) {
                    std::iter_swap(it, pivot_index);
                    pivot_index++;
                }
            }

            if (pivot_index == begin || end - pivot_index <= 1 || end - begin < MIN_TRIANGLES) {
                return;
            }

            left = new Node(begin, pivot_index);
            right = new Node(pivot_index, end);
        }

        ~Node() {
            if (left != nullptr)
                delete left;
            if (right != nullptr)
                delete right;
        }
    };

    class BVH final {
    private:
        Node *root;
        const Geometry *src;

        bool intersects_top(const Ray& ray, double &tmin, double &tmax) const {
            return root->bbox.intersects(ray, tmin, tmax);
        }
    public:
        BVH(const Geometry *src, triangle_index begin, triangle_index end): src(src) {
            root = new Node(begin, end);
        }

        ~BVH() {
            delete root;
        }

        std::optional<Intersection> intersects(const Ray& ray) const {
            double tmin, tmax;
            if (!intersects_top(ray, tmin, tmax))
                return std::nullopt;

            std::optional<Intersection> closest_intersection = std::nullopt;
            std::list<Node*> nodes;
            nodes.push_front(root);
            while(!nodes.empty()) {
                Node *node = nodes.back();
                nodes.pop_back();

                if (node->left == nullptr && node->right == nullptr) {
                    for(auto it = node->begin_triangle; it != node->end_triangle; it++) {
                        std::optional<Intersection> intersection = it->intersect(ray, src);
                        if (intersection.has_value() && 
                        (
                            !closest_intersection.has_value() || 
                            (intersection.value().t < closest_intersection.value().t)
                        )) {
                            closest_intersection.emplace(std::move(intersection.value()));
                        }
                    }
                } 
                else {
                    double tmin_left, tmax_left, tmin_right, tmax_right;
                    bool intersects_left = node->left == nullptr ? false : node->left->bbox.intersects(ray, tmin_left, tmax_left);
                    bool intersects_right = node->right == nullptr ? false : node->right->bbox.intersects(ray, tmin_right, tmax_right);
                    if (intersects_left && (!closest_intersection.has_value() || tmin_left < closest_intersection.value().t)) {
                        nodes.push_back(node->left);
                    } 
                    if (intersects_right && (!closest_intersection.has_value() || tmin_right < closest_intersection.value().t)) {
                        nodes.push_back(node->right);
                    }
                }
            }

            return closest_intersection;
        }
    };
}