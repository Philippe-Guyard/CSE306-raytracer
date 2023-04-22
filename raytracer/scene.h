#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "rng.hpp"
#include "geometry.h"

struct LightSource {
    Vector3 position;
    double intensity;

    LightSource(const Vector3& position, double intensity): position(position), intensity(intensity) {};
};

class Scene {
private:
    // TODO: Proper storage of this pointers is needed
    // NOTE: Storing pointers is required by polymorphism
    std::vector<Geometry*> objects;
    LightSource light_source;
    Vector3 bg_color;

    std::optional<Intersection> find_closest_intersection(const Ray& r) const {
        std::optional<Intersection> closest_intersection = std::nullopt;
        for(auto object: objects) {
            std::optional<Intersection> intersection = object->intersect(r);
            if (intersection.has_value() && 
                (
                    !closest_intersection.has_value() || 
                    (intersection.value().t < closest_intersection.value().t)
                )) {
                closest_intersection.emplace(std::move(intersection.value()));
            }
        }

        return closest_intersection;
    }

    Vector3 reflected_get_color(const Ray& r, const Vector3& P_inter, const Vector3& N_inter, int max_bounces) const {
        Vector3 reflected_direction = r.get_direction() - 2 * N_inter.dot(r.get_direction()) * N_inter;
        Ray reflected_ray = Ray(P_inter + 0.001 * reflected_direction, reflected_direction);
        return get_color(reflected_ray, max_bounces - 1);
    }

    Vector3 refracted_get_color(const Ray& r, const Vector3& P_inter, const Vector3& N_inter, int max_bounces) const {
        // TODO: Refraction indices should be stored in Geometry                 
        double n1 = 1;
        double n2 = 1.4;
        double n = n1 / n2;
        Vector3 correctedN = N_inter;
        if (r.get_direction().dot(N_inter) > 0) {
            n = n2 / n1;
            correctedN = -N_inter;
        }

        double cosI = correctedN.dot(r.get_direction());
        Vector3 tangent = n * (r.get_direction() - cosI * correctedN);
        double rad = 1 - n * n * (1 - cosI * cosI);
        if (rad < 0) {
            Vector3 reflected_direction = r.get_direction() - 2 * correctedN.dot(r.get_direction()) * correctedN;
            Ray reflected_ray = Ray(P_inter - 0.001 * reflected_direction, reflected_direction);
            return get_color(reflected_ray, max_bounces - 1);
        }
        Vector3 normal = -correctedN * sqrt(rad);
        Ray refracted_ray = Ray(P_inter - 0.001 * correctedN, tangent + normal);
        return get_color(refracted_ray, max_bounces - 1);
    }
public:
    Scene(const LightSource& light_source, 
          const Vector3& background_color = Vector3(0, 0, 0)): light_source(light_source),
                                                               bg_color(background_color) {};

    ~Scene() {
        objects.clear();
    };

    void add_object(Geometry* object) {
        objects.push_back(object);
    }

    Vector3 get_color(const Ray& r, int max_bounces = 5) const {
        if (max_bounces < 0)
            return bg_color;

        auto intersection = find_closest_intersection(r);
        if (!intersection.has_value())
            return bg_color;
        
        if (intersection.value().object == nullptr) {
            std::cout << "Error: object is nullptr" << std::endl;
            return bg_color;
        }
        if (intersection.value().t < 0) {
            std::cout << "Error: t < 0" << std::endl;
            return bg_color;
        }

        const Vector3& N = intersection.value().normal;
        const Vector3& P = intersection.value().point;
        const Vector3& rho = intersection.value().object->get_color();

        if (intersection.value().object->is_mirror()) 
            return reflected_get_color(r, P, N, max_bounces);
        else if (intersection.value().object->is_transparent()) 
            return refracted_get_color(r, P, N, max_bounces);

        Vector3 random_cos = Rng::random_cos(N);
        Ray random_ray = Ray(P + 0.001 * N, random_cos);
        Vector3 indirect_light = get_color(random_ray, max_bounces - 1) * rho;

        Vector3 light_direction = (light_source.position - P).normalized();
        Ray shadowRay = Ray(P + 0.001 * light_direction, light_direction);
        auto shadow_intersection = find_closest_intersection(shadowRay);
        double d2 = P.dist2(light_source.position);
        if (shadow_intersection.has_value()) {
            double t = shadow_intersection.value().t;
            if (t * t < d2)
                return indirect_light;
        }

        Vector3 direct_light = light_source.intensity / (4 * M_PI * d2) * (rho / M_PI) * std::max(0., N.dot(light_direction));
        return direct_light + indirect_light;
    }
};

#endif 