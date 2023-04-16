#ifndef RNG_HPP
#define RNG_HPP

#include <random>
#include "vector.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);
class Rng {
public:
    static double random_real() {
        return uniform(engine);
    }

    static Vector3 random_cos(const Vector3& N) {
        double r1 = random_real();
        double r2 = random_real();
        double x = std::cos(2 * M_PI * r1) * std::sqrt(1 - r2);
        double y = std::sin(2 * M_PI * r1) * std::sqrt(1 - r2);
        double z = std::sqrt(r2);
        // Find the smallest component of N and set it to zero
        double new_x = N.x, new_y = N.y, new_z = N.z;
        if (std::abs(N.x) < std::abs(N.y) && std::abs(N.x) < std::abs(N.z)) {
            new_x = 0;
            new_y = N.z;
            new_z = -N.y;
        }
        else if (std::abs(N.y) < std::abs(N.x) && std::abs(N.y) < std::abs(N.z)) {
            new_y = 0;
            new_x = N.z;
            new_z = -N.x;
        }
        else {
            new_z = 0;
            new_x = N.y;
            new_y = -N.x;
        }
        Vector3 T1(new_x, new_y, new_z);
        T1 = T1.normalized();
        Vector3 T2 = N.cross(T1);
        return x * T1 + y * T2 + z * N;
    }
};

#endif // RNG_HPP