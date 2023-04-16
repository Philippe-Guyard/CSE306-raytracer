#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray {
private:
    Vector3 origin, direction;
public:
    Ray(const Vector3& origin, const Vector3& direction): origin(origin), direction(direction.normalized()) {};
    
    const Vector3& get_origin() const { return origin; }
    const Vector3& get_direction() const { return direction; }
};

#endif