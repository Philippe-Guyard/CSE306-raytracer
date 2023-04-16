#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

class Vector3 {
public:
    double x, y, z;
    Vector3(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z) {};
    // Not necessary since C++ has a default copy constructor
    // Vector3(const Vector3& v): x(v.x), y(v.y), z(v.z) {};

	double norm2() const;
	double norm() const;

    double dot(const Vector3& b) const;
    double dist2(const Vector3& b) const;
    double dist(const Vector3& b) const;
    Vector3 cross(const Vector3& b) const;

    Vector3 normalized() const;
    Vector3 operator-() const;

    void normalize();
    void operator+=(const Vector3& b);
};

Vector3 operator+(const Vector3& a, const Vector3& b);
Vector3 operator-(const Vector3& a, const Vector3& b);
Vector3 operator*(const double a, const Vector3& b);
Vector3 operator*(const Vector3& a, const double b);
Vector3 operator/(const Vector3& a, const double b);
Vector3 operator*(const Vector3& a, const Vector3& b);

Vector3 gamma_correct(const Vector3& v);

#endif