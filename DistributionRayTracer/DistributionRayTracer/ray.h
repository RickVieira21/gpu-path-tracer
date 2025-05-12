#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray
{
public:
	Ray() {};
	Ray(const Vector& o, const Vector& dir ) : origin(o), direction(dir) {};

	Vector origin;
	Vector direction;
};
#endif