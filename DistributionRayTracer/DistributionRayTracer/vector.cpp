#include <cmath>
#include "vector.h"

float Vector::length()
{
	return sqrt( x * x + y * y + z * z );
}

float Vector::getAxisValue(int axis) {
	return (axis == 0) ? x : (axis == 1) ? y : z;
}

// --------------------------------------------------------------------- copy constructor
Vector::Vector(const Vector& v)
{
	x = v.x; y = v.y; z = v.z;
}

// --------------------------------------------------------------------- assignment operator
Vector Vector::operator= (const Vector& rhs) {
	if (this == &rhs)
		return (*this);
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return (*this);
}

Vector Vector::operator+(const  Vector& v )
{
	return Vector( x + v.x, y + v.y, z + v.z );
}


Vector Vector::operator-(const Vector& v )
{
	return Vector( x - v.x, y - v.y, z - v.z );
}


Vector Vector::operator*( float f )
{
	return Vector( x * f, y * f, z * f );
}

float Vector::operator*(const  Vector& v) //inner product
{
	return x * v.x + y * v.y + z * v.z;
}

Vector Vector::operator/( float f )
{
	return Vector( x / f, y / f, z / f );
}

bool Vector::operator!=(const Vector& v) {
	return (x != v.x && y != v.y && z != v.z);
}

bool Vector::operator==(const Vector& v) {
	return (x == v.x && y == v.y && z == v.z);
}

Vector&	Vector::normalize	()
{
				   float l=1.0/this->length();
				   x *= l; y *= l; z *= l;
				   return *this;
}

Vector&	Vector::operator-=(const Vector& v)
{ x-=v.x; y-=v.y; z-=v.z; return *this; }

Vector&	Vector::operator-=(const float v)
{ x-=v; y-=v; z-=v; return *this; }

Vector&	Vector::operator+=(const float v)
{ x+=v; y+=v; z+=v; return *this; }

Vector&	Vector::operator*=(const float v)
{ x*=v; y*=v; z*=v; return *this; }

Vector Vector::operator%( const Vector& v)
{
	float uX = x;
	float uY = y;
	float uZ = z;

	float vX = v.x;
	float vY = v.y;
	float vZ = v.z;

	float sX = uY * vZ - uZ * vY;
	float sY = uZ * vX - uX * vZ;
	float sZ = uX * vY - uY * vX;

	return Vector( sX, sY, sZ );
}