#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include <cfloat>
using namespace std;

class Vector
{
public:
	Vector() : x(0.0), y(0.0), z(0.0) {};
	Vector(float a_x, float a_y, float a_z) : x(a_x), y(a_y), z(a_z) {};
	Vector::Vector(float a_x) : x(a_x), y(a_x), z(a_x) {};
	Vector(const Vector& v);

	float length();

	float getAxisValue(int axis);

	Vector&	normalize();
	Vector operator=(const Vector& v);
	Vector operator+( const Vector& v );
	Vector operator-( const Vector& v );
	Vector operator*( float f );
	float  operator*(const Vector& v);   //inner product
	Vector operator/( float f );
	Vector operator%( const Vector& v); //external product
	Vector&	operator-=	(const Vector& v);
	Vector&	operator-=	(const float v);
	Vector&	operator*=	(const float v);
	Vector&	operator+=	(const float v);
	bool operator!=	(const Vector& v);
	bool operator==(const Vector& v);

	float x;
	float y;
	float z;

     friend inline
  istream&	operator >>	(istream& s, Vector& v)
	{ return s >> v.x >> v.y >> v.z; }
  
};

#endif