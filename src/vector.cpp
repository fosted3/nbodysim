#include "vector.h"
#include <math.h>
#include <iostream>

#ifdef DOUBLE
#ifndef datatype
#define datatype double
#endif
#endif
#ifdef FLOAT
#ifndef datatype
#define datatype float
#endif
#endif

vector::vector()
{
	this -> x = 0;
	this -> y = 0;
	this -> z = 0;
}

vector::vector(datatype a, datatype b, datatype c)
{
	this -> x = a;
	this -> y = b;
	this -> z = c;
}

vector::~vector()
{
}

void vector::print(void)
{
	std::cout << this -> x << ", " << this -> y << ", " << this -> z << std::endl;
}

void vector::print_inline(void)
{
	std::cout << this -> x << ", " << this -> y << ", " << this -> z;
}

vector& vector::operator += (const vector &v)
{
	this -> x += v.x;
	this -> y += v.y;
	this -> z += v.z;
	return *this;
}

vector& vector::operator -= (const vector &v)
{
	this -> x -= v.x;
	this -> y -= v.y;
	this -> z -= v.z;
	return *this;
}

vector& vector::operator *= (const datatype &a)
{
	this -> x *= a;
	this -> y *= a;
	this -> z *= a;
	return *this;
}

vector& vector::operator /= (const datatype &a)
{
	this -> x /= a;
	this -> y /= a;
	this -> z /= a;
	return *this;
}

vector& vector::operator = (const vector &v)
{
	this -> x = v.x;
	this -> y = v.y;
	this -> z = v.z;
	return *this;
}

datatype vector::get_x(void)
{
	return this -> x;
}

datatype vector::get_y(void)
{
	return this -> y;
}

datatype vector::get_z(void)
{
	return this -> z;
}

datatype distance(vector *a, vector *b)
{
	return sqrt(pow(a -> get_x() - b -> get_x(), 2) + pow(a -> get_y() - b -> get_y(), 2) + pow(a -> get_z() - b -> get_z(), 2));
}

datatype vector::magnitude(void)
{
	return sqrt(pow(this -> x, 2) + pow(this -> y, 2) + pow(this -> z, 2));
}

void vector::normalize(void)
{
	datatype mag = this -> magnitude();
	this -> x /= mag;
	this -> y /= mag;
	this -> z /= mag;
}

vector cross(vector &a, vector &b)
{
	return vector(a.get_y()*b.get_z() - a.get_z()*b.get_y(), a.get_z()*b.get_x() - a.get_x()*b.get_z(), a.get_x()*b.get_y() - a.get_y()*b.get_x());
}

void vector::scale(datatype scale_x, datatype scale_y, datatype scale_z)
{
	this -> x *= scale_x;
	this -> y *= scale_y;
	this -> z *= scale_z;
}
