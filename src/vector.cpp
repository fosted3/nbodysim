#include "vector.h"
#include <math.h>
#include <iostream>

vector::vector()
{
	this -> x = 0;
	this -> y = 0;
	this -> z = 0;
}

vector::vector(double a, double b, double c)
{
	this -> x = a;
	this -> y = b;
	this -> z = c;
}

vector::~vector()
{
}

void vector::print()
{
	std::cout << this -> x << ", " << this -> y << ", " << this -> z << std::endl;
}

void vector::print_inline()
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

vector& vector::operator *= (const double &a)
{
	this -> x *= a;
	this -> y *= a;
	this -> z *= a;
	return *this;
}

vector& vector::operator /= (const double &a)
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

const double vector::get_x()
{
	return this -> x;
}

const double vector::get_y()
{
	return this -> y;
}

const double vector::get_z()
{
	return this -> z;
}

const double distance(vector *a, vector *b)
{
	return sqrt(pow(a -> get_x() - b -> get_x(), 2) + pow(a -> get_y() - b -> get_y(), 2) + pow(a -> get_z() - b -> get_z(), 2));
}

const double vector::magnitude()
{
	return sqrt(pow(this -> x, 2) + pow(this -> y, 2) + pow(this -> z, 2));
}

void vector::normalize()
{
	double mag = this -> magnitude();
	this -> x /= mag;
	this -> y /= mag;
	this -> z /= mag;
}
