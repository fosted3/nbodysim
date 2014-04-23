#include "vector.h"
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

vector& vector::operator += (const vector& v)
{
	this -> x += v.x;
	this -> y += v.y;
	this -> z += v.z;
	return *this;
}

vector& vector::operator -= (const vector& v)
{
	this -> x -= v.x;
	this -> y -= v.y;
	this -> z -= v.z;
	return *this;
}

vector& vector::operator = (const vector& v)
{
	this -> x = v.x;
	this -> y = v.y;
	this -> z = v.z;
	return *this;
}
