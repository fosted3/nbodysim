#include "vector.h"
#include <iostream>

vector::vector()
{
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
