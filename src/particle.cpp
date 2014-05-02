#include "particle.h"
#include <iostream>

particle::particle()
{
}

particle::particle(vector *pos, vector *vel, vector *acc, double mas)
{
	this -> position = *pos;
	this -> velocity = *vel;
	this -> acceleration = *acc;
	this -> mass = mas;
}

double particle::get_mass(void)
{
	return this -> mass;
}

vector* particle::get_pos(void)
{
	return &(this -> position);
}

void particle::set_acc_zero(void)
{
	vector temp = vector(0, 0, 0);
	this -> acceleration = temp;
}

void particle::set_acc_offset(vector *off)
{
	this -> acceleration += *off;
}

void particle::update(double dt)
{
	vector temp = this -> acceleration;
	temp *= dt;
	this -> velocity += temp;
	temp = this -> velocity;
	temp *= dt;
	this -> position += temp;
}

void particle::print(void)
{
	std::cout << "Particle @ " << this << ", pos: ";
	this -> position.print_inline();
	std::cout << ", vel: ";
	this -> velocity.print_inline();
	std::cout << ", acc: ";
	this -> acceleration.print();
}
