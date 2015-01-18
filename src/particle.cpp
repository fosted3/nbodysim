#include "particle.h"
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

particle::particle()
{
}

particle::particle(vector *pos, vector *vel, vector *acc, datatype mas)
{
	this -> position = *pos;
	this -> velocity = *vel;
	this -> acceleration = *acc;
	this -> mass = mas;
}

particle::particle(particle &p)
{
	this -> position = p.position;
	this -> velocity = p.velocity;
	this -> acceleration = p.acceleration;
	this -> mass = p.mass;
}

datatype particle::get_mass(void)
{
	return this -> mass;
}

vector* particle::get_pos(void)
{
	return &(this -> position);
}

vector* particle::get_vel(void)
{
	return &(this -> velocity);
}

vector* particle::get_acc(void)
{
	return &(this -> acceleration);
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

void particle::update(datatype dt)
{
	vector temp = this -> acceleration;
	temp *= dt;
	this -> velocity += temp;
	temp = this -> velocity;
	temp *= dt;
	this -> position += temp;
}

void particle::update(datatype &dv, datatype &dp)
{
	dv = this -> acceleration.magnitude();
	dp = this -> velocity.magnitude();
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
