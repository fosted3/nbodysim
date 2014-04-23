#include "particle.h"
#include <iostream>

particle::particle()
{
}

particle::particle(vector *pos, vector *vel, vector *acc, double mas, quadtree *con)
{
	this -> position = *pos;
	this -> velocity = *vel;
	this -> acceleration = *acc;
	this -> mass = mas;
	this -> container = con;
}

void particle::set_container(quadtree *con)
{
	this -> container = con;
}

double particle::get_mass()
{
	return this -> mass;
}

vector* particle::get_pos()
{
	return &(this -> position);
}

void particle::set_acc_zero()
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

void particle::print()
{
	std::cout << "Particle @ " << this << ", pos: ";
	this -> position.print_inline();
	std::cout << ", vel: ";
	this -> velocity.print_inline();
	std::cout << ", acc: ";
	this -> acceleration.print();
}
