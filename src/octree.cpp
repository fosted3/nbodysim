#include "octree.h"
#include "thread_functions.h"
#include <cstddef>
#include <cassert>
#include <iostream>
#include <pthread.h>
//#include <stdlib.h>
//#include <cmath>


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

octree::octree()
{
}

octree::octree(vector *cen, datatype sid, octree *par)
{
	this -> center = *cen;
	this -> side = sid;
	this -> parent = par;
	for (int i = 0; i < 8; i++)
	{
		this -> children[i] = NULL;
	}
	this -> p = NULL;
	this -> com = *cen;
	this -> mass = 0;
}

octree::octree(vector *cen, datatype sid)
{
	this -> center = *cen;
	this -> side = sid;
	this -> parent = NULL;
	for (int i = 0; i < 8; i++)
	{
		this -> children[i] = NULL;
	}
	this -> p = NULL;
	this -> com = *cen;
	this -> mass = 0;
}

void octree::allocate_child(unsigned int i)
{	/*
	 z  y  x  i
	-1 -1 -1  0
	-1 -1  1  1
	-1  1 -1  2
	-1  1  1  3
	 1 -1 -1  4
	 1 -1  1  5
	 1  1 -1  6
	 1  1  1  7*/
	int x = -1;
	int y = -1;
	int z = -1;
	datatype half = (this -> side) / 2;
	if (i & 1)
	{
		x = 1;
	}
	if (i & 2)
	{
		y = 1;
	}
	if (i & 4)
	{
		z = 1;
	}
	assert(this -> children[i] == NULL);
	vector temp = this -> center;
	temp += vector(x * half / 2, y * half / 2, z * half / 2);
	children[i] = new octree(&temp, half, this);
}

octree::~octree()
{
	for (int i = 0; i < 8; i++)
	{
		if (children[i] != NULL)
		{
			delete children[i];
		}
	}
	/*if (this -> p != NULL)
	{
		delete p;
	}*/
}

void octree::print_info(void)
{
	std::cout << "octree node @ " << this << ", center @ ";
	this -> center.print_inline();
	std::cout << ", size of " << this -> side;
	if (this -> p != NULL)
	{
		std::cout << ", ";
		this -> p -> print();
	}
	else
	{
		std::cout << std::endl;
	}
}

void octree::print_info(unsigned int depth)
{
	for (unsigned int i = 0; i < depth; i++)
	{
		std::cout << "\t";
	}
	this -> print_info();
	for (unsigned int i = 0; i < 8; i ++)
	{
		if (this -> children[i] != NULL)
		{
			children[i] -> print_info(depth + 1);
		}
	}
}

void octree::add_particle(particle *par)
{
	if (this -> p == NULL)
	{
		bool leaf = true;
		for (int i = 0; i < 8; i++)
		{
			if (this -> children[i] != NULL)
			{
				leaf = false;
				break;
			}
		}
		if (leaf)
		{
			this -> p = par;
		}
		else
		{
			int i = 0;
			if (par -> get_pos() -> get_x() > this -> center.get_x())
			{
				i += 1;
			}
			if (par -> get_pos() -> get_y() > this -> center.get_y())
			{
				i += 2;
			}
			if (par -> get_pos() -> get_z() > this -> center.get_z())
			{
				i += 4;
			}
			if (this -> children[i] == NULL)
			{
				this -> allocate_child(i);
			}
			this -> children[i] -> add_particle(par);
		}
	}
	else
	{
		int i = 0;
		if (this -> p -> get_pos() -> get_x() > this -> center.get_x())
		{
			i += 1;
		}
		if (this -> p -> get_pos() -> get_y() > this -> center.get_y())
		{
			i += 2;
		}
		if (this -> p -> get_pos() -> get_z() > this -> center.get_z())
		{
			i += 4;
		}
		if (this -> children[i] == NULL)
		{
			this -> allocate_child(i);
		}
		this -> children[i] -> add_particle(this -> p);
		this -> p = NULL;
		this -> add_particle(par);
	}
}

datatype octree::get_mass(void)
{
	return this -> mass;
}

void octree::calc_mass(void) //call on root *only*
{
	if (this -> p == NULL)
	{
		this -> mass = 0;
		for (unsigned int i = 0; i < 8; i++)
		{
			if (this -> children[i] != NULL)
			{
				this -> children[i] -> calc_mass();
				this -> mass += this -> children[i] -> get_mass();
			}
		}
	}
	else
	{
		this -> mass = this -> p -> get_mass();
	}
}

void *calc_mass_thread(void *obj)
{
	octree* target = (octree*) obj;
	target -> calc_mass();
	pthread_exit(NULL);
}

void octree::calc_mass_threaded(void) //call on root *only*
{
	pthread_t threads[8];
	octree* objs[8];
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		objs[i] = this -> children[i];
		create_thread(&threads[i], NULL, calc_mass_thread, (void*) objs[i]);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		pthread_join(threads[i], NULL);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		this -> mass += this -> children[i] -> get_mass();
	}
}

void octree::calc_com(void) //call on root *only*
{
	if (this -> p != NULL)
	{
		this -> com = *(this -> p -> get_pos());
	}
	else
	{
		vector temp = vector(0, 0, 0);
		for (int i = 0; i < 8; i++)
		{
			if (this -> children[i] != NULL)
			{
				this -> children[i] -> calc_com();
				vector temp2 = *(this -> children[i] -> get_com());
				temp2 *= this -> children[i] -> get_mass();
				temp += temp2;
			}
		}
		temp /= this -> mass;
		this -> com = temp;
	}
	//this -> com.print();
}

void *calc_com_thread(void *obj)
{
	octree *target = (octree*) obj;
	target -> calc_com();
	pthread_exit(NULL);
}

void octree::calc_com_threaded(void) //call on root *only*
{
	pthread_t threads[8];
	octree* objs[8];
	vector temp = vector(0, 0, 0);
	vector temp2;
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		objs[i] = this -> children[i];
		create_thread(&threads[i], NULL, calc_com_thread, (void*) objs[i]);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		pthread_join(threads[i], NULL);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		if (children[i] == NULL) { continue; }
		temp2 = *(this -> children[i] -> get_com());
		temp2 *= this -> children[i] -> get_mass();
		temp += temp2;
	}
	temp /= this -> mass;
	this -> com = temp;
}

vector* octree::get_com(void)
{
	return &(this -> com);
}

bool octree::inside(particle* par)
{
	vector *temp = par -> get_pos();
	if ((temp -> get_x() - this -> center.get_x()) > ((this -> side) / 2.0) || (temp -> get_x() - this -> center.get_x()) < ((-1.0 * this -> side) / 2.0))
	{
		return false;
	}
	if ((temp -> get_y() - this -> center.get_y()) > ((this -> side) / 2.0) || (temp -> get_y() - this -> center.get_y()) < ((-1.0 * this -> side) / 2.0))
	{
		return false;
	}
	if ((temp -> get_z() - this -> center.get_z()) > ((this -> side) / 2.0) || (temp -> get_z() - this -> center.get_z()) < ((-1.0 * this -> side) / 2.0))
	{
		return false;
	}
	return true;
}

void octree::release_particle(void)
{
	this -> p = NULL;
}

octree* octree::get_parent(void)
{
	return this -> parent;
}

datatype octree::get_side(void)
{
	return this -> side;
}

octree* octree::get_child(unsigned int i)
{
	return (this -> children[i]);
}

particle* octree::get_particle(void)
{
	return this -> p;
}

void octree::release_child(unsigned int i) //you better make sure you've deleted it first...
{
	this -> children[i] = NULL;
}
