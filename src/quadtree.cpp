#include "quadtree.h"
#include <cassert>
#include <cstddef>
#include <iostream>

quadtree::quadtree()
{
}

quadtree::quadtree(vector *cen, unsigned long sid, quadtree *par)
{
	assert(((sid % 2 == 0) || (sid == 1)) && sid > 0);
	this -> center = *cen;
	this -> side = sid;
	this -> parent = par;
	for (int i = 0; i < 8; i++)
	{
		this -> children[i] = NULL;
	}
	this -> p = NULL;
}

quadtree::quadtree(vector *cen, unsigned long sid)
{
	assert(((sid % 2 == 0) || (sid == 1)) && sid > 0);
	this -> center = *cen;
	this -> side = sid;
	this -> parent = NULL;
	for (int i = 0; i < 8; i++)
	{
		this -> children[i] = NULL;
	}
	this -> p = NULL;
}

void quadtree::allocate_child(int i)
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
	long half = (this -> side) / 2;
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
	children[i] = new quadtree(&temp, half, this);
}

quadtree::~quadtree()
{
	for (int i = 0; i < 8; i++)
	{
		if (children[i] != NULL)
		{
			delete children[i];
		}
	}
}

void quadtree::print_info()
{
	std::cout << "Quadtree node @ " << this << ", center @ ";
	this -> center.print_inline();
	std::cout << ", size of " << this -> side << std::endl;
}	

void quadtree::print_info(int depth)
{
	for (int i = 0; i < depth; i++)
	{
		std::cout << "\t";
	}
	this -> print_info();
	for (int i = 0; i < 8; i ++)
	{
		if (this -> children[i] != NULL)
		{
			children[i] -> print_info(depth + 1);
		}
	}
}

void quadtree::add_particle(particle *par)
{
	if (this -> p == NULL)
	{
	}
}

double quadtree::get_mass()
{
	return this -> mass;
}

void quadtree::update_mass()
{
	if (this -> p == NULL)
	{
		this -> mass = 0;
		for (int i = 0; i < 8; i++)
		{
			if (this -> children[i] != NULL)
			{
				this -> children[i] -> update_mass();
				this -> mass += this -> children[i] -> get_mass();
			}
		}
	}
	else
	{
		this -> mass = this -> p -> get_mass();
	}
}
