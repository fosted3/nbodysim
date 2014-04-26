#include "quadtree.h"
#include <cstddef>
#include <cassert>
#include <iostream>
#include <stdlib.h>

quadtree::quadtree()
{
}

quadtree::quadtree(vector *cen, double sid, quadtree *par)
{
	this -> center = *cen;
	this -> side = sid;
	this -> parent = par;
	for (int i = 0; i < 8; i++)
	{
		this -> children[i] = NULL;
	}
	this -> p = NULL;
}

quadtree::quadtree(vector *cen, double sid)
{
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
	double half = (this -> side) / 2;
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
	/*if (this -> p != NULL)
	{
		delete p;
	}*/
}

void quadtree::print_info(void)
{
	std::cout << "Quadtree node @ " << this << ", center @ ";
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

bool quadtree::add_particle(particle *par)
{
	if (!this -> inside(par))
	{
		if (this -> parent == NULL)
		{
			return false;
		}
		return (this -> parent -> add_particle(par));
	}
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
			par -> set_container(this);
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
			return (this -> children[i] -> add_particle(par));
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
	return true;
}

double quadtree::get_mass(void)
{
	return this -> mass;
}

void quadtree::calc_mass(void) //call on root
{
	if (this -> p == NULL)
	{
		this -> mass = 0;
		for (int i = 0; i < 8; i++)
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

void quadtree::calc_com(void) //call on root
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

vector* quadtree::get_com(void)
{
	return &(this -> com);
}

bool quadtree::clean(void) //only call on root
{
	if (this -> p != NULL)
	{
		return false;
	}
	bool empty = true;
	for (int i = 0; i < 8; i++)
	{
		if (this -> children[i] != NULL)
		{
			if (this -> children[i] -> clean())
			{
				std::cout << "Child at " << children[i] << " allocated but empty. Cleaning." << std::endl;
				delete this -> children[i];
				this -> children[i] = NULL;
			}
			else
			{
				empty = false;
			}
		}
	}
	return empty;
}

bool quadtree::inside(particle* par)
{
	vector *temp = par -> get_pos();
	if (abs(temp -> get_x() - this -> center.get_x()) > (this -> side) / 2)
	{
		return false;
	}
	if (abs(temp -> get_y() - this -> center.get_y()) > (this -> side) / 2)
	{
		return false;
	}
	if (abs(temp -> get_z() - this -> center.get_z()) > (this -> side) / 2)
	{
		return false;
	}
	return true;
}

void quadtree::release_particle(void)
{
	this -> p = NULL;
}

quadtree* quadtree::get_parent(void)
{
	return this -> parent;
}

double quadtree::get_side(void)
{
	return this -> side;
}

quadtree* quadtree::get_child(int i)
{
	return (this -> children[i]);
}

particle* quadtree::get_particle(void)
{
	return this -> p;
}

void quadtree::remove_redundancy(void)
{
	if (this -> p != NULL)
	{
		return;
	}
	bool single_node_particle = false;
	int first_index = -1;
	int last_index = -1;
	for (int i = 0; i < 8; i++)
	{
		if (this -> children[i] != NULL)
		{
			if (single_node_particle == false)
			{
				first_index = i;
				single_node_particle = true;
			}
			else
			{
				single_node_particle = false;
				break;
			}
			if (this -> children[i] -> get_particle() != NULL)
			{
				last_index = i;
			}
		}
	}
	if (single_node_particle && last_index != -1)
	{
		particle *temp = this -> children[last_index] -> get_particle();
		this -> children[last_index] -> release_particle();
		delete (this -> children[last_index]);
		this -> add_particle(temp);
	}
	else if(first_index != -1)
	{
		for (int i = first_index; i < 8; i++)
		{
			if (this -> children[i] != NULL)
			{
				this -> children[i] -> remove_redundancy();
			}
		}
	}
}
