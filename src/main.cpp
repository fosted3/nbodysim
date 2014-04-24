#include "vector.h"
#include "particle.h"
#include "quadtree.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <queue>
#include <math.h>

double random_double(double low, double high)
{
	double r = rand();
	r /= RAND_MAX;
	r *= high - low;
	r += low;
	return r;
}

vector random_vector(double low, double high)
{
	vector rv = vector(random_double(low, high), random_double(low, high), random_double(low, high));
	return rv;
}

void generate_particle(double size, std::vector<particle*> &particles, quadtree *root)
{
	vector temp = random_vector(-1*(size/2), size/2);
	vector null = vector(0, 0, 0);
	particle* par = new particle(&temp, &null, &null, 1.0);
	particles.push_back(par);
	root -> add_particle(par);
}

void check_tree(std::vector<particle*> &particles, quadtree *root)
{
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		if (!particles[i] -> get_container() -> inside(particles[i]))
		{	
			particles[i] -> get_container() -> release_particle();
			if (!particles[i] -> get_container() -> get_parent() -> add_particle(particles[i]))
			{
				delete particles[i];
				particles.erase(particles.begin() + i);
				i--;
			}
		}
	}
}	

void update_all(std::vector<particle*> &particles, double dt)
{
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		particles[i] -> update(dt);
	}
}

vector gravity(particle* par, quadtree* node)
{
	if (par == node -> get_particle())
	{
		return vector(0, 0, 0);
	}
	vector acc = *(node -> get_com());
	acc -= *(par -> get_pos());
	double r_sq = pow(acc.magnitude(), -2);
	acc.normalize();
	acc *= node -> get_mass();
	acc *= 6.67384e-11;
	acc *= r_sq;
	return acc;
}

void barnes_hut(std::vector<particle*> &particles, quadtree *root, double theta)
{
	std::queue<quadtree*> nodes;
	quadtree* node;
	particle* curr;
	vector grav_to;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		curr = particles[i];
		curr -> set_acc_zero();
		nodes.push(root);
		while (!nodes.empty())
		{
			node = nodes.front();
			nodes.pop();
			if (node -> get_side() / distance(node -> get_com(), curr -> get_pos()) > theta && node -> get_particle() == NULL)
			{
				for (unsigned int i = 0; i < 8; i++)
				{
					if (node -> get_child(i) != NULL)
					{
						nodes.push(node -> get_child(i));
					}
				}
			}
			else
			{
				grav_to = gravity(curr, node);
				curr -> set_acc_offset(&grav_to);
			}
		}
	}
}

int main()
{
	srand(time(NULL));
	vector origin = vector(0, 0, 0);
	double size = 1024;
	double theta = 0.5;
	quadtree* root = new quadtree(&origin, size);
	std::vector<particle*> particles;
	vector a = vector(-256, -256, -256);
	vector b = vector(-256, -256,  256);
	vector c = vector(-256,  256, -256);
	vector d = vector(-256,  256,  256);
	vector e = vector( 256, -256, -256);
	vector f = vector( 256, -256,  256);
	vector g = vector( 256,  256, -256);
	vector h = vector( 256,  256,  256);
	vector q = vector( 384,  384,  384);
	particle *i = new particle(&a, &origin, &origin, 1.0);
	particle *j = new particle(&b, &origin, &origin, 1.0);
	particle *k = new particle(&c, &origin, &origin, 1.0);
	particle *l = new particle(&d, &origin, &origin, 1.0);
	particle *m = new particle(&e, &origin, &origin, 1.0);
	particle *n = new particle(&f, &origin, &origin, 1.0);
	particle *o = new particle(&g, &origin, &origin, 1.0);
	particle *p = new particle(&h, &origin, &origin, 1e16);
	root -> add_particle(i);
	root -> add_particle(j);
	root -> add_particle(k);
	root -> add_particle(l);
	root -> add_particle(m);
	root -> add_particle(n);
	root -> add_particle(o);
	root -> add_particle(p);
	particles.push_back(j);
	particles.push_back(k);
	particles.push_back(l);
	particles.push_back(m);
	particles.push_back(n);
	particles.push_back(o);
	particles.push_back(p);
	for (int i = 0; i < 50; i++)
	{
		root -> print_info(0);
		root -> calc_mass();
		root -> calc_com();
		barnes_hut(particles, root, theta);
		update_all(particles, 1);
		check_tree(particles, root);
		root -> clean();
	}
	delete root;
	return 0;
}
