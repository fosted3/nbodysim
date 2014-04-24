#include "vector.h"
#include "particle.h"
#include "quadtree.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>

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

int main()
{
	srand(time(NULL));
	vector origin = vector(0, 0, 0);
	double size = 1024;
	quadtree* root = new quadtree(&origin, size);
	std::vector<particle*> particles;
	/*for (int i = 0; i < 100; i++)
	{
		generate_particle(size, particles, root);
	}*/
	vector a = vector(-256, -256, -256);
	vector b = vector(-256, -256,  256);
	vector c = vector(-256,  256, -256);
	vector d = vector(-256,  256,  256);
	vector e = vector( 256, -256, -256);
	vector f = vector( 256, -256,  256);
	vector g = vector( 256,  256, -256);
	vector h = vector( 256,  256,  256);
	vector q = vector( 384,  384,  384);
	particle *i = new particle(&a, &q, &origin, 1.0);
	particle *j = new particle(&b, &origin, &origin, 1.0);
	particle *k = new particle(&c, &origin, &origin, 1.0);
	particle *l = new particle(&d, &origin, &origin, 1.0);
	particle *m = new particle(&e, &origin, &origin, 1.0);
	particle *n = new particle(&f, &origin, &origin, 1.0);
	particle *o = new particle(&g, &origin, &origin, 1.0);
	particle *p = new particle(&h, &origin, &origin, 10.0);
	root -> add_particle(i);
	root -> add_particle(j);
	root -> add_particle(k);
	root -> add_particle(l);
	root -> add_particle(m);
	root -> add_particle(n);
	root -> add_particle(o);
	root -> add_particle(p);
	particles.push_back(i);
	root -> calc_mass();
	root -> calc_com();
	root -> clean();
	root -> print_info(0);	
	check_tree(particles, root);
	update_all(particles, 1);
	check_tree(particles, root);
	root -> print_info(0);
	root -> clean();
	delete root;
	return 0;
}
