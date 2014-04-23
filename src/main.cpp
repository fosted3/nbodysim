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

int main()
{
	srand(time(NULL));
	vector origin = vector(0, 0, 0);
	double size = 8192;
	quadtree* root = new quadtree(&origin, size);
	std::vector<particle*> particles;
	for (int i = 0; i < 100; i++)
	{
		generate_particle(size, particles, root);
	}
	root -> clean();
	delete root;
	return 0;
}
