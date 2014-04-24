//#define THREADED
#include "vector.h"
#include "particle.h"
#include "quadtree.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <queue>
#include <math.h>
#include <stdio.h>
#include <fstream>

#ifdef THREADED

#define NUM_THREADS 4
#include <pthread.h>

struct thread_data
{
	int thread_id;
	std::vector<particle*> *particles;
	quadtree *root;
	unsigned int start;
	unsigned int end;
	double theta;
};

#endif

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
	particle* par = new particle(&temp, &null, &null, 1e10);
	particles.push_back(par);
	root -> add_particle(par);
}

void check_tree(std::vector<particle*> &particles)
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
	double percent;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		percent = (double) i * 100;
		percent /= particles.size();
		printf("%3.2f%%", percent);
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
		std::cout << "\b\b\b\b\b\b\b";
	}
}

#ifdef THREADED

void *barnes_hut_thread(void *data)
{
	struct thread_data *args;
	args = (struct thread_data*) data;
	std::queue<quadtree*> nodes;
	quadtree* node;
	quadtree* root = args -> root;
	particle* curr;
	std::vector<particle*> *particles = args -> particles;
	vector grav_to;
	double theta = args -> theta;
	for (unsigned int i = (args -> start); i < (args -> end); i++)
	{
		curr = (*particles)[i];
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
	pthread_exit(NULL);
}

#endif

void dump(unsigned int frame, std::vector<particle*> &particles)
{
	unsigned int size = sizeof(particle);
	std::string filename = std::to_string(frame);
	while (filename.length() < 4)
	{
		filename.insert(0, "0");
	}
	filename += ".dat";
	std::fstream outfile(filename, std::ios::out | std::ios::binary);
	outfile.seekp(0);
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		outfile.write((char*)particles[i], size);
	}
	outfile.close();
}

void print_particles(std::vector<particle*> &particles)
{
	std::cout << "~" << std::endl;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		particles[i] -> print_compact();
	}
}

int main()
{
	//sizeof(particle) = 88
	//sizeof(quadtree) = 144
	srand(time(NULL));
	vector origin = vector(0, 0, 0);
	double size = 2048;
	double theta = 0.5;
	double dt = 0.03333;
	unsigned int frames = 1;
	unsigned int n_p = 100000;
	quadtree* root = new quadtree(&origin, size);
	std::vector<particle*> particles;
#ifdef THREADED
	pthread_t threads[NUM_THREADS];
	struct thread_data td[NUM_THREADS];
	int rc;
#endif
	std::cout << "Generating particles..." << std::endl;
	for (unsigned int i = 0; i < n_p; i++)
	{
		generate_particle(size, particles, root);
	}
	for (unsigned int i = 0; i < frames; i++)
	{
		std::cout << "Frame " << i + 1 << std::endl;
		//root -> print_info(0);
		root -> calc_mass();
		root -> calc_com();
#ifndef THREADED
		barnes_hut(particles, root, theta);
#endif
#ifdef THREADED
		for (unsigned int i = 0; i < NUM_THREADS; i++)
		{
			td[i].thread_id = i;
			td[i].particles = &particles;
			td[i].root = root;
			td[i].start = (n_p / NUM_THREADS) * i;
			td[i].end = (n_p / NUM_THREADS) * (i + 1);
			td[i].theta = theta;
			rc = pthread_create(&threads[i], NULL, barnes_hut_thread, (void*) &td[i]);
			if (rc)
			{
				std::cout << "Could not create thread." << std::endl;
				exit(-1);
			}
		}
#endif
		update_all(particles, dt);
		check_tree(particles);
		root -> clean();
		root -> remove_redundancy();
		dump(i, particles);
	}
	delete root;
	pthread_exit(NULL);
	return 0;
}
