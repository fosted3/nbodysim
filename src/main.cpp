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
#include <cassert>

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

struct settings
{
	bool read_existing;
	bool use_seed;
	bool display_progress;
	bool dump_binary;
	bool dump_plaintext;
	bool overwrite_data;
	unsigned int num_particles;
	unsigned int num_frames;
	double size;
	double theta;
	double dt;
	unsigned long seed;
	double min_mass;
	double max_mass;
	double min_vel;
	double max_vel;
};

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

void generate_particle(settings &s, std::vector<particle*> &particles, quadtree *root)
{
	vector temp = random_vector(-1*(s.size/2), s.size/2);
	vector null = vector(0, 0, 0);
	double mass = random_double(s.min_mass, s.max_mass);
	vector vec = random_vector(s.min_vel, s.max_vel);
	particle* par = new particle(&temp, &vec, &null, mass);
	particles.push_back(par);
	root -> add_particle(par);
}

void check_tree(std::vector<particle*> &particles, quadtree *root)
{
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		if (!root -> inside(particles[i]))
		{
			delete particles[i];
			particles.erase(particles.begin() + i);
			i--;
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

void barnes_hut(std::vector<particle*> &particles, quadtree *root, double theta, bool print)
{
	std::queue<quadtree*> nodes;
	quadtree* node;
	particle* curr;
	vector grav_to;
	double percent;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		if (print)
		{
			percent = (double) i * 100;
			percent /= particles.size();
			printf("%3.3f%%", percent);
		}
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
		if (print)
		{
			std::cout << "\b\b\b\b\b\b\b\b";
		}
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

bool file_exists(const char *filename)
{
	std::ifstream ifile(filename);
	return ifile;
}

void dump(unsigned int frame, std::vector<particle*> &particles, bool overwrite)
{
	unsigned int size = sizeof(particle);
	std::string filename = std::to_string(frame);
	while (filename.length() < 4)
	{
		filename.insert(0, "0");
	}
	filename += ".dat";
	filename.insert(0, "./data/");
	if (!overwrite && file_exists(filename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << filename << ". Aborting." << std::endl;
		exit(1);
	}
	if (file_exists(filename.c_str()))
	{
		assert(remove(filename.c_str()) == 0);
	}
	std::fstream outfile(filename, std::ios::out | std::ios::binary);
	outfile.seekp(0);
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		outfile.write((char*)particles[i], size);
	}
	outfile.close();
}

void dump_plaintext(unsigned int frame, std::vector<particle*> &particles, bool overwrite)
{
	std::string filename = std::to_string(frame);
	vector temp;
	while (filename.length() < 4)
	{
		filename.insert(0, "0");
	}
	filename += ".txt";
	filename.insert(0, "./data/");
	if (!overwrite && file_exists(filename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << filename << ". Aborting." << std::endl;
		exit(1);
	}
	if (file_exists(filename.c_str()))
	{
		assert(remove(filename.c_str()) == 0);
	}
	std::fstream outfile(filename, std::ios::out);
	outfile.seekp(0);
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		temp = *(particles[i] -> get_pos());
		outfile << temp.get_x() << ", " << temp.get_y() << ", " << temp.get_z() << "\n";
	}
	outfile.close();
}

void read_settings(settings &s, const char* sfile)
{
	std::cout << "Reading config file " << sfile << std::endl;
	std::string var;
	if (file_exists(sfile))
	{
		std::ifstream cfg(sfile);
		while (cfg)
		{
			cfg >> var;
			if (var.compare("read_existing") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.read_existing = true;
				}
				else
				{
					s.read_existing = false;
				}
			}
			else if (var.compare("use_seed") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.use_seed = true;
				}
				else
				{
					s.use_seed = false;
				}
			}
			else if (var.compare("display_progress") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.display_progress = true;
				}
				else
				{
					s.display_progress = false;
				}
			}
			else if (var.compare("dump_binary") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.dump_binary = true;
				}
				else
				{
					s.dump_binary = false;
				}
			}
			else if (var.compare("dump_plaintext") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.dump_plaintext = true;
				}
				else
				{
					s.dump_plaintext = false;
				}
			}
			else if (var.compare("overwrite_data") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0)
				{
					s.overwrite_data = true;
				}
				else
				{
					s.overwrite_data = false;
				}
			}
			else if (var.compare("num_particles") == 0)
			{
				cfg >> s.num_particles;
			}
			else if (var.compare("num_frames") == 0)
			{
				cfg >> s.num_frames;
			}
			else if (var.compare("size") == 0)
			{
				cfg >> s.size;
			}
			else if (var.compare("theta") == 0)
			{
				cfg >> s.theta;
			}
			else if (var.compare("dt") == 0)
			{
				cfg >> s.dt;
			}
			else if (var.compare("seed") == 0)
			{
				cfg >> s.seed;
			}
			else if (var.compare("min_mass") == 0)
			{
				cfg >> s.min_mass;
			}
			else if (var.compare("max_mass") == 0)
			{
				cfg >> s.max_mass;
			}
			else if (var.compare("min_vel") == 0)
			{
				cfg >> s.min_vel;
			}
			else if (var.compare("max_vel") == 0)
			{
				cfg >> s.max_vel;
			}
			else
			{
				std::cerr << "Unrecognized variable " << var << std::endl;
			}
		}
	}
	else
	{
		std::cerr << sfile << " doesn't exist or can't be read." << std::endl;
		exit(1);
	}
}

void set_default(settings &s)
{
	s.read_existing = false;
	s.use_seed = false;
	s.display_progress = true;
	s.dump_binary = false;
	s.dump_plaintext = true;
	s.overwrite_data = false;
	s.num_particles = 5000;
	s.num_frames = 300;
	s.size = 512;
	s.theta = 0.5;
	s.dt = 0.0333;
	s.seed = 0xDEADBEEF;
	s.min_mass = 5e10;
	s.max_mass = 5e11;
	s.min_vel = 0;
	s.max_vel = 0;
}

int main(int argc, char **argv)
{
	settings config;
	set_default(config);
	if (argc == 1)
	{
		read_settings(config, "settings.cfg");
	}
	else if (argc == 2)
	{
		read_settings(config, argv[1]);
	}
	else
	{
		std::cerr << "Usage: " << argv[0] << " [settings file]" << std::endl;
		exit(1);
	}
	unsigned long seed;
	if (config.use_seed)
	{
		seed = config.seed;
	}
	else
	{
		seed = time(NULL);
	}
	std::cout << "Seed: " << seed << std::endl;
	srand(seed);
	vector origin = vector(0, 0, 0);
	quadtree* root = new quadtree(&origin, config.size);
	std::vector<particle*> particles;
#ifdef THREADED
	pthread_t threads[NUM_THREADS];
	struct thread_data td[NUM_THREADS];
	int rc;
#endif
	std::cout << "Generating particles..." << std::endl;
	for (unsigned int i = 0; i < config.num_particles; i++)
	{
		generate_particle(config, particles, root);
	}
	for (unsigned int i = 0; i < config.num_frames; i++)
	{
		std::cout << "Frame " << i + 1 << "/" << config.num_frames << std::endl;
		root -> calc_mass();
		root -> calc_com();
#ifndef THREADED
		barnes_hut(particles, root, config.theta, config.display_progress);
#endif
#ifdef THREADED
		for (unsigned int i = 0; i < NUM_THREADS; i++)
		{
			td[i].thread_id = i;
			td[i].particles = &particles;
			td[i].root = root;
			td[i].start = (config.num_particles / NUM_THREADS) * i;
			td[i].end = (config.num_particles / NUM_THREADS) * (i + 1);
			td[i].theta = config.theta;
			rc = pthread_create(&threads[i], NULL, barnes_hut_thread, (void*) &td[i]);
			if (rc)
			{
				std::cout << "Could not create thread." << std::endl;
				exit(-1);
			}
		}
#endif
		update_all(particles, config.dt);
		check_tree(particles, root);
		delete root;
		root = new quadtree(&origin, config.size);
		for (unsigned int i = 0; i < particles.size(); i++)
		{
			root -> add_particle(particles[i]);
		}
		if (config.dump_binary)
		{
			dump(i, particles, config.overwrite_data);
		}
		if (config.dump_plaintext)
		{
			dump_plaintext(i, particles, config.overwrite_data);
		}
	}
	delete root;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		delete particles[i];
	}
#ifdef THREADED
	pthread_exit(NULL);
#endif
	return 0;
}
