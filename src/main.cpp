#define CUBE 0
#define SPHERE 1
#define SHELL 2
#define LINEAR 0
#define EXP 1
#define NORMAL 2

#include "thread_functions.h"
#include "vector.h"
#include "particle.h"
#include "octree.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <queue>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <pthread.h>

struct thread_data
{
	int thread_id;
	std::vector<particle*> *particles;
	octree *root;
	unsigned int start;
	unsigned int end;
	double theta;
	bool damping;
	bool print;
	unsigned int num_particles;
};

struct settings
{
	bool read_existing;
	bool use_seed;
	bool display_progress;
	bool dump_binary;
	bool dump_text;
	bool overwrite_data;
	bool keep_previous_binary;
	bool keep_previous_text;
	bool verbose;
	bool damping;
	unsigned int num_particles;
	unsigned int num_frames;
	unsigned int threads;
	double size;
	double theta;
	double dt;
	unsigned long seed;
	double min_mass;
	double max_mass;
	double min_vel;
	double max_vel;
	double max_vel_change;
	double max_pos_change;
	double r_sphere;
	double rotation_magnitude;
	double scale_x;
	double scale_y;
	double scale_z;
	vector rotation_vector;
	unsigned int gen_type;
	unsigned int mass_dist;
	unsigned int vel_dist;
};

struct file_write_data
{
	std::vector<particle*> *particles;
	unsigned int frame;
	bool binary;
	bool text;
	bool overwrite;
	bool keep_prev_binary;
	bool keep_prev_text;
};

struct update_data
{
	std::vector<particle*> *particles;
	unsigned int start;
	unsigned int end;
	double dt;
};

struct generate_data
{
	octree *target;
	std::vector<particle*> *particles;
};

double clamp(double min, double x, double max)
{
	if (x < min)
	{
		return min;
	}
	if (x > max)
	{
		return max;
	}
	return x;
}

double random_double(double low, double high, unsigned int dist)
{
	assert(dist == LINEAR || dist == EXP);
	double r = -1;
	if (dist == LINEAR)
	{
		r = rand();
		r /= RAND_MAX;
		r *= high - low;
		r += low;
	}
	else if (dist == EXP)
	{
		do
		{
			r = rand();
			r /= RAND_MAX;
			r = log(r) * -0.5 * (high - low);
		} while (r > high || r < low);
	}
	else
	{
		std::cerr << "This shouldn't happen." << std::endl;
		exit(1);
	}
	return r;
}

vector random_vector(double low, double high)
{
	vector rv = vector(random_double(low, high, LINEAR), random_double(low, high, LINEAR), random_double(low, high, LINEAR));
	return rv;
}

void generate_particle(settings &s, std::vector<particle*> &particles)
{
	vector null = vector(0, 0, 0);
	double mass = random_double(s.min_mass, s.max_mass, s.mass_dist);
	assert(mass != -1);
	particle *par = NULL;
	assert(s.scale_x >= 0 && s.scale_x <= 1);
	assert(s.scale_y >= 0 && s.scale_y <= 1);
	assert(s.scale_z >= 0 && s.scale_z <= 1);
	if (s.gen_type == CUBE)
	{
		vector temp = random_vector(-1*(s.size/2), s.size/2);
		temp.scale(s.scale_x, s.scale_y, s.scale_z);
		vector vel = random_vector(-1, 1);
		vel.normalize();
		vel *= random_double(s.min_vel, s.max_vel, s.vel_dist);
		par = new particle(&temp, &vel, &null, mass);
	}
	else if (s.gen_type == SPHERE || s.gen_type == SHELL)
	{
		double radius = -1;
		if (s.gen_type == SPHERE)
		{
			radius = random_double(0, 1, LINEAR);
			radius = sqrt(radius) * s.r_sphere;
		}
		else if (s.gen_type == SHELL)
		{
			radius = s.r_sphere;
		}
		assert(radius != -1);
		double theta = random_double(0, 6.28318530718, LINEAR);
		double azimuth = acos(random_double(-1, 1, LINEAR));
		vector point = vector(radius * sin(azimuth) * cos(theta), radius * sin(azimuth) * sin(theta), radius * cos(azimuth));
		point.scale(s.scale_x, s.scale_y, s.scale_z);
		vector rotation = s.rotation_vector;
		rotation.normalize();
		rotation *= s.rotation_magnitude;
		vector velocity = cross(rotation, point);
		par = new particle(&point, &velocity, &null, mass);
	}
	assert(par != NULL);
	particles.push_back(par);
}

void check_tree(std::vector<particle*> &particles, octree *root)
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
	for (unsigned int i = 0; i < particles.size(); i++) { particles[i] -> update(dt); }
}

void *update_all(void *data)
{
	struct update_data *args;
	args = (struct update_data*) data;
	for (unsigned int i = args -> start; i < args -> end; i++) { (*(args -> particles))[i] -> update(args -> dt); }
	pthread_exit(NULL);
}

double update_all(std::vector<particle*> &particles, double max_vel_change, double max_pos_change)
{
	double max_vel = 0;
	double max_pos = 0;
	double temp_vel = -1;
	double temp_pos = -1;
	double dt = -1;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		particles[i] -> update(temp_vel, temp_pos);
		if (temp_vel > max_vel) { max_vel = temp_vel; }
		if (temp_pos > max_pos) { max_pos = temp_pos; }
	}
	if ((max_vel_change / max_vel) < (max_pos_change / max_pos)) { dt = (max_vel_change / max_vel); }
	else { dt = (max_pos_change / max_pos); }
	return dt;
}

vector gravity(particle* par, octree* node, bool damping)
{
	if (par == node -> get_particle()) { return vector(0, 0, 0); }
	vector acc = *(node -> get_com());
	acc -= *(par -> get_pos());
	double r_sq;
	if (!damping)
	{
		r_sq = pow(acc.magnitude(), -2);
	}
	else
	{
		r_sq = 1.0 / (pow(acc.magnitude(), 2) + exp(-1 * acc.magnitude()));
	}
	acc.normalize();
	acc *= node -> get_mass();
	acc *= 6.67384e-11;
	acc *= r_sq;
	return acc;
}

void *barnes_hut_thread(void *data)
{
	struct thread_data *args;
	args = (struct thread_data*) data;
	std::queue<octree*> nodes;
	octree* node;
	octree* root = args -> root;
	particle* curr;
	std::vector<particle*> *particles = args -> particles;
	vector grav_to;
	double theta = args -> theta;
	bool damping = args -> damping;
	double percent;
	unsigned int completed = 0;
	for (unsigned int i = (args -> start); i < (args -> end); i++)
	{
		if (args -> thread_id == 0 && args -> print && i % 10 == 0)
		{
			completed += 10 * particles -> size() / args -> end;
			percent = completed * 100;
			percent /= args ->  num_particles;
			printf("\b\b\b\b\b\b\b%3.2f%%", percent);
		}
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
					if (node -> get_child(i) != NULL) { nodes.push(node -> get_child(i)); }
				}
			}
			else
			{
				grav_to = gravity(curr, node, damping);
				curr -> set_acc_offset(&grav_to);
			}
		}
	}
	pthread_exit(NULL);
}

void barnes_hut_threaded(struct settings &config, std::vector<particle*> &particles, octree* root)
{
	pthread_t *threads = new pthread_t[config.threads];
	struct thread_data *td = new thread_data[config.threads];
	for (unsigned int i = 0; i < config.threads; i++)
	{
		td[i].thread_id = i;
		td[i].particles = &particles;
		td[i].root = root;
		td[i].start = (particles.size() / config.threads) * i;
		td[i].end = (particles.size() / config.threads) * (i + 1);
		td[i].theta = config.theta;
		td[i].damping = config.damping;
		td[i].print = config.display_progress;
		td[i].num_particles = particles.size();
		create_thread(&threads[i], NULL, barnes_hut_thread, (void*) &td[i]);
	}
	for (unsigned int i = 0; i < config.threads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	if (config.display_progress)
	{
		printf("\b\b\b\b\b\b\b");
	}
	delete[] threads;
	delete[] td;
}

bool file_exists(const char *filename)
{
	std::ifstream ifile(filename);
	return ifile;
}

std::string gen_filename(unsigned int frame, bool binary)
{
	std::string filename = std::to_string(frame);
	while (filename.length() < 4) { filename.insert(0, "0"); }
	filename.insert(0, "./data/");
	if (binary) { filename += ".dat"; }
	else { filename += ".txt"; }
	return filename;
}

void *dump(void *data)
{
	struct file_write_data *args;
	args = (struct file_write_data*) data;
	unsigned int size = sizeof(particle);
	unsigned int frame = args -> frame;
	bool overwrite = args -> overwrite;
	bool keep_prev_binary = args -> keep_prev_binary;
	bool keep_prev_text = args -> keep_prev_text;
	bool binary = args -> binary;
	bool text = args -> text;
	std::vector<particle*> *particles = args -> particles;
	std::string bfilename = gen_filename(frame, true);
	std::string tfilename = gen_filename(frame, false);
	if (!overwrite && file_exists(bfilename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << bfilename << ". Aborting." << std::endl;
		exit(1);
	}
	if (!overwrite && file_exists(tfilename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << tfilename << ". Aborting." << std::endl;
		exit(1);
	}
	if (file_exists(bfilename.c_str())) { assert(remove(bfilename.c_str()) == 0); }
	if (file_exists(tfilename.c_str())) { assert(remove(tfilename.c_str()) == 0); }
	if (!keep_prev_binary && frame > 0)
	{
		if (file_exists(gen_filename(frame - 1, true).c_str())) { assert(remove(gen_filename(frame - 1, true).c_str()) == 0); }
	}
	if (!keep_prev_text && frame > 0)
	{
		if (file_exists(gen_filename(frame - 1, false).c_str())) { assert(remove(gen_filename(frame - 1, false).c_str()) == 0); }
	}
	if (binary)
	{
		std::fstream boutfile(bfilename, std::ios::out | std::ios::binary);
		boutfile.seekp(0);
		for (unsigned int i = 0; i < particles -> size(); i++) { boutfile.write((char*)(*particles)[i], size); }
		boutfile.flush();
		boutfile.close();
	}
	if (text)
	{
		vector temp;
		std::fstream toutfile(tfilename, std::ios::out);
		toutfile.seekp(0);
		for (unsigned int i = 0; i < particles -> size(); i++)
		{
			temp = *((*particles)[i] -> get_pos());
			toutfile << temp.get_x() << ", " << temp.get_y() << ", " << temp.get_z() << "\n";
		}
		toutfile.flush();
		toutfile.close();
	}
	pthread_exit(NULL);
}

unsigned int read_data(std::vector<particle*> &particles, unsigned int num_frames)
{
	assert (particles.size() == 0);
	unsigned int size = sizeof(particle);
	int frame = num_frames;
	unsigned int num_particles;
	while (!file_exists(gen_filename(frame, true).c_str()) && frame >= 0) { frame--; }
	if (frame == -1) { return 0; }
	std::fstream infile(gen_filename(frame, true), std::ios::in | std::ios::binary);
	particle temp;
	particle *to_add;
	infile.seekg(0, infile.end);
	num_particles = infile.tellg() / size;
	infile.seekg(0, infile.beg);
	for (unsigned int i = 0; i < num_particles; i++)
	{
		infile.read((char*) &temp, size);
		to_add = new particle(temp);
		particles.push_back(to_add);
	}
	return frame + 1;
}

octree* gen_root(std::vector<particle*> &particles)
{
	double min_x =  1e6;
	double min_y =  1e6;
	double min_z =  1e6;
	double max_x = -1e6;
	double max_y = -1e6;
	double max_z = -1e6;
	double size;
	vector origin;
	vector temp;
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		temp = *(particles[i] -> get_pos());
		if (temp.get_x() < min_x) { min_x = temp.get_x(); }
		if (temp.get_x() > max_x) { max_x = temp.get_x(); }
		if (temp.get_y() < min_y) { min_y = temp.get_y(); }
		if (temp.get_y() > max_y) { max_y = temp.get_y(); }
		if (temp.get_z() < min_z) { min_z = temp.get_z(); }
		if (temp.get_z() > max_z) { max_z = temp.get_z(); }
	}
	origin = vector((min_x + max_x) / 2.0, (min_y + max_y) / 2.0, (min_z + max_z) / 2.0);
	if ((max_x - min_x) > (max_y - min_y) && (max_x - min_x) > (max_z - min_z))
	{
		size = (max_x - min_x) + 2.0;
	}
	else if ((max_y - min_y) > (max_x - min_x) && (max_y - min_y) > (max_z - min_z))
	{
		size = (max_y - min_y) + 2.0;
	}
	else
	{
		size = (max_z - min_z) + 2.0;
	}
	octree* root = new octree(&origin, size);
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		root -> add_particle(particles[i]);
	}
	return root;
}

void *gen_root_thread(void *data)
{
	struct generate_data *args;
	args = (struct generate_data*) data;
	std::vector<particle*> *particles = args -> particles;
	octree *target = args -> target;
	assert (target != NULL);
	unsigned int added = 0;
	for (unsigned int i = 0; i < particles -> size(); i++)
	{
		if (target -> inside((*particles)[i]))
		{
			target -> add_particle((*particles)[i]);
			added ++;
		}
	}
	if (added == 0)
	{
		delete target;
	}
	pthread_exit(NULL);
}

octree* gen_root_threaded(std::vector<particle*> &particles)
{
	double min_x =  1e6;
	double min_y =  1e6;
	double min_z =  1e6;
	double max_x = -1e6;
	double max_y = -1e6;
	double max_z = -1e6;
	double size;
	vector origin;
	vector temp;
	pthread_t threads[8];
	struct generate_data data[8];
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		temp = *(particles[i] -> get_pos());
		if (temp.get_x() < min_x) { min_x = temp.get_x(); }
		if (temp.get_x() > max_x) { max_x = temp.get_x(); }
		if (temp.get_y() < min_y) { min_y = temp.get_y(); }
		if (temp.get_y() > max_y) { max_y = temp.get_y(); }
		if (temp.get_z() < min_z) { min_z = temp.get_z(); }
		if (temp.get_z() > max_z) { max_z = temp.get_z(); }
	}
	origin = vector((min_x + max_x) / 2.0, (min_y + max_y) / 2.0, (min_z + max_z) / 2.0);
	if ((max_x - min_x) > (max_y - min_y) && (max_x - min_x) > (max_z - min_z))
	{
		size = (max_x - min_x) + 2.0;
	}
	else if ((max_y - min_y) > (max_x - min_x) && (max_y - min_y) > (max_z - min_z))
	{
		size = (max_y - min_y) + 2.0;
	}
	else
	{
		size = (max_z - min_z) + 2.0;
	}
	octree* root = new octree(&origin, size);
	for (unsigned int i = 0; i < 8; i++)
	{
		root -> allocate_child(i);
		data[i].target = root -> get_child(i);
		data[i].particles = &particles;
		create_thread(&threads[i], NULL, gen_root_thread, (void*) &data[i]);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		pthread_join(threads[i], NULL);
	}
	return root;
}

void update_all_threaded(struct settings &config, std::vector<particle*> &particles)
{
	pthread_t *threads = new pthread_t[config.threads];
	struct update_data *ud = new update_data[config.threads];
	for (unsigned int i = 0; i < config.threads; i++)
	{
		ud[i].particles = &particles;
		ud[i].start = (particles.size() / config.threads) * i;
		ud[i].end = (particles.size() / config.threads) * (i + 1);
		ud[i].dt = config.dt;
		create_thread(&threads[i], NULL, update_all, (void*) &ud[i]);
	}
	for (unsigned int i = 0; i < config.threads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	delete[] threads;
	delete[] ud;
}

void dump_threaded(struct settings &config, unsigned int frame, std::vector<particle*> &particles, pthread_t &file_thread)
{
	struct file_write_data fd; 
	fd.particles = &particles;
	fd.frame = frame;
	fd.binary = config.dump_binary;
	fd.text = config.dump_text;
	fd.overwrite = config.overwrite_data;
	fd.keep_prev_binary = config.keep_previous_binary;
	fd.keep_prev_text = config.keep_previous_text;
	create_thread(&file_thread, NULL, dump, (void*) &fd);
}

void *delete_root_thread(void* obj)
{
	octree* target = (octree*) obj;
	delete target;
	pthread_exit(NULL);
}

void delete_root_threaded(octree *root)
{
	pthread_t threads[8];
	octree* objs[8];
	for (unsigned int i = 0; i < 8; i++)
	{
		if (root -> get_child(i) == NULL) { continue; }
		objs[i] = root -> get_child(i);
		create_thread(&threads[i], NULL, delete_root_thread, (void*) objs[i]);
	}
	for (unsigned int i = 0; i < 8; i++)
	{
		if (root -> get_child(i) == NULL) { continue; }
		pthread_join(threads[i], NULL);
		root -> release_child(i);
	}
	delete root;
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
				if (var.compare("true") == 0) { s.read_existing = true; }
				else { s.read_existing = false; }
			}
			else if (var.compare("use_seed") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.use_seed = true; }
				else { s.use_seed = false; }
			}
			else if (var.compare("display_progress") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.display_progress = true; }
				else { s.display_progress = false; }
			}
			else if (var.compare("dump_binary") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.dump_binary = true; }
				else { s.dump_binary = false; }
			}
			else if (var.compare("dump_text") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.dump_text = true; }
				else { s.dump_text = false; }
			}
			else if (var.compare("overwrite_data") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.overwrite_data = true; }
				else { s.overwrite_data = false; }
			}
			else if (var.compare("keep_previous_binary") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.keep_previous_binary = true; }
				else { s.keep_previous_binary = false; }
			}
			else if (var.compare("keep_previous_text") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.keep_previous_text = true; }
				else { s.keep_previous_text = false; }
			}
			else if (var.compare("verbose") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.verbose = true; }
				else { s.verbose = false; }
			}
			else if (var.compare("damping") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.damping = true; }
				else { s.damping = false; }
			}
			else if (var.compare("rotation_vector") == 0)
			{
				double x;
				double y;
				double z;
				cfg >> x;
				cfg >> y;
				cfg >> z;
				s.rotation_vector = vector(x, y, z);
			}
			else if (var.compare("gen_type") == 0)
			{
				cfg >> var;
				s.gen_type = 128;
				if (var.compare("cube") == 0) { s.gen_type = CUBE; }
				else if (var.compare("sphere") == 0) { s.gen_type = SPHERE; }
				else if (var.compare("shell") == 0) { s.gen_type = SHELL; }
				assert(s.gen_type != 128);
			}
			else if (var.compare("mass_dist") == 0)
			{
				cfg >> var;
				s.mass_dist = 128;
				if (var.compare("linear") == 0) { s.mass_dist = LINEAR; }
				else if (var.compare("exp") == 0) { s.mass_dist = EXP; }
				assert(s.mass_dist != 128);
			}
			else if (var.compare("vel_dist") == 0)
			{
				cfg >> var;
				s.vel_dist = 128;
				if (var.compare("linear") == 0) { s.vel_dist = LINEAR; }
				else if (var.compare("exp") == 0) { s.vel_dist = EXP; }
				assert(s.vel_dist != 128);
			}
			else if (var.compare("threads") == 0) { cfg >> s.threads; }
			else if (var.compare("scale_x") == 0) { cfg >> s.scale_x; }
			else if (var.compare("scale_y") == 0) { cfg >> s.scale_y; }
			else if (var.compare("scale_z") == 0) { cfg >> s.scale_z; }
			else if (var.compare("r_sphere") == 0) { cfg >> s.r_sphere; }
			else if (var.compare("rotation_magnitude") == 0) { cfg >> s.rotation_magnitude; }
			else if (var.compare("num_particles") == 0) { cfg >> s.num_particles; }
			else if (var.compare("num_frames") == 0) { cfg >> s.num_frames; }
			else if (var.compare("size") == 0) { cfg >> s.size; }
			else if (var.compare("theta") == 0) { cfg >> s.theta; }
			else if (var.compare("dt") == 0) { cfg >> s.dt; }
			else if (var.compare("seed") == 0) { cfg >> s.seed; }
			else if (var.compare("min_mass") == 0) { cfg >> s.min_mass; }
			else if (var.compare("max_mass") == 0) { cfg >> s.max_mass; }
			else if (var.compare("min_vel") == 0) { cfg >> s.min_vel; }
			else if (var.compare("max_vel") == 0) { cfg >> s.max_vel; }
			else if (var.compare("max_vel_change") == 0) { cfg >> s.max_vel_change; }
			else if (var.compare("max_pos_change") == 0) { cfg >> s.max_pos_change; }
			else if (var.compare("brightness") == 0) { cfg >> var; }
			else if (var.compare("projection") == 0) { cfg >> var; }
			else if (var.compare("img_w") == 0) { cfg >> var; }
			else if (var.compare("img_h") == 0) { cfg >> var; }
			else if (var.compare("scale") == 0) { cfg >> var; }
			else { std::cerr << "Unrecognized variable " << var << std::endl; }
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
	s.dump_text = true;
	s.overwrite_data = false;
	s.keep_previous_binary = false;
	s.keep_previous_text = true;
	s.verbose = false;
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
	s.max_vel_change = 3;
	s.max_pos_change = 3;
	s.r_sphere = 100;
	s.rotation_magnitude = 0.1;
	s.rotation_vector = vector(0, 0, 1);
	s.gen_type = CUBE;
	s.mass_dist = LINEAR;
	s.vel_dist = LINEAR;
	s.scale_x = 1;
	s.scale_y = 1;
	s.scale_z = 1;
	s.damping = false;
	s.threads = 1;
}

int main(int argc, char **argv)
{
	unsigned int start_frame = 0;
	unsigned int frame = 0;
	bool first = true;
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
	unsigned long seed = 0;
	if (config.use_seed)
	{
		seed = config.seed;
	}
	else
	{
		seed = time(NULL);
	}
	srand(seed);
	octree* root;
	std::vector<particle*> particles;
	pthread_t file_thread;
	if (config.read_existing)
	{
		start_frame = read_data(particles, config.num_frames);
		if (start_frame == 0)
		{
			std::cout << "No data to resume from." << std::endl;
		}
		else
		{
			std::cout << "Resuming from " << start_frame << std::endl;
		}
	}
	if (!config.read_existing || start_frame == 0)
	{
		std::cout << "Seed: " << seed << std::endl;
		std::cout << "Generating particles..." << std::endl;
		for (unsigned int i = 0; i < config.num_particles; i++)
		{
			generate_particle(config, particles);
		}
	}
	frame = start_frame;
	std::cout << "Frame " << frame << "/" << config.num_frames << std::endl;
	assert(config.threads > 0);
	while (frame < config.num_frames)
	{
		if (config.verbose) { std::cout << "Generating root..." << std::endl; }
		root = gen_root_threaded(particles);
		if (config.verbose) { std::cout << "Calculating masses of nodes..." << std::endl; }
		root -> calc_mass_threaded();
		if (config.verbose) { std::cout << "Calculating COMs of nodes..." << std::endl; }
		root -> calc_com_threaded();
		if (config.verbose) { std::cout << "Starting Barnes-Hut algorithm..." << std::endl; }
		barnes_hut_threaded(config, particles, root);
		if (first)
		{
			first = false;
		}
		else
		{
			pthread_join(file_thread, NULL);
		}
		if (config.verbose) { std::cout << "Updating particles..." << std::endl; }
		update_all_threaded(config, particles);
		if (config.verbose) { std::cout << "Deallocating nodes..." << std::endl; }
		delete_root_threaded(root);
		if (config.dump_binary || config.dump_text)
		{
			if (config.verbose) { std::cout << "Dumping data..." << std::endl; }
			dump_threaded(config, frame, particles, file_thread);
		}
		frame++;
		std::cout << "Frame " << frame << "/" << config.num_frames << std::endl;
	}
	if (!first)
	{
		pthread_join(file_thread, NULL);
	}
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		delete particles[i];
	}
	return 0;
}
