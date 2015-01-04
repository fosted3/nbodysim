#define CUBE 0
#define SPHERE 1
#define SHELL 2
#define LINEAR 0
#define EXP 1
#define NORMAL 2
#define FRONT 0
#define SIDE 1
#define TOP 2
#define ISO 3
#define BW 0
#define HEAT 1

#include "thread_functions.h"
#include "vector.h"
#include "particle.h"
#include "octree.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <queue>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <pthread.h>
#include <FreeImage.h>
#include <limits>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>

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

typedef std::pair<particle*, particle*> particle_pair;
typedef std::unordered_set<particle*> particle_set;
typedef std::unordered_set<std::pair<particle*, particle*> > particle_pair_set;
typedef std::unordered_map<particle*, particle*> particle_map;

/****************************************\
| Structs for passing thread information |
\****************************************/

struct thread_data //Used for Barnes-Hut thread
{
	unsigned int thread_id;
	particle_set *particles;
	octree *root;
	unsigned int modulus;
	datatype theta;
	bool damping;
	bool print;
	unsigned int num_particles;
	particle_pair_set *collision_data;
	datatype range;
};

struct file_write_data //Used for writing data (binary, text, and image)
{
	particle_set *particles;
	unsigned int frame;
	unsigned int img_w;
	unsigned int img_h;
	unsigned int projection;
	unsigned int color;
	datatype scale;
	datatype brightness;
	bool binary;
	bool text;
	bool image;
	bool overwrite;
	bool keep_prev_binary;
	bool keep_prev_text;
	bool adaptive;
	bool nonlinear;

};

struct update_data //Used for updating particles
{
	particle_set *particles;
	unsigned int start;
	unsigned int end;
	datatype dt;
};

struct generate_data //Used for generating root
{
	octree *target;
	particle_set *particles;
};

/********************************\
| Struct for all config settings |
\********************************/

struct settings
{
	bool read_existing;			//Read data from existing binary file
	bool use_seed;				//Use user specified seed
	bool display_progress;		//Display estimated percentage completion during Barnes-Hut calculations
	bool dump_binary;			//Write binary (resume) data to file
	bool dump_text;				//Write text (position) data to plaintext
	bool dump_image;			//Write image to file
	bool overwrite_data;		//If file exists, overwrite it
	bool keep_previous_binary;	//Keep previous binary file after writing
	bool keep_previous_text;	//Keep previous text file after writing
	bool verbose;				//Display extra information
	bool damping;				//Close interaction damping
	bool adaptive_brightness;	//Adaptive brightness
	bool nonlinear_brightness;	//Nonlinear brightness
	bool collide;				//Calculate collisions
	unsigned int num_particles; //How many particles to generate if no resume present
	unsigned int num_frames;	//How many frames to compute
	unsigned int threads;		//How many threads to use for Barnes-Hut calculation
	unsigned int img_w;			//Image width
	unsigned int img_h;			//Image height
	unsigned int projection;	//Image projection
	datatype size;				//Size of cube if cubic generation specified
	datatype theta;				//Theta value for Barnes-Hut calculations
	datatype dt;				//Time between frames
	unsigned long seed;			//User specified seed value
	datatype min_mass;			//Minimum particle mass
	datatype max_mass;			//Maximum particle mass
	datatype min_vel;			//Minimum particle velocity
	datatype max_vel;			//Maximum particle velocity
	datatype r_sphere;			//Radius of sphere if sphere / shell generation specified
	datatype rotation_magnitude;//Rotation of sphere / shell
	datatype scale_x;			//Scale particle x coordinate
	datatype scale_y;			//Scale particle y coordinate
	datatype scale_z;			//Scale particle z coordinate
	datatype brightness;		//How much white to add to location on image if particle lands in pixel
	datatype scale;				//Image scale (zoom)
	datatype collision_range;	//Distance at which particles collide
	vector rotation_vector;		//Vector of rotation (use in file as rotation_vector x y z)
	unsigned int gen_type;		//How simulation generates particles (cube, sphere, shell)
	unsigned int mass_dist;		//How mass is distributed (normal, exp, linear)
	unsigned int vel_dist;		//How velocity is distrubuted in cubic mode (normal, exp, linear)
	unsigned int color;			//Color mapping
};

/*******************\
| Utility functions |
\*******************/

namespace std
{
	template<> struct hash<particle_pair>
	{
		size_t operator() (const particle_pair &p) const
		{
			return std::hash<int>()(p.first -> get_pos() -> get_x()) ^ std::hash<int>()(p.first -> get_pos() -> get_y()) ^ std::hash<int>()(p.first -> get_pos() -> get_z()) ^ std::hash<int>()(p.second -> get_pos() -> get_x()) ^ std::hash<int>()(p.second -> get_pos() -> get_y()) ^ std::hash<int>()(p.second -> get_pos() -> get_z());
		}
	};
}

particle_set::iterator& operator+=(particle_set::iterator &itr, const unsigned int &inc)
{
	for (unsigned int x = 0; x < inc; x++)
	{
		*itr ++;
	}
	return itr;
}

datatype clamp(datatype a, datatype x, datatype b) //Make x such that a < x < b
{
	if (x < a) { return a; }
	if (x > b) { return b; }
	return x;
}

datatype random_datatype(datatype low, datatype high, unsigned int dist, std::default_random_engine &generator) //Random float between low and high with distribution dist
{
	//assert(dist == LINEAR || dist == EXP || dist == NORMAL);
	datatype r = -1;
	if (dist == LINEAR)
	{
		std::uniform_real_distribution<datatype> distribution(low, high);
		r = distribution(generator);
	}
	else if (dist == EXP)
	{
		datatype lambda = -1.0 * log(0.001) / (high - low); //Make it such that 99.9% of numbers lie from 0 to (high - low)
		std::exponential_distribution<datatype> distribution(lambda);
		do
		{
			r = distribution(generator) + low;
		} while (r > high || r < low); //Boundary check
	}
	else if (dist == NORMAL)
	{
		datatype center = (high + low) / 2.0;
		datatype deviation = (high - low) / 6.0; // +/- 3 std. devs
		std::normal_distribution<datatype> distribution(center, deviation);
		do
		{
			r = distribution(generator);
		} while (r > high || r < low); //Boundary check
	}
	else
	{
		std::cerr << "Invalid random distribution specified." << std::endl;
		exit(1);
	}
	return r;
}

vector random_vector(datatype low, datatype high, std::default_random_engine &generator) //Random vector, used for cubic generation
{
	vector rv = vector(random_datatype(low, high, LINEAR, generator), random_datatype(low, high, LINEAR, generator), random_datatype(low, high, LINEAR, generator));
	return rv;
}

void generate_particle(settings &s, particle_set *particles, std::default_random_engine &generator) //Generates a particle with the specified settings
{
	//assert(s.gen_type == SPHERE || s.gen_type == SHELL || s.gen_type == CUBE);
	vector null = vector(0, 0, 0); //Used to set acceleration to zero
	datatype mass = random_datatype(s.min_mass, s.max_mass, s.mass_dist, generator);
	assert(mass != -1);
	particle *par = NULL;
	assert(s.scale_x >= 0);
	assert(s.scale_y >= 0);
	assert(s.scale_z >= 0);
	if (s.gen_type == CUBE)
	{
		vector temp = random_vector(-1*(s.size/2), s.size/2, generator);
		temp.scale(s.scale_x, s.scale_y, s.scale_z);
		vector vel = random_vector(-1, 1, generator);
		vel.normalize();
		vel *= random_datatype(s.min_vel, s.max_vel, s.vel_dist, generator);
		par = new particle(&temp, &vel, &null, mass);
	}
	else if (s.gen_type == SPHERE || s.gen_type == SHELL)
	{
		datatype radius = -1;
		if (s.gen_type == SPHERE)
		{
			radius = random_datatype(0, 1, LINEAR, generator);
			radius = cbrt(radius) * s.r_sphere; //volume is proportional to the cube of the radius
		}
		else if (s.gen_type == SHELL)
		{
			radius = s.r_sphere;
		}
		assert(radius != -1);
		datatype theta = random_datatype(0, 6.28318530718, LINEAR, generator);
		datatype azimuth = acos(random_datatype(-1, 1, LINEAR, generator));
		vector point = vector(radius * sin(azimuth) * cos(theta), radius * sin(azimuth) * sin(theta), radius * cos(azimuth));
		point.scale(s.scale_x, s.scale_y, s.scale_z);
		vector rotation = s.rotation_vector;
		rotation.normalize();
		rotation *= s.rotation_magnitude;
		vector velocity = cross(rotation, point);
		par = new particle(&point, &velocity, &null, mass);
	}
	else
	{
		std::cerr << "Invalid generation type specified." << std::endl;
		exit(1);
	}
	assert(par != NULL);
	particles -> insert(par);
}

bool file_exists(const char *filename) //Check if a file exists
{
	std::ifstream ifile(filename);
	return ifile;
}

std::string gen_filename(unsigned int frame, bool binary) //Generate text / binary filename
{
	std::string filename = std::to_string(frame);
	while (filename.length() < 4) { filename.insert(0, "0"); }
	filename.insert(0, "./data/");
	if (binary) { filename += ".dat"; }
	else { filename += ".txt"; }
	return filename;
}

std::string gen_image(unsigned int frame) //Generate image filename
{
	std::string filename = std::to_string(frame);
	while (filename.length() < 4) { filename.insert(0, "0"); }
	filename.insert(0, "./img/");
	filename += ".png";
	return filename;
}

/*******************\
| Gravity functions |
\*******************/

vector gravity(particle* par, octree* node, bool damping) //Calculate the force acting on a particle from a node
{
	if (par == node -> get_particle()) { return vector(0, 0, 0); } //Particles exert no force on themselves
	vector acc = *(node -> get_com()); //This vector is multipurpose
	acc -= *(par -> get_pos()); //Here it is used as vector fromp particle to node (radius)
	datatype r_sq;
	if (!damping)
	{
		r_sq = pow(acc.magnitude(), -2); // 1/r^2 no damping
	}
	else
	{
		r_sq = 1.0 / (pow(acc.magnitude(), 2) + exp(-1 * acc.magnitude())); //1 / (r^2 + exp(-r)) damped
	}
	acc.normalize(); //Unit vector indicating direction of force
	acc *= node -> get_mass(); //Don't need particle mass, f=ma=GmM/r^2 cancels this out
	acc *= 6.67384e-11; //Newton's gravitational constant
	acc *= r_sq;
	return acc;
}

void *barnes_hut_thread(void *data) //Thread that calculates Barnes-Hut algorithm on a subset of particles
{
	struct thread_data *args;
	args = (struct thread_data*) data;
	std::queue<octree*> nodes;
	octree* node;
	octree* root = args -> root;
	particle* curr;
	particle_set *particles = args -> particles;
	assert(particles != NULL);
	vector grav_to;
	datatype theta = args -> theta;
	bool damping = args -> damping;
	datatype percent;
	unsigned int completed = 0;
	double collide_distance = args -> range;
	particle_pair_set *collision_data = args -> collision_data;
	particle_set::iterator itr = particles -> begin();
	itr += args -> thread_id;
	for (unsigned int pos = args -> thread_id; pos < particles -> size(); pos += args -> modulus)
	{
		if (args -> thread_id == 0 && args -> print && (pos - args -> thread_id) / (args -> modulus) % 25 == 0) //Thread 0 displays its progress because mutex locks
		{
			completed += 25 * (args -> modulus);
			percent = completed * 100;
			percent /= args ->  num_particles;
			printf("\b\b\b\b\b\b\b%3.2f%%", percent);
		}
		curr = *itr;
		curr -> set_acc_zero();
		nodes.push(root);
		while (!nodes.empty()) //Read wikipedia if you want to know the details of how this works
		{
			node = nodes.front();
			assert(node != NULL);
			assert(curr != NULL);
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
				if (collision_data != NULL && node -> get_particle() != NULL && distance(node -> get_particle() -> get_pos(), curr -> get_pos()) < collide_distance && node -> get_particle() != curr)
				{
					if (node -> get_particle() -> get_pos() -> get_x() < curr -> get_pos() -> get_x())
					{
						collision_data -> insert(std::make_pair(node -> get_particle(), curr));
					}
					else
					{
						collision_data -> insert(std::make_pair(curr, node -> get_particle()));
					}
				}
				grav_to = gravity(curr, node, damping);
				curr -> set_acc_offset(&grav_to);
			}
		}
		for (unsigned int i = 0; i < args -> modulus && itr != particles -> end(); i++)
		{
			itr ++;
		}
		if (itr == particles -> end()) { break; }
	}
	pthread_exit(NULL);
}

particle* collide(particle *a, particle *b)
{
	vector pos = weighted_average(a -> get_pos(), b -> get_pos(), a -> get_mass(), b -> get_mass());
	vector vel = weighted_average(a -> get_vel(), b -> get_vel(), a -> get_mass(), b -> get_mass());
	vector null = vector(0, 0, 0);
	datatype mass = a -> get_mass() + b -> get_mass();
	return new particle(&pos, &vel, &null, mass);
}
	
void barnes_hut_threaded(struct settings &config, particle_set *particles, octree *root, particle_set *added, particle_set *removed)
{
	pthread_t *threads = new pthread_t[config.threads];
	struct thread_data *td = new thread_data[config.threads];
	particle_pair_set collision_data;
	particle_map update_table;
	unsigned long count;
	for (unsigned int i = 0; i < config.threads; i++)
	{
		td[i].thread_id = i;
		td[i].particles = particles;
		td[i].root = root;
		td[i].modulus = config.threads;
		td[i].theta = config.theta;
		td[i].damping = config.damping;
		td[i].print = config.display_progress;
		td[i].num_particles = particles -> size();
		td[i].range = config.collision_range;
		if (added != NULL && removed != NULL) { td[i].collision_data = new particle_pair_set; }
		else { td[i].collision_data = NULL; }
		create_thread(&threads[i], NULL, barnes_hut_thread, (void*) &td[i]);
	}
	for (unsigned int i = 0; i < config.threads; i++)
	{
		pthread_join(threads[i], NULL);
		if (config.display_progress && i == 0) //Get rid of that last print statement
		{
			printf("\b\b\b\b\b\b\b\b");
		}
		if (added != NULL && removed != NULL)
		{
			count = 0;
			collision_data.reserve(collision_data.size() + td[i].collision_data -> size());
			for (particle_pair_set::const_iterator itr = td[i].collision_data -> begin(); itr != td[i].collision_data -> end(); itr++)
			{
				assert(&(*itr) != NULL);
				collision_data.insert(*itr);
				if (config.display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu/%lu", count, td[i].collision_data -> size()); }
				count++;
			}
			if (config.display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); }
		}
	}
	if (added != NULL && removed != NULL)
	{
		for (unsigned int i = 0; i < config.threads; i++)
		{
			delete td[i].collision_data;
		}
		removed -> reserve(collision_data.size() * 2);
		added -> reserve(collision_data.size());
		count = 0;
		if (config.display_progress) { std::cout << "Culling collision data..." << std::endl; }
		for (particle_pair_set::const_iterator collision_iter = collision_data.begin(); collision_iter != collision_data.end(); collision_iter++)
		{
			particle *a = collision_iter -> first;
			particle *b = collision_iter -> second;
			particle_map::const_iterator first_find = update_table.find(a);
			particle_map::const_iterator second_find = update_table.find(b);
			while (first_find != update_table.end())
			{
				a = first_find -> second;
				first_find = update_table.find(a);
			
			}
			while (second_find != update_table.end())
			{
				b = second_find -> second;
				second_find = update_table.find(b);
			}
			if (a == b) { continue; }
			assert(a != NULL);
			assert(b != NULL);
			particle* newp = collide(a, b);
			update_table[a] = newp;
			update_table[b] = newp;
			added -> erase(a);
			added -> erase(b);
			added -> insert(newp);
			removed -> insert(a);
			removed -> insert(b);
			count ++;
			if (config.display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu/%lu", count, collision_data.size()); }
		}
		if (config.display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); }
	}
	delete[] threads; //Memory management
	delete[] td;
}

void update_collision(particle_set *particles, particle_set *added, particle_set *removed, bool display_progress)
{
	unsigned long count = 0;
	particle_set::iterator particle_itr;
	particle_set::const_iterator particle_find;
	particle_set::iterator added_itr;
	particle_set::iterator removed_itr;
	if (display_progress) { std::cout << "Processing removed..." << std::endl; } 
	for (removed_itr = removed -> begin(); removed_itr != removed -> end(); removed_itr++)
	{
		particle_find = particles -> find(*removed_itr);
		if (particle_find != particles -> end())
		{
			particles -> erase(*removed_itr);
		}
		delete *removed_itr;
		count ++;
		if (display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu/%lu", count, removed -> size()); }
	}
	count = 0;
	if (display_progress) { std::cout << "Processing added..." << std::endl; }
	particles -> reserve(particles -> size() + added -> size());
	for (added_itr = added -> begin(); added_itr != added -> end(); added_itr++)
	{
		count++;
		particles -> insert(*added_itr);
		if (display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu/%lu", count, added -> size()); }
	}
	if (display_progress) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); }
}

/****************\
| Data functions |
\****************/

void write_image(unsigned int img_w, unsigned int img_h, unsigned int projection, unsigned int color, bool adaptive, bool nonlinear, datatype scale, datatype brightness, unsigned int frame, particle_set *particles)
{
	FIBITMAP *image = FreeImage_Allocate(img_w, img_h, 24); //Image allocation, wxh, 24bpp
	RGBQUAD pixel; //Color variable
	if (!image)
	{
		std::cerr << "Can't allocate memory for image. Exiting." << std::endl;
		exit(1);
	}
	std::vector<std::vector<datatype> > temp(img_w, std::vector<datatype>(img_h)); //Temporary datatype array for more precise coloring
	for (unsigned int i = 0; i < img_w; i++)
	{
		for (unsigned int j = 0; j < img_h; j++)
		{
			temp[i][j] = 0;
		}
	}
	assert(projection == FRONT || projection == SIDE || projection == TOP || projection == ISO);
	datatype x = 0; //XY position of current particle
	datatype y = 0;
	int v = 0; //Value variable
	datatype max = 0;
	if (projection == ISO) { scale *= 2.0; } //Isometric projection shrinks scale by 2, compensate for this
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
	{
		if (projection == FRONT)
		{
			x = (*itr) -> get_pos() -> get_x() + (img_w / 2);
			y = (*itr) -> get_pos() -> get_z() + (img_h / 2);
		}
		else if (projection == SIDE)
		{
			x = (*itr) -> get_pos() -> get_y() + (img_w / 2);
			y = (*itr) -> get_pos() -> get_z() + (img_h / 2);
		}
		else if (projection == TOP)
		{
			x = (*itr) -> get_pos() -> get_x() + (img_w / 2);
			y = (*itr) -> get_pos() -> get_y() + (img_h / 2);
		}
		else
		{
			x = ((sqrt(3) / 2.0) * ((*itr) -> get_pos() -> get_x() - (*itr) -> get_pos() -> get_y()) + img_w) / 2.0;
			y = ((-0.5) * ((*itr) -> get_pos() -> get_x() + (*itr) -> get_pos() -> get_y() + 2 * (*itr) -> get_pos() -> get_z()) + img_h) / 2.0;
		}
		x += (x - (img_w / 2)) * (scale - 1); //Scaling
		y += (y - (img_h / 2)) * (scale - 1);
		if (x < 0 || y < 0 || x > img_w - 1 || y > img_w - 1) { continue; }
		x = clamp(0, x, img_w - 1); //Make sure it's inside the image
		y = clamp(0, y, img_h - 1);
		if (projection == ISO) { y = (img_h - 1) - y; } //I'm not quite sure why but the image flipped upside down in isometric
		temp[(int) x][(int) y] += brightness; //Store value into array
		if (temp[(int) x][(int) y] > max) { max = temp[(int) x][(int) y]; }
	}
	if (nonlinear)
	{
		max = sqrt(max);
		for (unsigned int i = 0; i < img_w; i++)
		{
			for (unsigned int j = 0; j < img_h; j++)
			{
				temp[i][j] = sqrt(temp[i][j]);
			}
		}
	}
	for (unsigned int x_i = 0; x_i < img_w; x_i++)
	{
		for (unsigned int y_i = 0; y_i < img_h; y_i++)
		{
			if (color == BW)
			{
				if (adaptive)
				{
					v = (int) clamp(0, temp[x_i][y_i] * 255.0 / max, 255);
				}
				else
				{
					v = (int) clamp(0, temp[x_i][y_i], 255); //Cast to int and clamp to sane values
				}
				pixel.rgbRed = v;
				pixel.rgbBlue = v;
				pixel.rgbGreen = v;
				FreeImage_SetPixelColor(image, x_i, y_i, &pixel); //Set pixel
			}
			else if (color == HEAT)
			{
				if (adaptive)
				{
					v = (int) clamp(0, temp[x_i][y_i] * 511.0 / max, 511);
				}
				else
				{
					v = (int) clamp(0, temp[x_i][y_i], 511);
				}
				if (v < 256)
				{
					pixel.rgbBlue = clamp(0, 255 - v, 255);
					pixel.rgbGreen = clamp(0, v, 255);
					pixel.rgbRed = 0;
				}
				else
				{
					v -= 256;
					pixel.rgbBlue = 0;
					pixel.rgbGreen = clamp(0, 255 - v, 255);
					pixel.rgbRed = clamp(0, v, 255);
				}
				FreeImage_SetPixelColor(image, x_i, y_i, &pixel);
			}
		}
	}
	if (!FreeImage_Save(FIF_PNG, image, gen_image(frame).c_str(), 0)) //Make sure the image is saved
	{
		std::cerr << "Cannot save " << gen_image(frame) << ". Exiting." << std::endl;
		exit(1);
	}
	FreeImage_Unload(image); //Deallocate memory from image
}

void *dump(void *data) //Data dumping thread
{
	struct file_write_data *args;
	args = (struct file_write_data*) data;
	unsigned int size = sizeof(particle);
	unsigned int frame = args -> frame; //Decompress arg pointers into local data for easier access
	unsigned int img_w = args -> img_w;
	unsigned int img_h = args -> img_h;
	unsigned int projection = args -> projection;
	unsigned int color = args -> color;
	bool overwrite = args -> overwrite;
	bool keep_prev_binary = args -> keep_prev_binary;
	bool keep_prev_text = args -> keep_prev_text;
	bool binary = args -> binary;
	bool text = args -> text;
	bool image = args -> image;
	bool adaptive = args -> adaptive;
	bool nonlinear = args -> nonlinear;
	datatype scale = args -> scale;
	datatype brightness = args -> brightness;
	particle_set *particles = args -> particles;
	std::string bfilename = gen_filename(frame, true);
	std::string tfilename = gen_filename(frame, false);
	std::string ifilename = gen_image(frame);
	if (!overwrite && file_exists(bfilename.c_str())) //Overwrite checks
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << bfilename << ". Aborting." << std::endl;
		exit(1);
	}
	else if (file_exists(bfilename.c_str())) { assert(remove(bfilename.c_str()) == 0); } //Delete if no overwrite
	if (!overwrite && file_exists(tfilename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << tfilename << ". Aborting." << std::endl;
		exit(1);
	}
	else if (file_exists(tfilename.c_str())) { assert(remove(tfilename.c_str()) == 0); }
	if (!overwrite && file_exists(ifilename.c_str()))
	{
		std::cerr << "Overwrite disabled, cannot overwrite " << ifilename << ". Aborting." << std::endl;
		exit(1);
	}
	else if (file_exists(ifilename.c_str())) { assert(remove(ifilename.c_str()) == 0); }
	if (!keep_prev_binary && frame > 0) //Delete last data if it exists
	{
		if (file_exists(gen_filename(frame - 1, true).c_str())) { assert(remove(gen_filename(frame - 1, true).c_str()) == 0); }
	}
	if (!keep_prev_text && frame > 0)
	{
		if (file_exists(gen_filename(frame - 1, false).c_str())) { assert(remove(gen_filename(frame - 1, false).c_str()) == 0); }
	}
	if (binary) //Dump binary data
	{
		std::fstream boutfile(bfilename, std::ios::out | std::ios::binary);
		boutfile.seekp(0);
		for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) { boutfile.write((char*)(*itr), size); }
		boutfile.flush();
		boutfile.close();
	}
	if (text) //Dump text data
	{
		vector temp;
		std::fstream toutfile(tfilename, std::ios::out);
		toutfile.seekp(0);
		for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
		{
			temp = *((*itr) -> get_pos());
			toutfile << temp.get_x() << ", " << temp.get_y() << ", " << temp.get_z() << "\n";
		}
		toutfile.flush();
		toutfile.close();
	}
	if (image) //Dump image
	{
		write_image(img_w, img_h, projection, color, adaptive, nonlinear, scale, brightness, frame, particles);
	}
	pthread_exit(NULL);
}

void dump_threaded(struct settings &config, unsigned int frame, particle_set *particles, pthread_t &file_thread) //Data dumping call
{
	struct file_write_data fd; //Set up args
	fd.particles = particles;
	fd.frame = frame;
	fd.binary = config.dump_binary;
	fd.text = config.dump_text;
	fd.overwrite = config.overwrite_data;
	fd.keep_prev_binary = config.keep_previous_binary;
	fd.keep_prev_text = config.keep_previous_text;
	fd.img_w = config.img_w;
	fd.img_h = config.img_h;
	fd.projection = config.projection;
	fd.scale = config.scale;
	fd.brightness = config.brightness;
	fd.image = config.dump_image;
	fd.color = config.color;
	fd.adaptive = config.adaptive_brightness;
	fd.nonlinear = config.nonlinear_brightness;
	create_thread(&file_thread, NULL, dump, (void*) &fd); //No thread joining here, thread is defined in main function and is joined there
}

unsigned int read_data(particle_set *particles, unsigned int num_frames) //Read data back into particles
{
	assert (particles -> size() == 0);
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
		particles -> insert(to_add);
	}
	return frame + 1;
}

/************************\
| Housekeeping functions |
\************************/

void *gen_root_thread(void *data) //Thread for generating root
{
	struct generate_data *args;
	args = (struct generate_data*) data;
	particle_set *particles = args -> particles;
	octree *target = args -> target;
	assert (target != NULL);
	unsigned int added = 0;
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
	{
		if (target -> inside(*itr)) //Make sure the particle is inside the thing you're adding it to
		{
			//std::cout << "Adding particle @ " << *itr << std::endl;
			target -> add_particle(*itr);
			added ++;
		}
	}
	pthread_exit(NULL);
}

octree* gen_root_threaded(particle_set *particles) //Root generation call
{
	datatype min_x = std::numeric_limits<datatype>::max(); //If your particles go beyond this then you have a problem
	datatype min_y = std::numeric_limits<datatype>::max();
	datatype min_z = std::numeric_limits<datatype>::max();
	datatype max_x = std::numeric_limits<datatype>::lowest(); //lowest is negative max, min is very close to zero
	datatype max_y = std::numeric_limits<datatype>::lowest();
	datatype max_z = std::numeric_limits<datatype>::lowest();
	datatype size;
	vector origin;
	vector temp;
	pthread_t threads[8];
	struct generate_data data[8];
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) //Find bounds of particles
	{
		temp = *((*itr) -> get_pos());
		if (temp.get_x() < min_x) { min_x = temp.get_x(); }
		if (temp.get_x() > max_x) { max_x = temp.get_x(); }
		if (temp.get_y() < min_y) { min_y = temp.get_y(); }
		if (temp.get_y() > max_y) { max_y = temp.get_y(); }
		if (temp.get_z() < min_z) { min_z = temp.get_z(); }
		if (temp.get_z() > max_z) { max_z = temp.get_z(); }
	}
	origin = vector((min_x + max_x) / 2.0, (min_y + max_y) / 2.0, (min_z + max_z) / 2.0);
	if ((max_x - min_x) > (max_y - min_y) && (max_x - min_x) > (max_z - min_z)) //Make sure the octree is cubic (this may be changed later)
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
	octree* root = new octree(&origin, size); //Allocate new root
	for (unsigned int i = 0; i < 8; i++) //Allocate the children & call new threads
	{
		root -> allocate_child(i);
		data[i].target = root -> get_child(i);
		data[i].particles = particles;
		create_thread(&threads[i], NULL, gen_root_thread, (void*) &data[i]);
	}
	for (unsigned int i = 0; i < 8; i++) //Wait for completion
	{
		pthread_join(threads[i], NULL);
	}
	return root;
}
octree* gen_root(particle_set *particles)
{
	datatype min_x = std::numeric_limits<datatype>::max(); //If your particles go beyond this then you have a problem
	datatype min_y = std::numeric_limits<datatype>::max();
	datatype min_z = std::numeric_limits<datatype>::max();
	datatype max_x = std::numeric_limits<datatype>::lowest(); //lowest is negative max, min is very close to zero
	datatype max_y = std::numeric_limits<datatype>::lowest();
	datatype max_z = std::numeric_limits<datatype>::lowest();
	datatype size;
	vector origin;
	vector temp;
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) //Find bounds of particles
	{
		temp = *((*itr) -> get_pos());
		if (temp.get_x() < min_x) { min_x = temp.get_x(); }
		if (temp.get_x() > max_x) { max_x = temp.get_x(); }
		if (temp.get_y() < min_y) { min_y = temp.get_y(); }
		if (temp.get_y() > max_y) { max_y = temp.get_y(); }
		if (temp.get_z() < min_z) { min_z = temp.get_z(); }
		if (temp.get_z() > max_z) { max_z = temp.get_z(); }
	}
	origin = vector((min_x + max_x) / 2.0, (min_y + max_y) / 2.0, (min_z + max_z) / 2.0);
	if ((max_x - min_x) > (max_y - min_y) && (max_x - min_x) > (max_z - min_z)) //Make sure the octree is cubic (this may be changed later)
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
	octree* root = new octree(&origin, size); //Allocate new root
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
	{
		root -> add_particle(*itr);
	}
	return root;
}

void update_all(struct settings &config, particle_set *particles)
{
	datatype dt = config.dt;
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) { (*itr) -> update(dt); }
}

void *update_all_thread(void *data) //Thread for updating particle positions & velocities
{
	struct update_data *args;
	args = (struct update_data*) data;
	particle_set *particles = args -> particles;
	particle_set::iterator itr = particles -> begin();
	itr += args -> start;
	particle_set::iterator end = particles -> begin();
	end += args -> end;
	for (; itr != end; itr++) { (*itr) -> update(args -> dt); }
	pthread_exit(NULL);
}

void update_all_threaded(struct settings &config, particle_set *particles) //Update call
{
	pthread_t *threads = new pthread_t[config.threads]; //Allocate threads & args
	struct update_data *ud = new update_data[config.threads];
	for (unsigned int i = 0; i < config.threads; i++) //Distribute work
	{
		ud[i].particles = particles;
		ud[i].start = (particles -> size() / config.threads) * i;
		if (i == config.threads - 1) { ud[i].end = particles -> size(); }
		else { ud[i].end = (particles -> size() / config.threads) * (i + 1); }
		ud[i].dt = config.dt;
		create_thread(&threads[i], NULL, update_all_thread, (void*) &ud[i]);
	}
	for (unsigned int i = 0; i < config.threads; i++) //Wait for completion
	{
		pthread_join(threads[i], NULL);
	}
	delete[] threads; //Memory management
	delete[] ud;
}

void *delete_root_thread(void* obj) //Thread for deleting root
{
	octree* target = (octree*) obj;
	delete target;
	pthread_exit(NULL);
}

void delete_root_threaded(octree *root) //Delete root call
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

/********************\
| Settings functions |
\********************/

void read_settings(settings &s, const char* sfile) //Read config file
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
			else if (var.compare("dump_image") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.dump_image = true; }
				else { s.dump_image = false; }
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
			else if (var.compare("adaptive_brightness") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.adaptive_brightness = true; }
				else { s.adaptive_brightness = false; }
			}
			else if (var.compare("nonlinear_brightness") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.nonlinear_brightness = true; }
				else { s.nonlinear_brightness = false; }
			}
			else if (var.compare("collide") == 0)
			{
				cfg >> var;
				if (var.compare("true") == 0) { s.collide = true; }
				else { s.collide = false; }
			}
			else if (var.compare("rotation_vector") == 0)
			{
				datatype x;
				datatype y;
				datatype z;
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
				else if (var.compare("normal") == 0) { s.mass_dist = NORMAL; }
				assert(s.mass_dist != 128);
			}
			else if (var.compare("vel_dist") == 0)
			{
				cfg >> var;
				s.vel_dist = 128;
				if (var.compare("linear") == 0) { s.vel_dist = LINEAR; }
				else if (var.compare("exp") == 0) { s.vel_dist = EXP; }
				else if (var.compare("normal") == 0) { s.vel_dist = NORMAL; }
				assert(s.vel_dist != 128);
			}
			else if (var.compare("projection") == 0)
			{
				cfg >> var;
				s.projection = 128;
				if (var.compare("front") == 0) { s.projection = FRONT; }
				else if (var.compare("side") == 0) { s.projection = SIDE; }
				else if (var.compare("top") == 0) { s.projection = TOP; }
				else if (var.compare("iso") == 0) { s.projection = ISO; }
				assert(s.projection != 128);
			}
			else if (var.compare("color") == 0)
			{
				cfg >> var;
				s.color = 128;
				if (var.compare("bw") == 0) { s.color = BW; }
				else if (var.compare("heat") == 0) { s.color = HEAT; }
				assert(s.color != 128);
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
			else if (var.compare("brightness") == 0) { cfg >> s.brightness; }
			else if (var.compare("img_w") == 0) { cfg >> s.img_w; }
			else if (var.compare("img_h") == 0) { cfg >> s.img_h; }
			else if (var.compare("scale") == 0) { cfg >> s.scale; }
			else if (var.compare("collision_range") == 0) { cfg >> s.collision_range; }
			else if (var.compare("false") == 0) { continue; }
			else if (var.compare("true") == 0) { continue; }
			else { std::cerr << "Unrecognized variable " << var << std::endl; }
		}
	}
	else
	{
		std::cerr << sfile << " doesn't exist or can't be read." << std::endl;
		exit(1);
	}
}

void set_default(settings &s) //Set settings to default values
{
	s.read_existing = false;
	s.use_seed = false;
	s.display_progress = true;
	s.dump_binary = true;
	s.dump_text = true;
	s.dump_image = true;
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
	s.brightness = 255;
	s.projection = ISO;
	s.img_w = 1920;
	s.img_h = 1080;
	s.scale = 1;
	s.color = BW;
	s.adaptive_brightness = false;
	s.nonlinear_brightness = false;
	s.collide = false;
	s.collision_range = 0.01;
}

int main(int argc, char **argv)
{
	unsigned int start_frame = 0; //Set some values, read settings
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
	unsigned long seed = 0; //Set up seed
	if (config.use_seed)
	{
		seed = config.seed;
	}
	else
	{
		seed = time(NULL);
	}
	std::default_random_engine generator;
	generator.seed(seed);
	octree* root;
	particle_set *particles = new particle_set;
	particle_set *added = NULL;
	particle_set *removed = NULL;
	pthread_t file_thread;
	if (config.read_existing) //Read existing data
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
	if (!config.read_existing || start_frame == 0) //Generate particles if not resuming or no data
	{
		std::cout << "Seed: " << seed << std::endl;
		std::cout << "Generating particles..." << std::endl;
		for (unsigned int i = 0; i < config.num_particles; i++)
		{
			generate_particle(config, particles, generator);
		}
	}
	frame = start_frame;
	std::cout << "Frame " << frame << "/" << config.num_frames << std::endl;
	assert(config.threads > 0);
	if (config.dump_image) { FreeImage_Initialise(); } //Init freeimage
	while (frame < config.num_frames) //Main loop
	{
		if (config.verbose) { std::cout << "Generating root..." << std::endl; }
		if (config.threads > 1) { root = gen_root_threaded(particles); }
		else { root = gen_root(particles); }
		//root = gen_root(particles);
		if (config.verbose) { std::cout << "Calculating masses of nodes..." << std::endl; }
		if (config.threads > 1) { root -> calc_mass_threaded(); }
		else { root -> calc_mass(); }
		if (config.verbose) { root -> print_info(); }
		if (config.verbose) { std::cout << "Calculating COMs of nodes..." << std::endl; }
		if (config.threads > 1) { root -> calc_com_threaded(); }
		else { root -> calc_com(); }
		if (config.verbose) { std::cout << "Starting Barnes-Hut algorithm..." << std::endl; }
		if (config.collide)
		{
			added = new particle_set;
			removed = new particle_set;
		}
		barnes_hut_threaded(config, particles, root, added, removed);
		if (config.verbose) { std::cout << "Deallocating nodes..." << std::endl; }
		if (config.threads > 1) { delete_root_threaded(root); }
		else { delete root; }
		if (first) //pthread_join can't be called on uninitialized value file_thread
		{
			first = false;
		}
		else if ((config.dump_binary || config.dump_text || config.dump_image) && config.threads > 1)
		{
			pthread_join(file_thread, NULL);
		}
		if (config.collide)
		{
			if (config.verbose) { std::cout << "Calculating collision results..." << std::endl; }
			update_collision(particles, added, removed, config.display_progress);
			delete added;
			delete removed;
		}
		if (config.verbose) { std::cout << "Number of particles: " << particles -> size() << std::endl; }
		if (config.verbose) { std::cout << "Updating particles..." << std::endl; }
		if (config.threads > 1) { update_all_threaded(config, particles); }
		else { update_all(config, particles); }
		if (config.dump_binary || config.dump_text || config.dump_image)
		{
			if (config.verbose) { std::cout << "Dumping data..." << std::endl; }
			dump_threaded(config, frame, particles, file_thread);
			if (config.threads == 1) { pthread_join(file_thread, NULL); }
		}
		frame++;
		std::cout << "Frame " << frame << "/" << config.num_frames << std::endl;
	}
	if (!first && (config.dump_binary || config.dump_text || config.dump_image)) //there's data to be written possibly
	{
		pthread_join(file_thread, NULL);
	}
	if (config.dump_image) { FreeImage_DeInitialise(); } //Deinit freeimage
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) //Deallocate particles
	{
		delete *itr;
	}
	delete particles;
	return 0;
}
