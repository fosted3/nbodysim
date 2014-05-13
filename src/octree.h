#ifndef octree_h_
#define octree_h_

#include "vector.h"
#include "particle.h"

class particle;

class octree
{
	public:
		octree();
		octree(vector*, double, octree*); //center, side, parent
		octree(vector*, double);
		~octree();
		void add_particle(particle*);
		void allocate_child(unsigned int);
		void print_info(void);
		void print_info(unsigned int);
		void calc_mass(void);
		double get_mass(void);
		void calc_com(void);
		vector* get_com(void);
		bool inside(particle*);
		void release_particle(void);
		octree* get_parent(void);
		double get_side(void);
		octree* get_child(unsigned int);
		particle* get_particle(void);
		void calc_mass_threaded(void);
		void calc_com_threaded(void);
		void release_child(unsigned int);
	private:
		octree* parent;
		vector center;
		vector com;
		double side;
		octree *children[8];
		particle* p;
		double mass;
};

void* calc_mass_thread(void*);
void* calc_com_thread(void*);


#endif
