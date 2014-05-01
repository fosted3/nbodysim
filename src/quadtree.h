#ifndef quadtree_h_
#define quadtree_h_

#include "vector.h"
#include "particle.h"

class particle;

class quadtree
{
	public:
		quadtree();
		quadtree(vector*, double, quadtree*); //center, side, parent
		quadtree(vector*, double);
		~quadtree();
		bool add_particle(particle*);
		void allocate_child(int);
		void print_info(void);
		void print_info(int);
		void calc_mass(void);
		double get_mass(void);
		void calc_com(void);
		vector* get_com(void);
		bool inside(particle*);
		void release_particle(void);
		quadtree* get_parent(void);
		double get_side(void);
		quadtree* get_child(int);
		particle* get_particle(void);
	private:
		quadtree* parent;
		vector center;
		vector com;
		double side;
		quadtree *children[8];
		particle* p;
		double mass;
};

#endif
