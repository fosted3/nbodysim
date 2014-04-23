#ifndef quadtree_h_
#define quadtree_h_

#include "vector.h"
#include "particle.h"

class particle;

class quadtree
{
	public:
		quadtree();
		quadtree(vector*, unsigned long, quadtree*); //center, side, parent
		quadtree(vector*, unsigned long);
		~quadtree();
		void add_particle(particle*);
		void allocate_child(int);
		void print_info();
		void print_info(int);
		void update_mass();
		double get_mass();
	private:
		quadtree* parent;
		vector center;
		unsigned long side;
		quadtree *children[8];
		particle* p;
		double mass;
};

#endif
