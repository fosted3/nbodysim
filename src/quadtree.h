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
		void print_info();
		void print_info(int);
		void calc_mass();
		double get_mass();
		void calc_com();
		vector* get_com();
		bool clean();
		bool inside(particle*);
		void release_particle();
		quadtree* get_parent();
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
