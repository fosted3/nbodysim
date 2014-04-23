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
		void add_particle(particle*);
		void allocate_child(int);
		void print_info();
		void print_info(int);
		void update_mass();
		double get_mass();
		void calc_com();
		vector* get_com();
		bool clean();
		bool inside(particle*);
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
