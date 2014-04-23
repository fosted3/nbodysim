#ifndef quadtree_h_
#define quadtree_h_

#include "vector.h"
#include "particle.h"

class quadtree
{
	public:
		quadtree();
		quadtree(vector*, unsigned long, quadtree*); //center, side, parent
		quadtree(vector*, unsigned long);
		~quadtree();
		void allocate_child(int);
		void print_info();
		void print_info(int);
	private:
		quadtree* parent;
		vector center;
		unsigned long side;
		quadtree *children[8];
		particle* p;
};

#endif
