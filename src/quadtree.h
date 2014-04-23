#ifndef quadtree_h_
#define quadtree_h_

#include "vector.h"
#include "particle.h"

class quadtree
{
	public:
		quadtree();
		quadtree(vector*, unsigned long, quadtree*); //center, side, parent
	private:
		quadtree* parent;
		vector center;
		unsigned long side;
		quadtree** children;
		//particle* p;
};

#endif
