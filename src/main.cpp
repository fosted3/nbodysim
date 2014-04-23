#include "vector.h"
#include "particle.h"
#include "quadtree.h"
#include <cstddef>

int main()
{
	vector origin = vector(0, 0, 0);
	unsigned long size = 512;
	quadtree* root = new quadtree(&origin, size);
	for (int i = 0; i < 8; i++)
	{
		root -> allocate_child(i);
	}
	delete root;
	return 0;
}
