#include "vector.h"
#include "particle.h"
#include "quadtree.h"
#include <vector>

int main()
{
	vector origin = vector(0, 0, 0);
	vector unit_vector = vector(1, 1, 1);
	unsigned long size = 512;
	quadtree* root = new quadtree(&origin, size);
	std::vector<particle*> particles;
	particle a = particle(&origin, &origin, &unit_vector, 1.0, root);
	for (int i = 0; i < 8; i++)
	{
		root -> allocate_child(i);
	}
	delete root;
	return 0;
}
