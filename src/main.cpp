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
	particle* a = new particle(&origin, &origin, &origin, 1.0, root);
	particle* b = new particle(&unit_vector, &origin, &origin, 1.0, root);
	//particles.push_back(a);
	for (int i = 0; i < 8; i++)
	{
		root -> allocate_child(i);
	}
	root -> add_particle(a);
	root -> add_particle(b);
	//root -> clean();
	delete root;
	delete a;
	return 0;
}
