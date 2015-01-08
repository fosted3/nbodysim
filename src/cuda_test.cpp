#include "cuda_helper.h"
#include "vector.h"
#include "particle.h"
#include "octree.h"
#include <queue>
#include <unordered_map>
#include <unordered_set>

typedef std::unordered_set<particle*> particle_set;

int main()
{

	particle_set *particles = new particle_set;
	particle* temp_par;
	vector temp_pos;
	vector null_vec = vector(0, 0, 0);
	octree* root = new octree(&null_vec, 4);
	for (int x = -1; x < 2; x += 1)
	{
		for (int y = -1; y < 2; y += 1)
		{
			for (int z = -1; z < 2; z += 1)
			{
				if ((x | y | z) == 0) { continue; }
				temp_pos = vector(x, y, z);
				temp_par = new particle(&temp_pos, &null_vec, &null_vec, 1e10);
				particles -> insert(temp_par);
			}
		}
	}
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
	{
		root -> add_particle(*itr);
	}
	root -> calc_mass();
	root -> calc_com();
	root -> print_info(0);
	barnes_hut_cuda(particles, root);
	root -> print_info(0);
	delete root;
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++) //Deallocate particles
	{
		delete *itr;
	}
	delete particles;
	return 0;
}
