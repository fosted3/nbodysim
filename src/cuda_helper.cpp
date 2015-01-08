#include "particle.h"
#include "octree.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "vector_types.h"
#include <cassert>
#include "cuda.h"
#include "cuda_code.h"
#include "cuda_types.h"
#include <iostream>

#ifdef DOUBLE
#ifndef datatype3
#define datatype3 double3
#endif
#ifndef datatype
#define datatype double
#endif
#endif
#ifdef FLOAT
#ifndef datatype3
#define datatype3 float3
#endif
#ifndef datatype
#define datatype float
#endif
#endif
#ifndef datatype3
#error "DOUBLE / FLOAT undefined, use -DDOUBLE or -DFLOAT"
#endif

typedef std::unordered_set<particle*> particle_set;

void resolve_dependancies(particle *par, cparticle *last, std::vector<cnode> &host_nodes, std::unordered_map<octree*, uint32_t> &indicies, octree *root)
{
	datatype theta = 1.0f;
	octree *node;
	std::queue<octree*> nodes;
	std::unordered_map<octree*, uint32_t>::const_iterator itr;
	nodes.push(root);
	unsigned int pos = 0;
	uint32_t index;
	while (!nodes.empty())
	{
		node = nodes.front();
		assert(node != NULL);
		nodes.pop();
		if (node -> get_side() / distance(node -> get_com(), par -> get_pos()) > theta && node -> get_particle() == NULL)
		{
			for (unsigned int i = 0; i < 8; i++)
			{
				if (node -> get_child(i) != NULL) { nodes.push(node -> get_child(i)); }
			}
		}
		else if (par != node -> get_particle())
		{
			itr = indicies.find(node);
			if (itr == indicies.end())
			{
				index = host_nodes.size();
				indicies.insert(std::make_pair(node, index));
				host_nodes.push_back(cnode());
				host_nodes.back().pos.x = node -> get_com() -> get_x();
				host_nodes.back().pos.y = node -> get_com() -> get_y();
				host_nodes.back().pos.z = node -> get_com() -> get_z();
				host_nodes.back().mass = node -> get_mass();
			}
			else
			{
				index = itr -> second;
			}
			last -> dependants[pos] = index;
			last -> size = pos;
			pos ++;
			assert(pos < shared_size);
		}
	}
}

void barnes_hut_cuda(std::unordered_set<particle*> *particles, octree *root)
{
	std::queue<particle*> to_be_processed;
	std::unordered_map<octree*, uint32_t> indicies;
	std::vector<particle*> host_ptrs;
	std::vector<cparticle> host_particles;
 	std::vector<cnode> host_nodes;
 	std::vector<datatype3> host_results;
 	vector temp_vec;
	for (particle_set::iterator itr = particles -> begin(); itr != particles -> end(); itr++)
	{
		to_be_processed.push(*itr);
	}

	while (to_be_processed.size() > 0)
	{
		while (sizeof(cparticle) * host_particles.size() + sizeof(cnode) * host_nodes.size() + sizeof(datatype3) * host_results.size() < cuda_mem && to_be_processed.size() > 0)
		{
			host_particles.push_back(cparticle());
			host_results.push_back(datatype3());
			host_ptrs.push_back(to_be_processed.front());
			host_particles.back().pos.x = to_be_processed.front() -> get_pos() -> get_x();
			host_particles.back().pos.y = to_be_processed.front() -> get_pos() -> get_y();
			host_particles.back().pos.z = to_be_processed.front() -> get_pos() -> get_z();
			host_results.back().x = 0;
			host_results.back().y = 0;
			host_results.back().z = 0;
			resolve_dependancies(to_be_processed.front(), &host_particles.back(), host_nodes, indicies, root);
			to_be_processed.pop();
		}
		run_compute(&host_particles, &host_nodes, &host_results);
		for (unsigned int i = 0; i < host_results.size(); i++)
		{
			temp_vec = vector(host_results[i].x, host_results[i].y, host_results[i].z);
			host_ptrs[i] -> set_acc_offset(&temp_vec);
		}
		indicies.clear();
		host_ptrs.clear();
		host_particles.clear();
		host_nodes.clear();
		host_results.clear();
	}
}
