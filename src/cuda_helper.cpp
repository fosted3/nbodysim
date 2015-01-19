#include "particle.h"
#include "octree.h"
#include "thread_functions.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "vector_types.h"
#include <cassert>
#include "cuda_code.h"
#include "cuda_types.h"
#include <iostream>
#include "data_structures.h"

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

struct cuda_thread_data
{
	cudaStream_t *stream;
	particle_set *particles;
	octree *root;
	unsigned int thread_id;
	unsigned int modulus;
	datatype theta;
	bool damping;
};

void *barnes_hut_cuda_thread(void *data)
{
	struct cuda_thread_data *args;
	args = (struct cuda_thread_data*) data;
	cudaStream_t *stream = args -> stream;
	particle_set *particles = args -> particles;
	octree *root = args -> root;
	datatype theta = args -> theta;
	octree *node;
	std::queue<octree*> nodes;
	bool damping = args -> damping;
	particle_set::iterator particle_itr = particles -> begin();
	for (unsigned int i = 0; i < args -> thread_id && particle_itr != particles -> end(); i++)
	{
		particle_itr ++;
	}
	uint16_t dependancy_index;
	cparticle *host_par = allocate_host_particles();
	cparticle *dev_par = allocate_dev_particles();
	cnode *host_nodes = allocate_host_nodes();
	cnode *dev_nodes = allocate_dev_nodes();
	datatype3 *host_res = allocate_host_results();
	datatype3 *dev_res = allocate_dev_results();
	particle **pars = new particle*[block_size];
	uint32_t num_nodes = 0;
	uint16_t par_loc = 0;
	vector temp_vec;
	datatype percent;
	uint64_t completed = 0;
	uint16_t thread_count = 0;
	bool is_cont = false;
	std::queue<octree*> cont_nodes;
	for (unsigned int pos = args -> thread_id; pos < particles -> size(); pos += compute_threads)
	{
		dependancy_index = 0;
		//std::cout << "Thread " << args -> thread_id << " at pos " << pos << std::endl;
		if (args -> thread_id == 0 && (pos - args -> thread_id) / (args -> modulus) % 100 == 0) //Thread 0 displays its progress because mutex locks
		{
			completed += 100 * (args -> modulus);
			percent = completed * 100;
			percent /= args ->  particles -> size();
			printf("\b\b\b\b\b\b\b%3.2f%%", percent);
		}
		assert(par_loc < block_size);
		pars[par_loc] = *particle_itr;
		pars[par_loc] -> set_acc_zero();
		host_par[par_loc].pos.x = pars[par_loc] -> get_pos() -> get_x();
		host_par[par_loc].pos.y = pars[par_loc] -> get_pos() -> get_y();
		host_par[par_loc].pos.z = pars[par_loc] -> get_pos() -> get_z();
		host_par[par_loc].size = 0;
		if (is_cont)
		{
			while(!cont_nodes.empty())
			{
				nodes.push(cont_nodes.front());
				cont_nodes.pop();
			}
			is_cont = false;
		}
		else
		{
			nodes.push(root);
		}
		while (!nodes.empty())
		{
			node = nodes.front();
			assert(node != NULL);
			nodes.pop();
			if (node -> get_side() / distance(node -> get_com(), (*particle_itr) -> get_pos()) > theta && node -> get_particle() == NULL)
			{
				for (unsigned int i = 0; i < 8; i++)
				{
					if (node -> get_child(i) != NULL) { nodes.push(node -> get_child(i)); }
				}
			}
			else if ((*particle_itr) != node -> get_particle()) //put node in host_cnodes
			{
				assert(num_nodes < shared_size * block_size);
				assert(par_loc < block_size);
				host_nodes[num_nodes].pos.x = node -> get_com() -> get_x();
				host_nodes[num_nodes].pos.y = node -> get_com() -> get_y();
				host_nodes[num_nodes].pos.z = node -> get_com() -> get_z();
				host_nodes[num_nodes].mass = node -> get_mass();
				assert(dependancy_index < shared_size);
				host_par[par_loc].dependants[dependancy_index] = num_nodes;
				dependancy_index++;
				//assert(dependancy_index < shared_size);
				host_par[par_loc].size = dependancy_index;
				if (dependancy_index == shared_size)
				{
					//std::cout << "Shared overflow, generating new unit" << std::endl;
					par_loc++;
					//assert(par_loc < block_size);
					if (par_loc == block_size && !nodes.empty())
					{
						while(!nodes.empty())
						{
							cont_nodes.push(nodes.front());
							nodes.pop();
						}
						is_cont = true;
					}
					else if (par_loc == block_size)
					{
						par_loc --;
					}
					else
					{
						assert(par_loc < block_size);
						assert(par_loc > 0);
						pars[par_loc] = *particle_itr;
						host_par[par_loc].pos.x = host_par[par_loc - 1].pos.x;
						host_par[par_loc].pos.y = host_par[par_loc - 1].pos.y;
						host_par[par_loc].pos.z = host_par[par_loc - 1].pos.z;
						host_par[par_loc].size = 0;
						dependancy_index = 0;
					}
				}
				num_nodes++;
			}
		}
		if (dependancy_index > thread_count) { thread_count = dependancy_index; }
		if (!is_cont)
		{
			par_loc++;
			for (unsigned int i = 0; i < args -> modulus && particle_itr != particles -> end(); i++)
			{
				particle_itr ++;
			}
		}
		if (par_loc == block_size || particle_itr == particles -> end())
		{
			if (thread_count % 32 && shared_size > 16)
			{
				thread_count += 32 - (thread_count % 32);
			}
			if (shared_size > 16) { assert(thread_count % 32 == 0); }
			run_compute(host_par, dev_par, host_nodes, dev_nodes, host_res, dev_res, par_loc, num_nodes, thread_count, damping, stream);
			//std::cout << par_loc << " " << num_nodes << std::endl;
			for (unsigned int i = 0; i < par_loc; i++)
			{
				temp_vec = vector(host_res[i].x, host_res[i].y, host_res[i].z);
				pars[i] -> set_acc_offset(&temp_vec);
			}
			num_nodes = 0;
			par_loc = 0;
			thread_count = 0;
		}
		if (particle_itr == particles -> end())
		{
			break;
		}
	}
	//std::cout << "Thread " << args -> thread_id << " exiting." << std::endl;
	free_host_particles(host_par);
	free_host_nodes(host_nodes);
	free_host_results(host_res);
	delete[] pars;
	free_dev_particles(dev_par);
	free_dev_nodes(dev_nodes);
	free_dev_results(dev_res);
	pthread_exit(NULL);
}

void barnes_hut_cuda(particle_set *particles, octree *root, bool damping, datatype theta)
{
	pthread_t *threads = new pthread_t[compute_threads];
	struct cuda_thread_data *td = new cuda_thread_data[compute_threads];
	cudaStream_t *streams = new cudaStream_t[compute_threads];
	init_streams(streams);
	//std::cout << "GPU data address is " << cache_addr << std::endl;
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		td[i].stream = &(streams[i]);
		td[i].particles = particles;
		td[i].root = root;
		td[i].thread_id = i;
		td[i].modulus = compute_threads;
		td[i].damping = damping;
		td[i].theta = theta;
		create_thread(&threads[i], NULL, barnes_hut_cuda_thread, (void*) &td[i]);
	}
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	printf("\b\b\b\b\b\b\b\b");
	delete[] threads;
	delete[] td;
	free_streams(streams);
	delete[] streams;
	//call_dev_reset();
}
