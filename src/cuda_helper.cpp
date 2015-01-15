#include "particle.h"
#include "octree.h"
#include "thread_functions.h"
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
#include <mutex>
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

struct cache_data
{
	uint32_t pos;
	bool in_use[compute_threads];
	cache_data(uint32_t gpu_addr, uint16_t thread_id)
	{
		pos = gpu_addr;
		for (uint16_t i = 0; i < compute_threads; i++)
		{
			if (i == thread_id) { in_use[i] = true; }
			else { in_use[i] = false; }
		}
	}
};

typedef std::unordered_map<octree*, cache_data*> active_map;
typedef std::unordered_map<octree*, uint32_t> stale_map;

struct cuda_thread_data
{
	cudaStream_t *stream;
	std::mutex *lock;
	particle_set *particles;
	octree *root;
	active_map *active;
	stale_map *stale;
	unsigned int thread_id;
	unsigned int modulus;
	uint32_t *last_free;
	cnode *cache_addr;
	cparticle *par_addr;
	datatype3 *res_addr;
};

void *barnes_hut_cuda_thread(void *data)
{
	struct cuda_thread_data *args;
	args = (struct cuda_thread_data*) data;
	cudaStream_t *stream = args -> stream;
	std::mutex *lock = args -> lock;
	particle_set *particles = args -> particles;
	octree *root = args -> root;
	active_map *active = args -> active;
	stale_map *stale = args -> stale;
	datatype theta = 1.0f;
	octree *node;
	particle **par = new particle*[block_size];
	std::queue<octree*> nodes;
	active_map::iterator active_itr;
	stale_map::iterator stale_itr;
	particle_set::iterator particle_itr = particles -> begin();
	for (unsigned int i = 0; i < args -> thread_id && particle_itr != particles -> end(); i++)
	{
		particle_itr ++;
	}
	datatype3 *results = new datatype3[block_size];
	uint32_t gpu_index;
	uint16_t dependancy_index;
	cnode temp_cnode;
	cparticle *temp_cparticles = new cparticle[block_size];
	vector temp_vec;
	datatype percent;
	uint64_t completed = 0;
	uint16_t in_queue = 0;
	for (unsigned int pos = args -> thread_id; pos < particles -> size(); pos += compute_threads)
	{
		while (!(args -> lock -> try_lock()));
		//std::cout << "Thread " << args -> thread_id << " at pos " << pos << std::endl;
		if (args -> thread_id == 0 && (pos - args -> thread_id) / (args -> modulus) % 25 == 0) //Thread 0 displays its progress because mutex locks
		{
			completed += 25 * (args -> modulus);
			percent = completed * 100;
			percent /= args ->  particles -> size();
			printf("\b\b\b\b\b\b\b%3.2f%%", percent);
		}
		args -> lock -> unlock();
		par[in_queue] = *particle_itr;
		par[in_queue] -> set_acc_zero();
		temp_cparticles[in_queue].pos.x = par[in_queue] -> get_pos() -> get_x();
		temp_cparticles[in_queue].pos.y = par[in_queue] -> get_pos() -> get_y();
		temp_cparticles[in_queue].pos.z = par[in_queue] -> get_pos() -> get_z();
		temp_cparticles[in_queue].size = 0;
		dependancy_index = 0;
		nodes.push(root);
		while (!nodes.empty())
		{
			node = nodes.front();
			assert(node != NULL);
			nodes.pop();
			if (node -> get_side() / distance(node -> get_com(), par[in_queue] -> get_pos()) > theta && node -> get_particle() == NULL)
			{
				for (unsigned int i = 0; i < 8; i++)
				{
					if (node -> get_child(i) != NULL) { nodes.push(node -> get_child(i)); }
				}
			}
			else if (par[in_queue] != node -> get_particle()) //node should be in cache
			{
				while (!(args -> lock -> try_lock()));
				//std::cout << "Thread " << args -> thread_id << " checking if " << node << " is in cache." << std::endl;
				active_itr = active -> find(node);
				stale_itr = stale -> find(node);
				assert(active_itr == active -> end() || stale_itr == stale -> end());
				if (active_itr == active -> end() && stale_itr == stale -> end()) // Not in active or stale, copy to GPU (protip - make sure there's enough VRAM)
				{
					//std::cout << node << " is not in cache, inserting." << std::endl;
					temp_cnode.pos.x = node -> get_com() -> get_x();
					temp_cnode.pos.y = node -> get_com() -> get_y();
					temp_cnode.pos.z = node -> get_com() -> get_z();
					temp_cnode.mass = node -> get_mass();
					if (active -> size() + stale -> size() == cache_size) //cache is full, stick in a stale spot
					{
						//std::cout << "Cache full, inserting into stale data" << std::endl;
						gpu_index = stale -> begin() -> second;
						stale -> erase(stale -> begin());
					}
					else
					{
						//std::cout << "Cache has space, inserting at " << *(args -> last_free) << std::endl;
						gpu_index = *(args -> last_free);
						(*(args -> last_free))++;
					}
					copy_to_gpu(args -> cache_addr, gpu_index, stream, &temp_cnode); //copy to gpu returning index for node
					active -> insert(std::make_pair(node, new cache_data(gpu_index, args -> thread_id)));
				}
				else if (active_itr != active -> end()) //Being used by a node
				{
					//std::cout << node << " in cache already." << std::endl;
					active_itr -> second -> in_use[(args -> thread_id)] = true;
					gpu_index = active_itr -> second -> pos;
				}
				else if (stale_itr != stale -> end()) //Not being used by a node, move to active
				{
					//std::cout << node << " is stale, moving to active" << std::endl;
					gpu_index = stale_itr -> second;
					active -> insert(std::make_pair(node, new cache_data(gpu_index, args -> thread_id)));
					stale -> erase(stale_itr);
				}
				args -> lock -> unlock();
				assert(dependancy_index < shared_size);
				temp_cparticles[in_queue].dependants[dependancy_index] = gpu_index;
				dependancy_index ++;
				temp_cparticles[in_queue].size = dependancy_index;
			}
		}
		in_queue++;
		for (unsigned int i = 0; i < args -> modulus && particle_itr != particles -> end(); i++)
		{
			particle_itr ++;
		}
		if (particle_itr == particles -> end() || in_queue == block_size)
		{
			run_compute(temp_cparticles, args -> par_addr, args -> cache_addr, stream, results, args -> res_addr, in_queue);
			//while(!(args -> lock -> try_lock()));
			//std::cout << "Ran compute on:" << std::endl;
			for (uint16_t i = 0; i < in_queue; i++)
			{
				//std::cout << "(" << temp_cparticles[i].pos.x << ", " << temp_cparticles[i].pos.y << ", " << temp_cparticles[i].pos.z << "), size " << temp_cparticles[i].size << std::endl;
				//for (unsigned int j = 0; j < temp_cparticles[i].size; j++)
				//{
				//	std::cout << j << ": " << temp_cparticles[i].dependants[j] << std::endl;
				//}
				temp_vec = vector(results[i].x, results[i].y, results[i].z);
				par[i] -> set_acc_offset(&temp_vec);
			}
			//args -> lock -> unlock();
			in_queue = 0;
		}
		if (particle_itr == particles -> end()) { break; }
	}
	delete[] par;
	delete[] results;
	delete[] temp_cparticles;
	//while (!(args -> lock -> try_lock()));
	//std::cout << "Thread " << args -> thread_id << " exiting." << std::endl;
	//args -> lock -> unlock();
	pthread_exit(NULL);
}

void barnes_hut_cuda(particle_set *particles, octree *root)
{
	std::mutex master_lock;
	active_map active;
	stale_map stale;
	pthread_t *threads = new pthread_t[compute_threads];
	struct cuda_thread_data *td = new cuda_thread_data[compute_threads];
	cudaStream_t *streams = new cudaStream_t[compute_threads];
	cnode *cache_addr = init_cache();
	init_streams(streams);
	//std::cout << "GPU data address is " << cache_addr << std::endl;
	uint32_t last_free = 0;
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		td[i].stream = &(streams[i]);
		td[i].lock = &master_lock;
		td[i].particles = particles;
		td[i].root = root;
		td[i].active = &active;
		td[i].stale = &stale;
		td[i].thread_id = i;
		td[i].modulus = compute_threads;
		td[i].last_free = &last_free;
		td[i].cache_addr = cache_addr;
		td[i].par_addr = allocate_particles();
		td[i].res_addr = allocate_results();
		create_thread(&threads[i], NULL, barnes_hut_cuda_thread, (void*) &td[i]);
	}
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		pthread_join(threads[i], NULL);
		free_particles(td[i].par_addr);
		free_results(td[i].res_addr);
	}
	delete[] threads;
	delete[] td;
	delete[] streams;
	free_streams(streams);
	free_cache(cache_addr);
	call_dev_reset();
}
