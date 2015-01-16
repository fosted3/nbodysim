#include "particle.h"
#include "octree.h"
#include "thread_functions.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "vector_types.h"
#include <cassert>
//#include "cuda.h"
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

/*struct cache_data
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
};*/

typedef std::unordered_map<octree*, uint32_t> active_map;
typedef std::unordered_map<octree*, uint32_t> stale_map;

struct cuda_thread_data
{
	//cudaStream_t *stream;
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
	particle **par;
	uint16_t *current_pos;
	datatype3 *host_res;
	cparticle *host_par;
	uint16_t *last_man_standing;
};

void *barnes_hut_cuda_thread(void *data)
{
	struct cuda_thread_data *args;
	args = (struct cuda_thread_data*) data;
	//cudaStream_t *stream = args -> stream;
	particle_set *particles = args -> particles;
	octree *root = args -> root;
	active_map *active = args -> active;
	stale_map *stale = args -> stale;
	datatype theta = 1.0f;
	octree *node;
	std::queue<octree*> nodes;
	active_map::iterator active_itr;
	stale_map::iterator stale_itr;
	particle_set::iterator particle_itr = particles -> begin();
	for (unsigned int i = 0; i < args -> thread_id && particle_itr != particles -> end(); i++)
	{
		particle_itr ++;
	}
	uint32_t gpu_index;
	uint16_t dependancy_index;
	cnode temp_cnode;
	vector temp_vec;
	datatype percent;
	uint64_t completed = 0;
	cparticle temp_cparticle;
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
		/*par[in_queue] = *particle_itr;
		par[in_queue] -> set_acc_zero();
		temp_cparticles[in_queue].pos.x = par[in_queue] -> get_pos() -> get_x();
		temp_cparticles[in_queue].pos.y = par[in_queue] -> get_pos() -> get_y();
		temp_cparticles[in_queue].pos.z = par[in_queue] -> get_pos() -> get_z();
		temp_cparticles[in_queue].size = 0;*/
		temp_cparticle.pos.x = 0;
		temp_cparticle.pos.y = 0;
		temp_cparticle.pos.z = 0;
		temp_cparticle.size = 0;
		dependancy_index = 0;
		nodes.push(root);
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
			else if ((*particle_itr) != node -> get_particle()) //node should be in cache
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
						assert(stale -> size() > 0);
						gpu_index = stale -> begin() -> second;
						stale -> erase(stale -> begin());
					}
					else
					{
						//std::cout << "Cache has space, inserting at " << *(args -> last_free) << std::endl;
						gpu_index = *(args -> last_free);
						(*(args -> last_free))++;
					}
					copy_to_gpu(args -> cache_addr, gpu_index, &temp_cnode); //copy to gpu returning index for node
					active -> insert(std::make_pair(node, gpu_index));
				}
				else if (active_itr != active -> end()) //Being used by a node
				{
					//std::cout << node << " in cache already." << std::endl;
					//active_itr -> second -> in_use[(args -> thread_id)] = true;
					//gpu_index = active_itr -> second -> pos;
					gpu_index = active_itr -> second;
				}
				else if (stale_itr != stale -> end()) //Not being used by a node, move to active
				{
					//std::cout << node << " is stale, moving to active" << std::endl;
					gpu_index = stale_itr -> second;
					active -> insert(std::make_pair(node, gpu_index));
					stale -> erase(stale_itr);
				}
				args -> lock -> unlock();
				assert(dependancy_index < shared_size);
				temp_cparticle.dependants[dependancy_index] = gpu_index;
				dependancy_index ++;
				temp_cparticle.size = dependancy_index;
			}
		}
		//completed calculating deps for *particle_itr
		///in_queue++;
		while (!(args -> lock -> try_lock()));
		(args -> par)[*(args -> current_pos)] = *particle_itr;
		(args -> host_par)[*(args -> current_pos)].pos.x = temp_cparticle.pos.x;
		(args -> host_par)[*(args -> current_pos)].pos.y = temp_cparticle.pos.y;
		(args -> host_par)[*(args -> current_pos)].pos.z = temp_cparticle.pos.z;
		(args -> host_par)[*(args -> current_pos)].size = temp_cparticle.size;
		for (uint16_t i = 0; i < temp_cparticle.size; i++)
		{
			(args -> host_par)[*(args -> current_pos)].dependants[i] = temp_cparticle.dependants[i];
		}		
		(*(args -> current_pos))++;
		if (*(args -> current_pos) == block_size)
		{
			run_compute(args -> host_par, args -> par_addr, args -> cache_addr, args -> host_res, args -> res_addr, block_size);
			for (uint16_t i = 0; i < block_size; i++)
			{
				temp_vec = vector((args -> host_res)[i].x, (args -> host_res)[i].y, (args -> host_res)[i].z);
				(args -> par)[i] -> set_acc_offset(&temp_vec);
			}
			(*(args -> current_pos)) = 0;
			/*active_itr = active -> begin();
			while (active_itr != active -> end())
			{
				stale_itr = stale -> find(active_itr -> first);
				assert(stale_itr == stale -> end());
				stale -> insert(std::make_pair(active_itr -> first, active_itr -> second));
				active -> erase(active_itr);
			}*/
		}
		args -> lock -> unlock();
		for (unsigned int i = 0; i < args -> modulus && particle_itr != particles -> end(); i++)
		{
			particle_itr ++;
		}
		/*
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
		}*/
		if (particle_itr == particles -> end())
		{
			while(!(args -> lock -> try_lock()));
			(*(args -> last_man_standing))--;
			if (*(args -> last_man_standing) == 0 && particles -> size() % block_size != 0)
			{
				run_compute(args -> host_par, args -> par_addr, args -> cache_addr, args -> host_res, args -> res_addr, *(args -> current_pos));
				for (uint16_t i = 0; i < *(args -> current_pos); i++)
				{
					temp_vec = vector((args -> host_res)[i].x, (args -> host_res)[i].y, (args -> host_res)[i].z);
					(args -> par)[i] -> set_acc_offset(&temp_vec);
				}
			}			
			args -> lock -> unlock();
			break;
		}
	}
	while (!(args -> lock -> try_lock()));
	//std::cout << "Thread " << args -> thread_id << " exiting." << std::endl;
	args -> lock -> unlock();
	pthread_exit(NULL);
}

void barnes_hut_cuda(particle_set *particles, octree *root)
{
	std::mutex master_lock;
	active_map active;
	stale_map stale;
	pthread_t *threads = new pthread_t[compute_threads];
	struct cuda_thread_data *td = new cuda_thread_data[compute_threads];
	particle **par = new particle*[block_size];
	cparticle *host_par = new cparticle[block_size];
	datatype3 *host_res = new datatype3[block_size];
	//cudaStream_t *streams = new cudaStream_t[compute_threads];
	cnode *cache_addr = init_cache();
	//init_streams(streams);
	//std::cout << "GPU data address is " << cache_addr << std::endl;
	uint32_t last_free = 0;
	uint16_t current_pos = 0;
	cparticle *par_addr = allocate_particles();
	datatype3 *res_addr = allocate_results();
	uint16_t last_man_standing = compute_threads;
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		//td[i].stream = &(streams[i]);
		td[i].lock = &master_lock;
		td[i].particles = particles;
		td[i].root = root;
		td[i].active = &active;
		td[i].stale = &stale;
		td[i].thread_id = i;
		td[i].modulus = compute_threads;
		td[i].last_free = &last_free;
		td[i].cache_addr = cache_addr;
		td[i].par_addr = par_addr;
		td[i].res_addr = res_addr;
		td[i].par = par;
		td[i].current_pos = &current_pos;
		td[i].host_par = host_par;
		td[i].host_res = host_res;
		td[i].last_man_standing = &last_man_standing;
		create_thread(&threads[i], NULL, barnes_hut_cuda_thread, (void*) &td[i]);
		
	}
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	free_particles(par_addr);
	free_results(res_addr);
	delete[] threads;
	delete[] td;
	delete[] par;
	delete[] host_par;
	delete[] host_res;
	//free_streams(streams);
	//delete[] streams;
	free_cache(cache_addr);
	call_dev_reset();
}
