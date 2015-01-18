#include "cuda.h"
#include <vector>
#include "vector_types.h"
#include <cassert>
#include "cuda_types.h"
#include <iostream>
#include "stdio.h"

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

#define handle_error(ans) { cuda_assert((ans), __FILE__, __LINE__); }
inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

__device__ datatype length(datatype3 v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

__device__ void mul_datatype3(datatype3 &a, const datatype &b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}

__global__ void compute(cparticle *particle, cnode *nodes, datatype3 *results)
{
	uint16_t tid = threadIdx.x;
	uint16_t bid = blockIdx.x;
	__shared__ datatype3 acc[shared_size];
	__shared__ cnode dep[shared_size];
	__shared__ cparticle par;
	__shared__ datatype r_sq[shared_size];
	if (tid == 0) { par = particle[bid]; }
	__syncthreads();
	if (tid < par.size)
	{
		dep[tid] = nodes[par.dependants[tid]];
		acc[tid].x = dep[tid].pos.x - par.pos.x;
		acc[tid].y = dep[tid].pos.y - par.pos.y;
		acc[tid].z = dep[tid].pos.z - par.pos.z;
		r_sq[tid] = pow(length(acc[tid]), -2.0f);
		r_sq[tid] *= 6.67384e-11f * dep[tid].mass / length(acc[tid]);
		mul_datatype3(acc[tid], r_sq[tid]);
	}
	__syncthreads();
	for (uint16_t s = 1; s < shared_size; s *= 2)
	{
		if (tid % (2 * s) == 0 && tid + s < par.size)
		{
			acc[tid].x = acc[tid].x + acc[tid + s].x;
			acc[tid].y = acc[tid].y + acc[tid + s].y;
			acc[tid].z = acc[tid].z + acc[tid + s].z;
		}
		__syncthreads();
	}
	results[bid] = acc[0];
}

void init_streams(cudaStream_t *streams)
{
	for(unsigned int i = 0; i < compute_threads; i++)
	{
		//std::cout << "Initializing stream " << i << std::endl;
		handle_error(cudaStreamCreate(&(streams[i])));
	}
}

void free_streams(cudaStream_t *streams)
{
	for (unsigned int i = 0; i < compute_threads; i++)
	{
		//std::cout << "Destroying stream " << i << std::endl;
		handle_error(cudaStreamDestroy(streams[i]));
	}
}

cparticle* allocate_particles(void)
{
	cparticle *addr = NULL;
	handle_error(cudaMalloc(&addr, sizeof(cparticle) * block_size));
	return addr;
}

cnode* allocate_nodes(void)
{
	cnode *addr = NULL;
	handle_error(cudaMalloc(&addr, sizeof(cnode) * block_size * shared_size));
	return addr;
}

datatype3* allocate_results(void)
{
	datatype3 *addr = NULL;
	handle_error(cudaMalloc(&addr, sizeof(datatype3) * block_size));
	return addr;
}

void free_particles(cparticle *addr)
{
	handle_error(cudaFree(addr));
}

void free_nodes(cnode *addr)
{
	handle_error(cudaFree(addr));
}

void free_results(datatype3 *addr)
{
	handle_error(cudaFree(addr));
}

void run_compute(cparticle *particles, cparticle *par_addr, cnode *node, cnode *node_addr, datatype3 *results, datatype3 *res_addr, uint32_t par_size, uint32_t node_size, uint16_t threads, cudaStream_t *stream)
{
	handle_error(cudaMemcpyAsync(par_addr, particles, sizeof(cparticle) * par_size, cudaMemcpyHostToDevice, *stream));
	handle_error(cudaMemcpyAsync(node_addr, node, sizeof(cnode) * node_size, cudaMemcpyHostToDevice, *stream));
	compute<<<par_size, threads, 0, *stream>>>(par_addr, node_addr, res_addr);
	handle_error(cudaMemcpyAsync(results, res_addr, sizeof(datatype3) * par_size, cudaMemcpyDeviceToHost, *stream));
	handle_error(cudaStreamSynchronize(*stream));
}

void call_dev_reset(void)
{
	cudaDeviceReset();
}
