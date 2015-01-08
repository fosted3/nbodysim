#include "cuda.h"
#include <vector>
#include "vector_types.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
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

__global__ void compute(cparticle *particles, unsigned int p_size, cnode *nodes, unsigned int n_size, datatype3 *results)
{
	uint16_t tid = threadIdx.x;
	uint32_t bid = blockIdx.x;
	__shared__ datatype3 acc[shared_size];
	__shared__ cnode dep[shared_size];
	__shared__ cparticle par;
	__shared__ datatype r_sq[shared_size];
	if (tid == 0)
	{
		par = particles[bid];
	}
	__syncthreads();
	if (tid < par.size)
	{
		dep[tid] = nodes[particles[bid].dependants[tid]];
		acc[tid].x = dep[tid].pos.x - par.pos.x;
		acc[tid].y = dep[tid].pos.y - par.pos.y;
		acc[tid].z = dep[tid].pos.z - par.pos.z;
		r_sq[tid] = pow(length(acc[tid]), -2.0f);
		r_sq[tid] *= 6.67384e-11f * dep[tid].mass / length(acc[tid]);
		mul_datatype3(acc[tid], r_sq[tid]);
	}
	__syncthreads();	
	for (unsigned int s = 1; s < shared_size; s *= 2)
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

void run_compute(std::vector<cparticle> *host_particles_vector, std::vector<cnode> *host_nodes_vector, std::vector<datatype3> *host_results_vector)
{
	thrust::device_vector<cparticle> device_particles(host_particles_vector -> begin(), host_particles_vector -> end());
	thrust::device_vector<cnode> device_nodes(host_nodes_vector -> begin(), host_nodes_vector -> end());
	thrust::device_vector<datatype3> device_results(host_results_vector -> begin(), host_results_vector -> end());
	dim3 grid(device_particles.size(),1,1);
	dim3 block(shared_size,1,1);
	cparticle* device_particles_ptr = thrust::raw_pointer_cast(&device_particles[0]);
	cnode* device_nodes_ptr = thrust::raw_pointer_cast(&device_nodes[0]);
	datatype3* device_results_ptr = thrust::raw_pointer_cast(&device_results[0]);
	assert(device_particles.size() == device_results.size());
	compute<<<grid, block>>>(device_particles_ptr, device_particles.size(), device_nodes_ptr, device_nodes.size(), device_results_ptr);
	assert(cudaThreadSynchronize() == cudaSuccess);
	host_results_vector -> resize(device_results.size());
	thrust::copy(device_results.begin(), device_results.end(), host_results_vector -> begin());
	device_particles.clear();
	device_nodes.clear();
	device_results.clear();	
}
