#ifndef cuda_h_
#define cuda_h_

#include "cuda_types.h"

void copy_to_gpu(cnode*, uint32_t, cudaStream_t *, cnode*);
cnode* init_cache(void);
void free_cache(cnode*);
void init_streams(cudaStream_t*);
void free_streams(cudaStream_t*);
void run_compute(cparticle*, cparticle*, cnode*, cudaStream_t*, datatype3*, datatype3*, uint16_t);
cparticle* allocate_particles(void);
datatype3* allocate_results(void);
void free_particles(cparticle*);
void free_results(datatype3*);

#endif
