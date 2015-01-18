#ifndef cuda_code_h_
#define cuda_code_h_

#include "cuda_types.h"

void init_streams(cudaStream_t*);
void free_streams(cudaStream_t*);
void run_compute(cparticle*, cparticle*, cnode*, cnode*, datatype3*, datatype3*, uint32_t, uint32_t, uint16_t, bool, cudaStream_t*);
cparticle* allocate_particles(void);
cnode* allocate_nodes(void);
datatype3* allocate_results(void);
void free_particles(cparticle*);
void free_nodes(cnode*);
void free_results(datatype3*);
void call_dev_reset(void);

#endif
