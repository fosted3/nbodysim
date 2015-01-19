#ifndef cuda_code_h_
#define cuda_code_h_

#include "cuda_types.h"

void init_streams(cudaStream_t*);
void free_streams(cudaStream_t*);
void run_compute(cparticle*, cparticle*, cnode*, cnode*, datatype3*, datatype3*, uint32_t, uint32_t, uint16_t, bool, cudaStream_t*);
cparticle* allocate_dev_particles(void);
cparticle* allocate_host_particles(void);
cnode* allocate_dev_nodes(void);
cnode* allocate_host_nodes(void);
datatype3* allocate_dev_results(void);
datatype3* allocate_host_results(void);
void free_dev_particles(cparticle*);
void free_host_particles(cparticle*);
void free_dev_nodes(cnode*);
void free_host_nodes(cnode*);
void free_dev_results(datatype3*);
void free_host_results(datatype3*);
void call_dev_reset(void);

#endif
