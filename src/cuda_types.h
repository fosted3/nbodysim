#ifndef cuda_types_h_
#define cuda_types_h_

#include <stdint.h>
#include "cuda.h"


#define compute_threads 8
#define block_size 1024
#define shared_size 64

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

struct cparticle
{
	datatype3 pos;
	uint32_t dependants[shared_size];
	uint16_t size;
};

struct cnode
{
	datatype3 pos;
	datatype mass;
};

#endif
