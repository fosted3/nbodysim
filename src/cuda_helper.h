#ifndef cuda_helper_h_
#define cuda_helper_h_

#include <unordered_set>
#include "particle.h"
#include "octree.h"

void barnes_hut_cuda(std::unordered_set<particle*>*, octree*, bool, datatype);

#endif
