#ifndef data_structures_h_
#define data_structures_h_

#include "particle.h"
#include <unordered_set>
#include <unordered_map>

typedef std::pair<particle*, particle*> particle_pair;
typedef std::unordered_set<particle*> particle_set;
typedef std::unordered_set<std::pair<particle*, particle*> > particle_pair_set;
typedef std::unordered_map<particle*, particle*> particle_map;

particle_set::iterator& operator+=(particle_set::iterator&, const unsigned int &);

#endif
