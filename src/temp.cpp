#include "data_structures.h"

particle_set::iterator& operator+=(particle_set::iterator &itr, const unsigned int &inc)
{
	for (unsigned int x = 0; x < inc; x++)
	{
		*itr ++;
	}
	return itr;
}
