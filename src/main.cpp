#include <iostream>
#include "quadtree.h"
#include "particle.h"
#include "vector.h"

int main()
{
	vector *a = new vector(1, 2, 3);
	a -> print();
	delete a;
	return 0;
}
