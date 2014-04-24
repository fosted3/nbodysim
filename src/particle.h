#ifndef particle_h_
#define particle_h_

#include "vector.h"
#include "quadtree.h"

class quadtree;

class particle
{
	public:
		particle();
		particle(vector*, vector*, vector*, double);
		void set_container(quadtree*);
		double get_mass(void);
		vector* get_pos(void);
		void set_acc_zero(void);
		void set_acc_offset(vector*);
		void update(double dt);
		void print(void);
		quadtree* get_container(void);
	private:
		vector position;
		vector velocity;
		vector acceleration;
		double mass;
		quadtree* container;
		//double radius;
};

#endif
