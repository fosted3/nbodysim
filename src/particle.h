#ifndef particle_h_
#define particle_h_

#include "vector.h"

class particle
{
	public:
		particle();
		particle(vector*, vector*, vector*, double);
		particle(particle&);
		double get_mass(void);
		vector* get_pos(void);
		void set_acc_zero(void);
		void set_acc_offset(vector*);
		void update(double);
		void update(double&, double&);
		void print(void);
	private:
		vector position;
		vector velocity;
		vector acceleration;
		double mass;
		//double radius;
};

#endif
