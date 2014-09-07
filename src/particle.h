#ifndef particle_h_
#define particle_h_

#ifdef DOUBLE
#ifndef datatype
#define datatype double
#endif
#endif
#ifdef FLOAT
#ifndef datatype
#define datatype float
#endif
#endif

#include "vector.h"

class particle
{
	public:
		particle();
		particle(vector*, vector*, vector*, datatype);
		particle(particle&);
		datatype get_mass(void);
		vector* get_pos(void);
		void set_acc_zero(void);
		void set_acc_offset(vector*);
		void update(datatype);
		void update(datatype&, datatype&);
		void print(void);
	private:
		vector position;
		vector velocity;
		vector acceleration;
		datatype mass;
		//double radius;
};

#endif
