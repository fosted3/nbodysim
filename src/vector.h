#ifndef vector_h_
#define vector_h_

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

class vector
{
	public:
		vector();
		vector(datatype, datatype, datatype);
		~vector();
		void print(void);
		void print_inline(void);
		vector& operator += (const vector&);
		vector& operator -= (const vector&);
		vector& operator *= (const datatype&);
		vector& operator /= (const datatype&);
		vector& operator = (const vector&);
		//const double operator[] (const int); //this doesn't seem to work properly (fix), temporarily replaced with get_x, etc.
		datatype get_x(void);
		datatype get_y(void);
		datatype get_z(void);
		datatype magnitude(void);
		void normalize(void);
		void scale(datatype, datatype, datatype);
	private:
		datatype x;
		datatype y;
		datatype z;
};

datatype distance(vector*, vector*);
vector cross(vector&, vector&);

#endif
