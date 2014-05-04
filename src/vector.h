#ifndef vector_h_
#define vector_h_

class vector
{
	public:
		vector();
		vector(double, double, double);
		~vector();
		void print(void);
		void print_inline(void);
		vector& operator += (const vector&);
		vector& operator -= (const vector&);
		vector& operator *= (const double&);
		vector& operator /= (const double&);
		vector& operator = (const vector&);
		//const double operator[] (const int); //this doesn't seem to work properly (fix), temporarily replaced with get_x, etc.
		double get_x(void);
		double get_y(void);
		double get_z(void);
		double magnitude(void);
		void normalize(void);
		void scale(double, double, double);
	private:
		double x;
		double y;
		double z;
};

double distance(vector*, vector*);
vector cross(vector&, vector&);

#endif
