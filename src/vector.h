#ifndef vector_h_
#define vector_h_

class vector
{
	public:
		vector();
		vector(double, double, double);
		~vector();
		void print();
		void print_inline();
		vector& operator += (const vector&);
		vector& operator -= (const vector&);
		vector& operator *= (const double&);
		vector& operator /= (const double&);
		vector& operator = (const vector&);
		//const double operator[] (const int); //this doesn't seem to work properly (fix), temporarily replaced with get_x...
		const double get_x();
		const double get_y();
		const double get_z();
	private:
		double x;
		double y;
		double z;
};

#endif
