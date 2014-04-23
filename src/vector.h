#ifndef vector_h_
#define vector_h_

class vector
{
	public:
		vector();
		vector(double, double, double);
		~vector();
		void print();
		vector& operator += (const vector&);
		vector& operator -= (const vector&);
		vector& operator = (const vector&);		
	private:
		double x;
		double y;
		double z;
};

#endif
