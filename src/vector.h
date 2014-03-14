#ifndef vector_h_
#define vector_h_

class vector
{
	public:
		vector();
		vector(double, double, double);
		~vector();
		void print();
	private:
		double x;
		double y;
		double z;
};

vector operator = (const vector&);
vector operator + (const vector&, const vector&);
vector operator * (const vector&, const vector&);
vector operator - (const vector&, const vector&);
vector operator += (const vector&);
vector operator -= (const vector&);

#endif
