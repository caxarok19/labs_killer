#include "funcsh.h"
#include <fstream>

double func(double x) { return std::sin(x); }

int main()
{
	const double start = 0;
	const double end = 10;

	std::ofstream data2("data2_RR.txt");
	data2.precision(16);
	 
	std::ofstream data5("data5_RR.txt");
	data5.precision(16);

	std::ofstream data2_1("data2.txt");
	data2_1.precision(16);

	std::ofstream data5_1("data5.txt");
	data5_1.precision(16);

	const double I = 1 - std::cos(10);
	for (double h = -7; h < 1; h += 0.001) {
		double step = std::exp(h);

		data2_1 << h << "		" << std::abs(I - integrate<decltype(func), double, 2>
			(func, start, end, nodes<double, 2>::p, nodes<double, 2>::w, step)) << std::endl;

		data5_1 << h << "		" << std::abs(I - integrate<decltype(func), double, 5>
			(func, start, end, nodes<double, 5>::p, nodes<double, 5>::w, step)) << std::endl;

		 
	}

	for (double err = -25; err < -1; err += 0.001) {
		data2 << err << "		" << std::abs(I - integrateRR<decltype(func), double, 2>
			(func, start, end, nodes<double, 2>::p, nodes<double, 2>::w, std::exp(err))) << std::endl;

		data5 << err << "		" << std::abs(I - integrateRR<decltype(func), double, 5>
			(func, start, end, nodes<double, 5>::p, nodes<double, 5>::w, std::exp(err))) << std::endl;

	 
} 

data2.close();
data5.close();

return 0;
}