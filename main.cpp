#include "NACA_Airfoils.hpp"
#include <iomanip>
using namespace std;

int main()
{
	naca::airfoil<float> wing{ "2412", 20, naca::airfoil_spacing::cosine,  naca::airfoil_trailing_edge::open };

	for(auto& point: wing.get_points()){
		std::cout << std::left << std::setw(15) << point.x << point.y << '\n';
	}
	
	return 0;
}
