#include <boost/math/cubature/transforms.hpp>
#include <iostream>
#include <iomanip>

using boost::math::cubature::duffy_transform;

int main() {
    double x[2] = {0.3, 0.2};
    double u[2];
    
    duffy_transform<double>::inverse(x, u, 2);
    
    std::cout << "Original x: [" << x[0] << ", " << x[1] << "]\n";
    std::cout << "Inverse u: [" << u[0] << ", " << u[1] << "]\n";
    
    double x_back[2];
    double jac = duffy_transform<double>::apply(u, x_back, 2);
    
    std::cout << "Forward x_back: [" << x_back[0] << ", " << x_back[1] << "]\n";
    std::cout << "Jacobian: " << jac << "\n";
    
    // Manual calculation
    // Forward: x[0] = u[0] * (1 - u[1]), x[1] = u[0] * u[1]
    // So with u[0] = 0.3, u[1] = 0.2/0.3 = 0.666...
    double u1_calc = 0.2 / 0.3;
    std::cout << "\nManual calc:\n";
    std::cout << "u[1] = " << u1_calc << "\n";
    std::cout << "x[0] should be: " << u[0] * (1 - u1_calc) << "\n";
    std::cout << "x[1] should be: " << u[0] * u1_calc << "\n";
    
    return 0;
}
