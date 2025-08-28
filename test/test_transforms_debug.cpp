#include <iostream>
#include <cmath>
#include <boost/math/cubature/transforms.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::cubature;
using namespace boost::math::constants;

int main() {
    using Real = double;
    
    std::cout << "Testing transforms..." << std::endl;
    
    // Test tangent transform
    tangent_transform<Real> t_tan;
    
    std::cout << "\nTangent transform [0,1] -> (-inf, inf):" << std::endl;
    Real test_points[] = {0.1, 0.25, 0.5, 0.75, 0.9};
    
    for (Real u : test_points) {
        auto [x, jac] = t_tan.forward(u);
        std::cout << "  u=" << u << " -> x=" << x << ", jacobian=" << jac << std::endl;
    }
    
    // Test rational transform
    rational_transform<Real> t_rat;
    
    std::cout << "\nRational transform [0,1] -> [0, inf):" << std::endl;
    for (Real u : test_points) {
        auto [x, jac] = t_rat.forward(u);
        std::cout << "  u=" << u << " -> x=" << x << ", jacobian=" << jac << std::endl;
    }
    
    // Test simple integral: exp(-x^2) from -inf to inf
    std::cout << "\nManual integration test of exp(-x^2) using tangent transform:" << std::endl;
    
    // Gauss-Legendre quadrature points for [0,1]
    const int n_points = 5;
    Real gl_nodes[] = {0.046910077, 0.230765345, 0.5, 0.769234655, 0.953089923};
    Real gl_weights[] = {0.118463443, 0.239314335, 0.284444444, 0.239314335, 0.118463443};
    
    Real sum = 0;
    for (int i = 0; i < n_points; ++i) {
        auto [x, jac] = t_tan.forward(gl_nodes[i]);
        Real integrand = std::exp(-x*x) * jac;
        sum += gl_weights[i] * integrand;
        std::cout << "  Node " << i << ": u=" << gl_nodes[i] 
                  << ", x=" << x << ", f(x)=" << std::exp(-x*x) 
                  << ", jacobian=" << jac << ", contribution=" << gl_weights[i] * integrand << std::endl;
    }
    
    std::cout << "  Numerical result: " << sum << std::endl;
    std::cout << "  Exact result: " << std::sqrt(pi<Real>()) << " = " << std::sqrt(pi<Real>()) << std::endl;
    std::cout << "  Error: " << std::abs(sum - std::sqrt(pi<Real>())) << std::endl;
    
    return 0;
}