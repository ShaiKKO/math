#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>
#include <vector>

using namespace boost::math::cubature;

int main() {
    // Test 1: Integrate constant function f(x,y) = 1 over [0,1]x[0,1]
    // Should give exactly 1.0
    std::cout << "Test 1: Constant function f(x,y) = 1 over unit square\n";
    
    auto constant = [](const std::vector<double>& x) -> double {
        return 1.0;
    };
    
    hypercube<double> unit_box(2);
    unit_box.lower = {0, 0};
    unit_box.upper = {1, 1};
    
    auto result1 = integrate_adaptive<double>(constant, unit_box, 1e-10, 1e-10);
    std::cout << "Result: " << result1.value << " (should be 1.0)\n";
    std::cout << "Error: " << result1.error << "\n";
    std::cout << "Evaluations: " << result1.evaluations << "\n\n";
    
    // Test 2: Integrate f(x,y) = x over [0,1]x[0,1]
    // Should give 0.5
    std::cout << "Test 2: Linear function f(x,y) = x over unit square\n";
    
    auto linear = [](const std::vector<double>& x) -> double {
        return x[0];
    };
    
    auto result2 = integrate_adaptive<double>(linear, unit_box, 1e-10, 1e-10);
    std::cout << "Result: " << result2.value << " (should be 0.5)\n";
    std::cout << "Error: " << result2.error << "\n";
    std::cout << "Evaluations: " << result2.evaluations << "\n\n";
    
    // Test 3: Integrate f(x,y) = x*y over [0,1]x[0,1]
    // Should give 0.25
    std::cout << "Test 3: Bilinear function f(x,y) = x*y over unit square\n";
    
    auto bilinear = [](const std::vector<double>& x) -> double {
        return x[0] * x[1];
    };
    
    auto result3 = integrate_adaptive<double>(bilinear, unit_box, 1e-10, 1e-10);
    std::cout << "Result: " << result3.value << " (should be 0.25)\n";
    std::cout << "Error: " << result3.error << "\n";
    std::cout << "Evaluations: " << result3.evaluations << "\n\n";
    
    // Test 4: The Gaussian from the main test
    std::cout << "Test 4: 2D Gaussian\n";
    auto gaussian = [](const std::vector<double>& x) -> double {
        double sigma = 0.5;
        double dx = x[0] - 0.5;
        double dy = x[1] - 0.5;
        return std::exp(-(dx*dx + dy*dy) / (2*sigma*sigma));
    };
    
    auto result4 = integrate_adaptive<double>(gaussian, unit_box, 1e-8, 1e-6);
    
    double sigma = 0.5;
    double analytical = 2 * M_PI * sigma * sigma * 
                       std::pow(std::erf(0.5 / (std::sqrt(2.0) * sigma)), 2);
    
    std::cout << "Result: " << result4.value << " (analytical: " << analytical << ")\n";
    std::cout << "Error: " << result4.error << "\n";
    std::cout << "Evaluations: " << result4.evaluations << "\n";
    
    return 0;
}