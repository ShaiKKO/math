#include <boost/math/cubature/policies.hpp>
#include <iostream>

using namespace boost::math::cubature;

int main() {
    std::cout << "Testing policy accumulator..." << std::endl;
    
    policy_accumulator<double> acc;
    
    // Test Kahan summation
    for (int i = 0; i < 1000000; ++i) {
        acc.add(0.1);
    }
    
    double result = acc.sum();
    double expected = 100000.0;
    double error = std::abs(result - expected);
    
    std::cout << "Sum of 1M × 0.1 = " << result << std::endl;
    std::cout << "Expected: " << expected << std::endl;
    std::cout << "Error: " << error << std::endl;
    
    if (error < 1e-10) {
        std::cout << "✓ Kahan summation working correctly" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Excessive numerical error" << std::endl;
        return 1;
    }
}
