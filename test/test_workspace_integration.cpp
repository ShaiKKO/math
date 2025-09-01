#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/workspace.hpp>
#include <iostream>
#include <chrono>

// Simple test integrand
template <typename Real>
Real test_function(const Real* x, std::size_t) {
    return x[0] * x[0] + x[1] * x[1];  // f(x,y) = x^2 + y^2
}

int main() {
    using Real = double;
    using namespace std::chrono;
    
    // Create integration domain [0,1]^2
    boost::math::cubature::hypercube<Real> box(2);
    box.lower[0] = 0.0; box.upper[0] = 1.0;
    box.lower[1] = 0.0; box.upper[1] = 1.0;
    
    Real abs_tol = 1e-8;
    Real rel_tol = 1e-8;
    
    // Test 1: Without workspace
    std::cout << "Test 1: Integration without workspace\n";
    auto start = high_resolution_clock::now();
    
    auto result1 = boost::math::cubature::integrate_adaptive(
        test_function<Real>, box, abs_tol, rel_tol);
    
    auto end = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(end - start).count();
    
    std::cout << "  Result: " << result1.value << "\n";
    std::cout << "  Error: " << result1.error << "\n";
    std::cout << "  Evaluations: " << result1.evaluations << "\n";
    std::cout << "  Time: " << duration1 << " μs\n\n";
    
    // Test 2: With workspace (multiple integrations)
    std::cout << "Test 2: Integration with workspace (3 runs)\n";
    boost::math::cubature::workspace<Real> ws;
    
    Real total_time_with_ws = 0;
    for (int i = 0; i < 3; ++i) {
        start = high_resolution_clock::now();
        
        auto result2 = boost::math::cubature::integrate_adaptive(
            test_function<Real>, box, abs_tol, rel_tol, ws);
        
        end = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(end - start).count();
        total_time_with_ws += duration2;
        
        std::cout << "  Run " << (i+1) << ": " << result2.value 
                  << " (" << duration2 << " μs)\n";
    }
    
    std::cout << "  Average time with workspace: " 
              << total_time_with_ws / 3.0 << " μs\n";
    
    // Expected result: ∫∫(x^2 + y^2) dx dy from 0 to 1 = 2/3
    Real expected = 2.0 / 3.0;
    std::cout << "\nExpected value: " << expected << "\n";
    std::cout << "Error: " << std::abs(result1.value - expected) << "\n";
    
    if (std::abs(result1.value - expected) < 1e-6) {
        std::cout << "Test PASSED!\n";
        return 0;
    } else {
        std::cout << "Test FAILED!\n";
        return 1;
    }
}