#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <chrono>
#include <iomanip>

// Simple constant integrand 
template <typename Real>
Real constant_func(const Real* /*x*/, std::size_t) {
    return Real(1);
}

template <std::size_t Dim>
bool test_dimension() {
    using Real = double;
    using namespace std::chrono;
    
    // Define integration domain [0, 1]^Dim
    boost::math::cubature::hypercube<Real> box(Dim);
    for (std::size_t i = 0; i < Dim; ++i) {
        box.lower[i] = 0.0;
        box.upper[i] = 1.0;
    }
    
    // Use loose tolerances for high dimensions
    Real abs_tol = (Dim <= 5) ? 1e-6 : 1e-2;
    Real rel_tol = (Dim <= 5) ? 1e-6 : 1e-2;
    
    boost::math::cubature::execution_options options;
    options.max_eval = (Dim <= 5) ? 100000 : 50000;
    
    auto start = high_resolution_clock::now();
    
    auto result = boost::math::cubature::integrate_adaptive(
        constant_func<Real>, box, abs_tol, rel_tol, options);
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    
    Real expected = 1.0;  // Volume of unit hypercube
    Real error = std::abs(result.value - expected);
    bool passed = error < std::max(abs_tol, rel_tol);
    
    std::cout << std::setw(3) << Dim 
              << " | " << std::setw(10) << std::fixed << std::setprecision(6) << result.value
              << " | " << std::setw(10) << std::scientific << std::setprecision(2) << result.error
              << " | " << std::setw(8) << result.evaluations
              << " | " << std::setw(6) << duration << " ms"
              << " | " << (passed ? "PASS" : "FAIL") << "\n";
    
    return passed;
}

int main() {
    std::cout << "Testing adaptive integration for dimensions 2-15:\n";
    std::cout << "Integrating f(x) = 1 over [0,1]^d (expected result = 1.0)\n\n";
    std::cout << "Dim |   Result   |   Error    |  Evals   |  Time  | Status\n";
    std::cout << "----+------------+------------+----------+--------+-------\n";
    
    bool all_passed = true;
    all_passed &= test_dimension<2>();
    all_passed &= test_dimension<3>();
    all_passed &= test_dimension<4>();
    all_passed &= test_dimension<5>();
    all_passed &= test_dimension<6>();
    all_passed &= test_dimension<7>();
    all_passed &= test_dimension<8>();
    all_passed &= test_dimension<9>();
    all_passed &= test_dimension<10>();
    all_passed &= test_dimension<11>();
    all_passed &= test_dimension<12>();
    all_passed &= test_dimension<13>();
    all_passed &= test_dimension<14>();
    all_passed &= test_dimension<15>();
    
    std::cout << "\n" << (all_passed ? "All tests PASSED!" : "Some tests FAILED!") << "\n";
    return all_passed ? 0 : 1;
}