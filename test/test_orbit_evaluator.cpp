#include <boost/math/cubature/detail/orbit_evaluator.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <boost/math/cubature/regions.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

// Test function: f(x,y) = exp(-x^2 - y^2) with sharp peak at origin
template <typename Real>
Real gaussian_2d(const Real* x) {
    return std::exp(-x[0]*x[0] - x[1]*x[1]);
}

// Test function with strong variation along x-axis
template <typename Real>
Real asymmetric_2d(const Real* x) {
    return std::exp(-10*x[0]*x[0] - x[1]*x[1]);
}

template <typename Real>
void test_orbit_evaluation() {
    std::cout << "Testing orbit-based evaluation..." << std::endl;
    
    // Create a region
    region<Real> reg(2);
    reg.a = {-1, -1};
    reg.b = {1, 1};
    
    // Test with Gaussian
    typename orbit_evaluator<Real>::template orbit_values<2> values;
    orbit_evaluator<Real>::template evaluate_orbits<decltype(gaussian_2d<Real>), 2, 9>(
        gaussian_2d<Real>, reg, values);
    
    std::cout << "  Center value: " << values.f_center << std::endl;
    std::cout << "  Total evaluations: " << values.total_evaluations << std::endl;
    
    // Compute fourth differences
    auto diffs = orbit_evaluator<Real>::compute_fourth_differences(values);
    std::cout << "  Fourth differences: [" << diffs[0] << ", " << diffs[1] << "]" << std::endl;
    
    // For symmetric Gaussian, differences should be nearly equal
    Real ratio = diffs[0] / diffs[1];
    std::cout << "  Ratio of differences (should be ~1 for symmetric): " << ratio << std::endl;
    assert(std::abs(ratio - 1.0) < 0.1);
    
    // The orbit evaluator is now used for computing fourth differences
    // Actual integration is done in genz_malik_evaluator using proper weights
    std::cout << "  Orbit evaluation successful" << std::endl;
    std::cout << "  Fourth differences computed for axis selection" << std::endl;
    
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_asymmetric_splitting() {
    std::cout << "\nTesting asymmetric function splitting..." << std::endl;
    
    // Create a region
    region<Real> reg(2);
    reg.a = {-1, -1};
    reg.b = {1, 1};
    
    // Test with asymmetric function
    typename orbit_evaluator<Real>::template orbit_values<2> values;
    orbit_evaluator<Real>::template evaluate_orbits<decltype(asymmetric_2d<Real>), 2, 9>(
        asymmetric_2d<Real>, reg, values);
    
    // Compute fourth differences
    auto diffs = orbit_evaluator<Real>::compute_fourth_differences(values);
    std::cout << "  Fourth differences: [" << diffs[0] << ", " << diffs[1] << "]" << std::endl;
    
    // For asymmetric function, x-direction should have larger difference
    Real ratio = diffs[0] / diffs[1];
    std::cout << "  Ratio of differences (should be >1 for x-asymmetric): " << ratio << std::endl;
    assert(ratio > 1.5);  // x-direction should have significantly larger variation
    
    // Select split dimension
    std::size_t split_dim = orbit_evaluator<Real>::select_split_dimension(diffs, reg);
    std::cout << "  Selected split dimension: " << split_dim << std::endl;
    assert(split_dim == 0);  // Should split along x-axis
    
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_genz_malik_enhanced() {
    std::cout << "\nTesting enhanced Genz-Malik evaluator..." << std::endl;
    
    // Create a region
    region<Real> reg(2);
    reg.a = {0, 0};
    reg.b = {1, 1};
    
    // Lambda wrapper for Gaussian
    auto f = [](const std::vector<Real>& x) -> Real {
        Real dx = x[0] - 0.5;
        Real dy = x[1] - 0.5;
        return std::exp(-10*(dx*dx + dy*dy));
    };
    
    // Evaluate using enhanced evaluator
    embedded_pair_result<Real> result;
    bool success = genz_malik_evaluator<Real>::template evaluate_embedded_pair<decltype(f), 2>(
        f, reg, result);
    
    assert(success);
    
    std::cout << "  Fine estimate: " << result.estimate_fine << std::endl;
    std::cout << "  Coarse estimate: " << result.estimate_coarse << std::endl;
    std::cout << "  Error estimate: " << result.embedded_error << std::endl;
    std::cout << "  Spread estimate: " << result.spread_estimate << std::endl;
    std::cout << "  Split dimension: " << result.split_dimension << std::endl;
    std::cout << "  Fourth diffs: [" << result.fourth_differences[0] 
              << ", " << result.fourth_differences[1] << "]" << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    // Check that we have cached values
    assert(result.cached_values != nullptr);
    
    std::cout << "  PASSED" << std::endl;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n=== Testing Orbit-Based Genz-Malik Evaluation ===" << std::endl;
    
    test_orbit_evaluation<double>();
    test_asymmetric_splitting<double>();
    test_genz_malik_enhanced<double>();
    
    // Test with higher precision if available
    if (sizeof(long double) > sizeof(double)) {
        std::cout << "\n=== Testing with long double ===" << std::endl;
        test_orbit_evaluation<long double>();
        test_asymmetric_splitting<long double>();
    }
    
    std::cout << "\nAll tests passed!" << std::endl;
    return 0;
}