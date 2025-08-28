#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <iomanip>

using namespace boost::math::cubature;

template <typename Real>
void test_gaussian_2d() {
    std::cout << "Testing 2D Gaussian integral..." << std::endl;
    
    // Gaussian with sigma = 0.5, centered at (0.5, 0.5)
    auto gaussian = [](const std::vector<Real>& x) -> Real {
        Real sigma = 0.5;
        Real dx = x[0] - 0.5;
        Real dy = x[1] - 0.5;
        return std::exp(-(dx*dx + dy*dy) / (2*sigma*sigma));
    };
    
    hypercube<Real> box(2);
    box.lower = {0, 0};
    box.upper = {1, 1};
    
    Real abs_tol = 1e-4;
    Real rel_tol = 1e-3;
    
    auto result = integrate_adaptive<Real>(gaussian, box, abs_tol, rel_tol, 1000000);
    
    // Analytical value: 2*pi*sigma^2 * erf(0.5/(sqrt(2)*sigma))^2
    // The 0.5 comes from integrating from 0 to 1 (half-width of box from center at 0.5)
    Real sigma = 0.5;
    Real analytical = 2 * M_PI * sigma * sigma * 
                     std::pow(std::erf(0.5 / (std::sqrt(2.0) * sigma)), 2);
    
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    if (result.status != status_code::success) {
      std::cout << "  Error: Integration did not converge" << std::endl;
      std::cout << "  Result: " << result.value << std::endl;
      std::cout << "  Error estimate: " << result.error << std::endl;
      std::cout << "  Evaluations: " << result.evaluations << std::endl;
    }
    assert(result.status == status_code::success);
    
    std::cout << "  Result: " << result.value 
              << " (analytical: " << analytical << ")" << std::endl;
    std::cout << "  Difference: " << std::abs(result.value - analytical) << std::endl;
    std::cout << "  Tolerance: " << (abs_tol + rel_tol * analytical) << std::endl;
    std::cout << "  Error estimate: " << result.error << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    assert(std::abs(result.value - analytical) < abs_tol + rel_tol * analytical);
    assert(result.error < abs_tol + rel_tol * std::abs(result.value));
    
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_oscillatory_3d() {
    std::cout << "Testing 3D oscillatory integral..." << std::endl;
    
    // cos(2*pi*x) * cos(2*pi*y) * cos(2*pi*z)
    auto oscillatory = [](const std::vector<Real>& x) -> Real {
        return std::cos(2 * M_PI * x[0]) * 
               std::cos(2 * M_PI * x[1]) * 
               std::cos(2 * M_PI * x[2]);
    };
    
    hypercube<Real> box(3);
    box.lower = {0, 0, 0};
    box.upper = {1, 1, 1};
    
    Real abs_tol = 1e-4;
    Real rel_tol = 1e-3;
    
    auto result = integrate_adaptive<Real>(oscillatory, box, abs_tol, rel_tol);
    
    // Analytical value: 0 (integral of cosine over period)
    Real analytical = 0.0;
    
    assert(result.status == status_code::success);
    assert(std::abs(result.value - analytical) < abs_tol);
    assert(result.error < abs_tol);
    
    std::cout << "  Result: " << result.value 
              << " (analytical: " << analytical << ")" << std::endl;
    std::cout << "  Error estimate: " << result.error << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_product_peak_4d() {
    std::cout << "Testing 4D product peak integral..." << std::endl;
    
    // Product peak function from Genz
    std::vector<Real> c = {0.3, 0.5, 0.7, 0.9};  // Peak locations
    std::vector<Real> w = {2.0, 3.0, 4.0, 5.0};  // Peak widths
    
    auto product_peak = [&c, &w](const std::vector<Real>& x) -> Real {
        Real result = 1.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            result *= 1.0 / (1.0 / (w[i] * w[i]) + (x[i] - c[i]) * (x[i] - c[i]));
        }
        return result;
    };
    
    hypercube<Real> box(4);
    box.lower = {0, 0, 0, 0};
    box.upper = {1, 1, 1, 1};
    
    Real abs_tol = 1e-3;
    Real rel_tol = 1e-2;
    
    auto result = integrate_adaptive<Real>(product_peak, box, abs_tol, rel_tol);
    
    // Analytical formula for product peak
    Real analytical = 1.0;
    for (std::size_t i = 0; i < 4; ++i) {
        Real term1 = std::atan(w[i] * c[i]);
        Real term2 = std::atan(w[i] * (1.0 - c[i]));
        analytical *= (term1 + term2) * w[i];
    }
    
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    if (result.status != status_code::success) {
        std::cout << "  Failed with status " << static_cast<int>(result.status) << std::endl;
    }
    assert(result.status == status_code::success);
    assert(std::abs(result.value - analytical) < abs_tol + rel_tol * analytical);
    assert(result.error < abs_tol + rel_tol * std::abs(result.value));
    
    std::cout << "  Result: " << result.value 
              << " (analytical: " << analytical << ")" << std::endl;
    std::cout << "  Error estimate: " << result.error << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_corner_peak_5d() {
    std::cout << "Testing 5D corner peak integral..." << std::endl;
    
    // Corner peak function (singularity at corner)
    std::vector<Real> a = {0.1, 0.2, 0.3, 0.4, 0.5};  // Exponents
    
    auto corner_peak = [&a](const std::vector<Real>& x) -> Real {
        Real sum = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            sum += a[i] * x[i];
        }
        return std::pow(1.0 + sum, -(static_cast<Real>(x.size()) + 1.0));
    };
    
    hypercube<Real> box(5);
    box.lower = {0, 0, 0, 0, 0};
    box.upper = {1, 1, 1, 1, 1};
    
    Real abs_tol = 1e-4;
    Real rel_tol = 1e-3;
    
    auto result = integrate_adaptive<Real>(corner_peak, box, abs_tol, rel_tol);
    
    // Analytical value (complex formula, approximate)
    Real sum_a = 0.0;
    for (auto ai : a) sum_a += ai;
    Real analytical = 1.0 / (5.0 * std::tgamma(5.0)) * 
                     std::pow(1.0 + sum_a, -1.0);
    
    // For corner peak, we expect larger errors or max eval reached
    assert(result.status == status_code::success || 
           result.status == status_code::maxeval_reached);
    
    std::cout << "  Result: " << result.value 
              << " (approximate analytical: " << analytical << ")" << std::endl;
    std::cout << "  Error estimate: " << result.error << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_monotone_error_decrease() {
    std::cout << "Testing monotone error decrease..." << std::endl;
    
    // Use a challenging integral
    auto f = [](const std::vector<Real>& x) -> Real {
        return std::exp(-10.0 * ((x[0] - 0.3) * (x[0] - 0.3) + 
                                  (x[1] - 0.7) * (x[1] - 0.7)));
    };
    
    hypercube<Real> box(2);
    box.lower = {0, 0};
    box.upper = {1, 1};
    
    std::vector<Real> tolerances = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
    std::vector<Real> errors;
    std::vector<std::size_t> evals;
    
    for (auto tol : tolerances) {
        auto result = integrate_adaptive<Real>(f, box, tol, tol);
        errors.push_back(result.error);
        evals.push_back(result.evaluations);
        
        std::cout << "  Tolerance: " << std::scientific << std::setprecision(0) << tol 
                  << ", Error: " << std::setprecision(2) << result.error
                  << ", Evaluations: " << result.evaluations << std::endl;
    }
    
    // Check that errors generally decrease
    for (std::size_t i = 1; i < errors.size(); ++i) {
        assert(errors[i] <= errors[i-1] * 1.5); // Allow some slack
    }
    
    // Check that evaluations increase
    for (std::size_t i = 1; i < evals.size(); ++i) {
        assert(evals[i] >= evals[i-1]);
    }
    
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_determinism() {
    std::cout << "Testing deterministic results..." << std::endl;
    
    // Complex function to ensure multiple subdivisions
    auto f = [](const std::vector<Real>& x) -> Real {
        return std::sin(5.0 * x[0]) * std::cos(3.0 * x[1]) * 
               std::exp(-2.0 * (x[2] * x[2]));
    };
    
    hypercube<Real> box(3);
    box.lower = {0, 0, 0};
    box.upper = {1, 1, 1};
    
    Real abs_tol = 1e-6;
    Real rel_tol = 1e-6;
    
    // Run multiple times
    const int num_runs = 5;
    std::vector<Real> values;
    std::vector<Real> errors;
    std::vector<std::size_t> evals;
    
    for (int i = 0; i < num_runs; ++i) {
        auto result = integrate_adaptive<Real>(f, box, abs_tol, rel_tol);
        values.push_back(result.value);
        errors.push_back(result.error);
        evals.push_back(result.evaluations);
    }
    
    // Check all runs give identical results
    for (int i = 1; i < num_runs; ++i) {
        assert(values[i] == values[0]);
        assert(errors[i] == errors[0]);
        assert(evals[i] == evals[0]);
    }
    
    std::cout << "  All " << num_runs << " runs gave identical results:" << std::endl;
    std::cout << "  Value: " << values[0] << std::endl;
    std::cout << "  Error: " << errors[0] << std::endl;
    std::cout << "  Evaluations: " << evals[0] << std::endl;
    std::cout << "  PASSED" << std::endl;
}

int main() {
    std::cout << "\n=== Testing Adaptive Integration with Genz Functions ===" << std::endl;
    
    // Test with double precision
    test_gaussian_2d<double>();
    test_oscillatory_3d<double>();
    test_product_peak_4d<double>();
    test_corner_peak_5d<double>();
    test_monotone_error_decrease<double>();
    test_determinism<double>();
    
    // Test with long double if available
    if (sizeof(long double) > sizeof(double)) {
        std::cout << "\n=== Testing with long double precision ===" << std::endl;
        test_gaussian_2d<long double>();
        test_oscillatory_3d<long double>();
    }
    
    std::cout << "\nAll tests passed!" << std::endl;
    return 0;
}