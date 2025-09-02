#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/math/cubature/cubature.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::cubature;
using namespace boost::math::constants;

// Example 1: Gaussian integral over entire real line
void example_gaussian_1d() {
    std::cout << "=== Example 1: Gaussian integral over (-∞, ∞) ===" << std::endl;
    
    using Real = double;
    
    // Integral of exp(-x²) from -∞ to ∞ = √π
    auto gaussian = [](const Real* x) -> Real {
        return std::exp(-x[0] * x[0]);
    };
    
    std::vector<Real> lower = {-std::numeric_limits<Real>::infinity()};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    auto result = integrate_adaptive_infinite<Real>(
        gaussian, lower, upper, 1e-10, 1e-10);
    
    Real exact = std::sqrt(pi<Real>());
    
    std::cout << "Computed: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "Exact:    " << exact << std::endl;
    std::cout << "Error:    " << std::abs(result.value - exact) << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    std::cout << std::endl;
}

// Example 2: Semi-infinite integral with exponential decay
void example_exponential_decay() {
    std::cout << "=== Example 2: Exponential decay over [0, ∞) ===" << std::endl;
    
    using Real = double;
    
    // Integral of x*exp(-x) from 0 to ∞ = 1
    auto integrand = [](const Real* x) -> Real {
        return x[0] * std::exp(-x[0]);
    };
    
    std::vector<Real> lower = {0};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    auto result = integrate_adaptive_infinite<Real>(
        integrand, lower, upper, 1e-10, 1e-10);
    
    Real exact = Real(1);
    
    std::cout << "Computed: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "Exact:    " << exact << std::endl;
    std::cout << "Error:    " << std::abs(result.value - exact) << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    std::cout << std::endl;
}

// Example 3: 2D Gaussian over infinite domain
void example_gaussian_2d() {
    std::cout << "=== Example 3: 2D Gaussian over (-∞, ∞)² ===" << std::endl;
    
    using Real = double;
    
    // Integral of exp(-(x²+y²)) over R² = π
    auto gaussian_2d = [](const Real* x) -> Real {
        return std::exp(-(x[0]*x[0] + x[1]*x[1]));
    };
    
    std::vector<Real> lower = {
        -std::numeric_limits<Real>::infinity(),
        -std::numeric_limits<Real>::infinity()
    };
    std::vector<Real> upper = {
        std::numeric_limits<Real>::infinity(),
        std::numeric_limits<Real>::infinity()
    };
    
    auto result = integrate_adaptive_infinite<Real>(
        gaussian_2d, lower, upper, 1e-8, 1e-8);
    
    Real exact = pi<Real>();
    
    std::cout << "Computed: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "Exact:    " << exact << std::endl;
    std::cout << "Error:    " << std::abs(result.value - exact) << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    std::cout << std::endl;
}

// Example 4: Mixed finite and infinite bounds
void example_mixed_bounds() {
    std::cout << "=== Example 4: Mixed bounds [0,1] × [0,∞) ===" << std::endl;
    
    using Real = double;
    
    // Integral of exp(-y) over [0,1] × [0,∞) = 1
    auto integrand = [](const Real* x) -> Real {
        return std::exp(-x[1]);
    };
    
    std::vector<Real> lower = {0, 0};
    std::vector<Real> upper = {1, std::numeric_limits<Real>::infinity()};
    
    auto result = integrate_adaptive_infinite<Real>(
        integrand, lower, upper, 1e-10, 1e-10);
    
    Real exact = Real(1);
    
    std::cout << "Computed: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "Exact:    " << exact << std::endl;
    std::cout << "Error:    " << std::abs(result.value - exact) << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    std::cout << std::endl;
}

// Example 5: Choosing between transform types
void example_transform_comparison() {
    std::cout << "=== Example 5: Transform comparison ===" << std::endl;
    
    using Real = double;
    
    // Integral with different decay rates
    // Fast decay: exp(-x²)
    auto fast_decay = [](const Real* x) -> Real {
        return std::exp(-x[0] * x[0]);
    };
    
    // Slow decay: 1/(1+x²)
    auto slow_decay = [](const Real* x) -> Real {
        return Real(1) / (Real(1) + x[0] * x[0]);
    };
    
    std::vector<Real> lower = {0};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    std::cout << "Fast decay (exp(-x²)):" << std::endl;
    auto result1 = integrate_adaptive_infinite<Real>(
        fast_decay, lower, upper, 1e-8, 1e-8);
    Real exact1 = std::sqrt(pi<Real>()) / Real(2);
    std::cout << "  Computed: " << result1.value << std::endl;
    std::cout << "  Exact: " << exact1 << std::endl;
    std::cout << "  Evaluations: " << result1.evaluations << std::endl;
    
    std::cout << "Slow decay (1/(1+x²)):" << std::endl;
    auto result2 = integrate_adaptive_infinite<Real>(
        slow_decay, lower, upper, 1e-8, 1e-8);
    Real exact2 = pi<Real>() / Real(2);
    std::cout << "  Computed: " << result2.value << std::endl;
    std::cout << "  Exact: " << exact2 << std::endl;
    std::cout << "  Evaluations: " << result2.evaluations << std::endl;
    std::cout << std::endl;
}

// Example 6: Using sparse grids for infinite domains
void example_sparse_grid_infinite() {
    std::cout << "=== Example 6: Sparse grid for infinite domain ===" << std::endl;
    
    using Real = double;
    
    // 3D Gaussian
    auto gaussian_3d = [](const Real* x) -> Real {
        return std::exp(-(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
    };
    
    std::vector<Real> lower(3, -std::numeric_limits<Real>::infinity());
    std::vector<Real> upper(3, std::numeric_limits<Real>::infinity());
    
    std::cout << "Using sparse grid level 3:" << std::endl;
    auto result = integrate_sparse_grid_infinite<Real>(
        gaussian_3d, lower, upper, 3);
    
    Real exact = std::pow(pi<Real>(), Real(1.5));
    
    std::cout << "Computed: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "Exact:    " << exact << std::endl;
    std::cout << "Error:    " << std::abs(result.value - exact) << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << std::fixed;
    
    example_gaussian_1d();
    example_exponential_decay();
    example_gaussian_2d();
    example_mixed_bounds();
    example_transform_comparison();
    example_sparse_grid_infinite();
    
    std::cout << "All examples completed!" << std::endl;
    
    return 0;
}