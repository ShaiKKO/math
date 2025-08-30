// Test QMC vector integration

#include <boost/math/cubature/qmc.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// Define this to use the placeholder Sobol implementation
#ifndef BOOST_MATH_HAS_BOOST_RANDOM
#define BOOST_MATH_HAS_BOOST_RANDOM
// Placeholder Sobol implementation for testing
namespace boost { namespace random {
    template <typename UIntType, std::size_t Dimension>
    class sobol_engine {
    public:
        explicit sobol_engine(std::size_t) {}
        std::vector<UIntType> operator()() {
            static std::mt19937_64 gen(42);
            std::uniform_int_distribution<UIntType> dist;
            std::vector<UIntType> result(Dimension);
            for (auto& val : result) {
                val = dist(gen);
            }
            return result;
        }
    };
}} // namespace boost::random
#endif

using namespace boost::math::cubature;
using namespace boost::math::constants;

// Vector test integrand: [cos(πx₀/2), sin(πx₁/2), exp(-||x||²)]
template <typename Real>
void vector_integrand(const Real* x, Real* out, std::size_t num_components) {
    if (num_components >= 1) {
        // Component 0: cos(π*x₀/2)
        out[0] = std::cos(pi<Real>() * x[0] / 2);
    }
    
    if (num_components >= 2) {
        // Component 1: sin(π*x₁/2) if 2D, else 1
        Real dim = 2; // Assuming 2D for this test
        out[1] = (dim >= 2) ? std::sin(pi<Real>() * x[1] / 2) : Real(1);
    }
    
    if (num_components >= 3) {
        // Component 2: exp(-||x||²) where x is mapped from [0,1] to [-1,1]
        Real sum = 0;
        Real dim = 2; // Assuming 2D for this test
        for (std::size_t i = 0; i < dim; ++i) {
            Real xi = 2*x[i] - 1;
            sum += xi * xi;
        }
        out[2] = std::exp(-sum);
    }
}

// Polynomial vector integrand for testing
template <typename Real>
void polynomial_vector(const Real* x, Real* out, std::size_t num_components) {
    if (num_components >= 1) out[0] = x[0] * x[0];           // x₀²
    if (num_components >= 2) out[1] = x[0] * x[1];           // x₀x₁
    if (num_components >= 3) out[2] = x[0] + x[1];           // x₀ + x₁
    if (num_components >= 4) out[3] = Real(1);               // constant
}

template <typename Real>
void test_qmc_vector(std::size_t dim, std::size_t n_points) {
    std::cout << "\n=== Testing QMC Vector Integration ===" << std::endl;
    std::cout << "Dimension: " << dim << ", Points: " << n_points << std::endl;
    
    // Create integration domain [0,1]^dim
    hypercube<Real> box(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        box.lower[i] = 0;
        box.upper[i] = 1;
    }
    
    // Test vector integrand with 3 components
    {
        std::cout << "\nVector integrand test (3 components):" << std::endl;
        
        // Plain QMC
        std::cout << "Plain QMC:" << std::endl;
        auto results = integrate_qmc_vector<Real>(
            vector_integrand<Real>, box, 3, n_points);
        
        for (std::size_t i = 0; i < 3; ++i) {
            std::cout << "  Component " << i << ": value=" << std::setprecision(10) 
                      << results[i].value << ", error=" << std::scientific 
                      << results[i].error << ", evals=" << results[i].evaluations << std::endl;
        }
        
        // RQMC with variance estimation
        std::cout << "\nRandomized QMC (10 replicates):" << std::endl;
        auto rqmc_results = integrate_rqmc_vector<Real>(
            vector_integrand<Real>, box, 3, n_points, 10);
        
        for (std::size_t i = 0; i < 3; ++i) {
            std::cout << "  Component " << i << ": value=" << std::setprecision(10) 
                      << rqmc_results[i].value << ", error=" << std::scientific 
                      << rqmc_results[i].error << ", evals=" << rqmc_results[i].evaluations << std::endl;
        }
    }
    
    // Test polynomial vector with 4 components
    {
        std::cout << "\nPolynomial vector test (4 components):" << std::endl;
        auto results = integrate_qmc_vector<Real>(
            polynomial_vector<Real>, box, 4, n_points);
        
        // Exact values for polynomials over [0,1]²
        // Component 0: ∫x₀² = 1/3
        // Component 1: ∫x₀x₁ = 1/4
        // Component 2: ∫(x₀ + x₁) = 1
        // Component 3: ∫1 = 1
        Real exact[4] = {Real(1)/3, Real(1)/4, Real(1), Real(1)};
        
        for (std::size_t i = 0; i < 4; ++i) {
            std::cout << "Component " << i << ":" << std::endl;
            std::cout << "  Result: " << std::setprecision(10) << results[i].value << std::endl;
            std::cout << "  Exact:  " << exact[i] << std::endl;
            std::cout << "  Error:  " << std::scientific 
                      << std::abs(results[i].value - exact[i]) << std::endl;
        }
        std::cout << "  Evaluations: " << results[0].evaluations << std::endl;
    }
}

template <typename Real>
void test_rqmc_convergence() {
    std::cout << "\n=== Testing RQMC Convergence ===" << std::endl;
    
    hypercube<Real> box(2);
    box.lower[0] = box.lower[1] = 0;
    box.upper[0] = box.upper[1] = 1;
    
    // Simple vector integrand: [x₀, x₁]
    auto f = [](const Real* x, Real* out, std::size_t m) {
        if (m >= 1) out[0] = x[0];
        if (m >= 2) out[1] = x[1];
    };
    
    std::cout << "\nExact values: [0.5, 0.5]" << std::endl;
    std::cout << "\nPoints\tComponent 0\t\tComponent 1\t\tError Est 0\tError Est 1" << std::endl;
    
    for (std::size_t n_points : {16, 64, 256, 1024}) {
        auto results = integrate_rqmc_vector<Real>(f, box, 2, n_points, 20);
        
        std::cout << n_points << "\t" 
                  << std::fixed << std::setprecision(10) << results[0].value << "\t"
                  << results[1].value << "\t"
                  << std::scientific << results[0].error << "\t"
                  << results[1].error << std::endl;
    }
}

// Test comparison between scalar and vector QMC
template <typename Real>
void test_scalar_vector_consistency() {
    std::cout << "\n=== Testing Scalar vs Vector Consistency ===" << std::endl;
    
    hypercube<Real> box(3);
    for (std::size_t i = 0; i < 3; ++i) {
        box.lower[i] = 0;
        box.upper[i] = 1;
    }
    
    // Scalar integrand
    auto scalar_f = [](const Real* x, std::size_t) {
        return std::exp(-x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
    };
    
    // Vector integrand with same function as first component
    auto vector_f = [](const Real* x, Real* out, std::size_t m) {
        if (m >= 1) out[0] = std::exp(-x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
        if (m >= 2) out[1] = x[0];  // Different component
    };
    
    std::size_t n_points = 1000;
    
    // Scalar QMC
    auto scalar_result = integrate_qmc<Real>(scalar_f, box, n_points);
    std::cout << "Scalar QMC: " << std::setprecision(10) << scalar_result.value << std::endl;
    
    // Vector QMC (first component should match)
    auto vector_results = integrate_qmc_vector<Real>(vector_f, box, 2, n_points);
    std::cout << "Vector QMC (component 0): " << vector_results[0].value << std::endl;
    std::cout << "Difference: " << std::scientific 
              << std::abs(scalar_result.value - vector_results[0].value) << std::endl;
    
    // RQMC comparison
    std::cout << "\nRQMC comparison (10 replicates):" << std::endl;
    auto scalar_rqmc = integrate_rqmc<Real>(scalar_f, box, n_points, 10);
    auto vector_rqmc = integrate_rqmc_vector<Real>(vector_f, box, 2, n_points, 10);
    
    std::cout << "Scalar RQMC: " << std::setprecision(10) << scalar_rqmc.value 
              << " ± " << std::scientific << scalar_rqmc.error << std::endl;
    std::cout << "Vector RQMC (component 0): " << std::setprecision(10) << vector_rqmc[0].value 
              << " ± " << std::scientific << vector_rqmc[0].error << std::endl;
}

int main() {
    std::cout << "===== QMC Vector Integration Tests =====\n" << std::endl;
    
    using Real = double;
    
    // Test basic QMC vector integration
    test_qmc_vector<Real>(2, 100);
    test_qmc_vector<Real>(3, 500);
    
    // Test RQMC convergence
    test_rqmc_convergence<Real>();
    
    // Test scalar vs vector consistency
    test_scalar_vector_consistency<Real>();
    
    std::cout << "\n===== Tests Complete =====\n" << std::endl;
    
    return 0;
}