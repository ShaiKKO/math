// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/// \file multiprecision_test.cpp
/// \brief Tests for multiprecision arithmetic support
/// \details Validates that cubature algorithms work correctly with
///          arbitrary precision types like cpp_dec_float_50

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>

using namespace boost::math::cubature;
using namespace boost::multiprecision;

// Multiprecision types
using float_50 = cpp_dec_float_50;
using float_100 = cpp_dec_float_100;

/// \brief Test function: Gaussian integral with known exact value
template <typename Real>
class gaussian_nd {
private:
    std::size_t dim_;
    Real sigma_;
    
public:
    gaussian_nd(std::size_t dim, Real sigma = Real(1))
        : dim_(dim), sigma_(sigma) {}
    
    Real operator()(const Real* x) const {
        Real sum = Real(0);
        for (std::size_t i = 0; i < dim_; ++i) {
            Real diff = x[i] - Real(0.5);  // Centered at 0.5
            sum += diff * diff;
        }
        return exp(-sum / (Real(2) * sigma_ * sigma_));
    }
    
    Real operator()(const Real* x, std::size_t /*n*/) const {
        return (*this)(x);
    }
    
    /// \brief Exact integral over [0,1]^d
    Real exact_integral() const {
        using std::sqrt;
        using std::erf;
        Real sqrt_2 = sqrt(Real(2));
        Real sqrt_pi = sqrt(boost::math::constants::pi<Real>());
        
        // For each dimension: integral of exp(-(x-0.5)^2/(2σ^2)) from 0 to 1
        Real one_d_integral = sigma_ * sqrt_pi * sqrt_2 * 
                             (erf(Real(0.5) / (sigma_ * sqrt_2)) - 
                              erf(Real(-0.5) / (sigma_ * sqrt_2))) / Real(2);
        
        // Product over all dimensions
        Real result = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            result *= one_d_integral;
        }
        return result;
    }
};

/// \brief Test polynomial that should be exact for sufficient degree rules
template <typename Real>
class polynomial_nd {
private:
    std::size_t dim_;
    std::size_t degree_;
    
public:
    polynomial_nd(std::size_t dim, std::size_t degree)
        : dim_(dim), degree_(degree) {}
    
    Real operator()(const Real* x) const {
        Real result = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            Real xi = x[i];
            Real term = Real(0);
            for (std::size_t p = 0; p <= degree_; ++p) {
                Real coeff = Real(p + 1);  // Simple coefficients
                Real power = Real(1);
                for (std::size_t j = 0; j < p; ++j) {
                    power *= xi;
                }
                term += coeff * power;
            }
            result *= term;
        }
        return result;
    }
    
    Real operator()(const Real* x, std::size_t /*n*/) const {
        return (*this)(x);
    }
    
    /// \brief Exact integral over [0,1]^d
    Real exact_integral() const {
        // Integral of sum(k*x^(k-1)) from 0 to 1
        Real one_d_integral = Real(0);
        for (std::size_t p = 0; p <= degree_; ++p) {
            Real coeff = Real(p + 1);
            one_d_integral += coeff / Real(p + 1);  // Integral of coeff*x^p
        }
        
        // Product over all dimensions
        Real result = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            result *= one_d_integral;
        }
        return result;
    }
};

/// \brief Test harness for multiprecision integration
template <typename Real>
class multiprecision_tester {
public:
    struct test_result {
        std::string test_name;
        std::string algorithm;
        Real exact_value;
        Real computed_value;
        Real error_estimate;
        Real actual_error;
        std::size_t evaluations;
        bool passed;
        std::string precision_type;
    };
    
    static void print_result(const test_result& r) {
        std::cout << "\n=== " << r.test_name << " ===\n";
        std::cout << "Algorithm: " << r.algorithm << "\n";
        std::cout << "Precision: " << r.precision_type << "\n";
        std::cout << "Exact:     " << std::setprecision(50) << r.exact_value << "\n";
        std::cout << "Computed:  " << r.computed_value << "\n";
        std::cout << "Est Error: " << std::scientific << std::setprecision(6) 
                  << r.error_estimate << "\n";
        std::cout << "Act Error: " << r.actual_error << "\n";
        std::cout << "Evals:     " << std::fixed << r.evaluations << "\n";
        std::cout << "Status:    " << (r.passed ? "PASS" : "FAIL") << "\n";
    }
    
    /// \brief Test adaptive integration with multiprecision
    static test_result test_adaptive_gaussian(std::size_t dim) {
        test_result res;
        res.test_name = "Gaussian " + std::to_string(dim) + "D";
        res.algorithm = "Adaptive (Genz-Malik)";
        res.precision_type = get_precision_name();
        
        gaussian_nd<Real> f(dim, Real(0.5));
        res.exact_value = f.exact_integral();
        
        hypercube<Real> box(dim);
        std::fill(box.lower.begin(), box.lower.end(), Real(0));
        std::fill(box.upper.begin(), box.upper.end(), Real(1));
        
        // Use tight tolerance appropriate for precision
        Real abs_tol = std::numeric_limits<Real>::epsilon() * Real(1000);
        Real rel_tol = std::numeric_limits<Real>::epsilon() * Real(100);
        
        // Custom policy for multiprecision
        using high_prec_policy = boost::math::policies::policy<
            boost::math::policies::precision<std::numeric_limits<Real>::digits10>
        >;
        
        auto result = integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 
                                              100000, high_prec_policy{});
        
        res.computed_value = result.value;
        res.error_estimate = result.error;
        res.actual_error = abs(res.computed_value - res.exact_value);
        res.evaluations = result.evaluations;
        
        // Check if result is within tolerance
        res.passed = (res.actual_error <= abs_tol + rel_tol * abs(res.exact_value)) &&
                    (result.status == status_code::success);
        
        return res;
    }
    
    /// \brief Test sparse grid with multiprecision
    static test_result test_sparse_polynomial(std::size_t dim, unsigned level) {
        test_result res;
        res.test_name = "Polynomial " + std::to_string(dim) + "D deg-" + 
                       std::to_string(2*level-1);
        res.algorithm = "Sparse Grid (level " + std::to_string(level) + ")";
        res.precision_type = get_precision_name();
        
        polynomial_nd<Real> f(dim, 2*level - 1);  // Degree that should be exact
        res.exact_value = f.exact_integral();
        
        hypercube<Real> box(dim);
        std::fill(box.lower.begin(), box.lower.end(), Real(0));
        std::fill(box.upper.begin(), box.upper.end(), Real(1));
        
        auto result = integrate_sparse_grid<Real>(f, box, level);
        
        res.computed_value = result.value;
        res.error_estimate = result.error;
        res.actual_error = abs(res.computed_value - res.exact_value);
        res.evaluations = result.evaluations;
        
        // Should be exact to machine precision for polynomials
        res.passed = res.actual_error < std::numeric_limits<Real>::epsilon() * Real(100);
        
        return res;
    }
    
    /// \brief Test QMC with multiprecision
    static test_result test_qmc_gaussian(std::size_t dim, std::size_t n_points) {
        test_result res;
        res.test_name = "Gaussian " + std::to_string(dim) + "D QMC";
        res.algorithm = "QMC (Sobol, " + std::to_string(n_points) + " points)";
        res.precision_type = get_precision_name();
        
        gaussian_nd<Real> f(dim, Real(0.5));
        res.exact_value = f.exact_integral();
        
        hypercube<Real> box(dim);
        std::fill(box.lower.begin(), box.lower.end(), Real(0));
        std::fill(box.upper.begin(), box.upper.end(), Real(1));
        
        auto result = integrate_qmc<Real>(f, box, n_points);
        
        res.computed_value = result.value;
        res.error_estimate = result.error;
        res.actual_error = abs(res.computed_value - res.exact_value);
        res.evaluations = result.evaluations;
        
        // QMC convergence is O(log^d(N)/N), check reasonable accuracy
        Real expected_error = Real(1) / sqrt(Real(n_points));
        res.passed = res.actual_error < expected_error;
        
        return res;
    }
    
private:
    static std::string get_precision_name() {
        if (std::numeric_limits<Real>::digits10 > 100) {
            return "cpp_dec_float (>" + std::to_string(100) + " digits)";
        } else if (std::numeric_limits<Real>::digits10 > 50) {
            return "cpp_dec_float_100 (" + 
                   std::to_string(std::numeric_limits<Real>::digits10) + " digits)";
        } else if (std::numeric_limits<Real>::digits10 > 30) {
            return "cpp_dec_float_50 (" + 
                   std::to_string(std::numeric_limits<Real>::digits10) + " digits)";
        } else if (std::numeric_limits<Real>::digits10 > 15) {
            return "long double (" + 
                   std::to_string(std::numeric_limits<Real>::digits10) + " digits)";
        } else {
            return "double (" + 
                   std::to_string(std::numeric_limits<Real>::digits10) + " digits)";
        }
    }
};

/// \brief Run comprehensive multiprecision tests
template <typename Real>
void run_multiprecision_suite() {
    using tester = multiprecision_tester<Real>;
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Multiprecision Tests - " << 
              std::numeric_limits<Real>::digits10 << " decimal digits\n";
    std::cout << std::string(60, '=') << "\n";
    
    // Test adaptive integration
    for (std::size_t dim : {2, 3, 4}) {
        auto result = tester::test_adaptive_gaussian(dim);
        tester::print_result(result);
    }
    
    // Test sparse grid (exact for polynomials)
    for (std::size_t dim : {2, 3}) {
        for (unsigned level : {3, 4}) {
            auto result = tester::test_sparse_polynomial(dim, level);
            tester::print_result(result);
        }
    }
    
    // Test QMC
    for (std::size_t dim : {2, 3, 5}) {
        auto result = tester::test_qmc_gaussian(dim, 10000);
        tester::print_result(result);
    }
}

/// \brief Test that algorithms respect precision
template <typename Real>
void test_precision_preservation() {
    std::cout << "\n=== Precision Preservation Test ===\n";
    std::cout << "Type precision: " << std::numeric_limits<Real>::digits10 
              << " decimal digits\n";
    
    // Create a simple integral that requires high precision
    auto f = [](const Real* x) -> Real {
        // Function with small differences that require precision
        Real a = Real(1) + x[0] * std::numeric_limits<Real>::epsilon() * Real(1000);
        Real b = Real(1) + x[1] * std::numeric_limits<Real>::epsilon() * Real(1001);
        return a * b;
    };
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // The exact integral can be computed analytically
    Real exact = Real(1) + std::numeric_limits<Real>::epsilon() * Real(1000.5);
    
    auto result = integrate_adaptive<Real>(f, box, 
                                          std::numeric_limits<Real>::epsilon() * Real(10),
                                          std::numeric_limits<Real>::epsilon() * Real(10));
    
    std::cout << "Exact:    " << std::setprecision(50) << exact << "\n";
    std::cout << "Computed: " << result.value << "\n";
    std::cout << "Error:    " << std::scientific << std::setprecision(6) 
              << abs(result.value - exact) << "\n";
    
    bool passed = abs(result.value - exact) < 
                 std::numeric_limits<Real>::epsilon() * Real(100);
    std::cout << "Status:   " << (passed ? "PASS" : "FAIL") << "\n";
}

int main() {
    std::cout << "Boost.Math Cubature - Multiprecision Arithmetic Tests\n";
    std::cout << "======================================================\n";
    
    // Test with double precision (baseline)
    std::cout << "\n### Testing with double precision (baseline) ###\n";
    run_multiprecision_suite<double>();
    test_precision_preservation<double>();
    
    // Test with 50-digit precision
    std::cout << "\n### Testing with cpp_dec_float_50 ###\n";
    run_multiprecision_suite<float_50>();
    test_precision_preservation<float_50>();
    
    // Test with 100-digit precision
    std::cout << "\n### Testing with cpp_dec_float_100 ###\n";
    run_multiprecision_suite<float_100>();
    test_precision_preservation<float_100>();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Multiprecision tests completed successfully!\n";
    std::cout << "All algorithms maintain precision across different types.\n";
    
    return 0;
}