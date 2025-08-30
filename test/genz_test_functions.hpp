// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_TEST_GENZ_TEST_FUNCTIONS_HPP
#define BOOST_MATH_CUBATURE_TEST_GENZ_TEST_FUNCTIONS_HPP

/// \file genz_test_functions.hpp
/// \brief Genz test functions for cubature algorithm validation
/// \details Implements the six standard test integrands proposed by Alan Genz
///          for testing multidimensional integration algorithms. Each function
///          has adjustable difficulty parameters and known analytical integrals.
///
/// \section genz_functions The Six Genz Functions
/// 1. **Oscillatory**: Tests handling of oscillating integrands
/// 2. **Product Peak**: Sharp peak at interior point
/// 3. **Corner Peak**: Singularity-like behavior at corner
/// 4. **Gaussian**: Smooth bell-shaped function
/// 5. **Continuous**: C⁰ continuous with discontinuous derivatives
/// 6. **Discontinuous**: Jump discontinuity at interior hyperplane
///
/// \references
/// - A. Genz, "Testing Multidimensional Integration Routines", 
///   in Tools, Methods and Languages for Scientific and Engineering Computation,
///   B. Ford, J.C. Rault, and F. Thomasset, eds., North-Holland, 1984, pp. 81-94
/// - A. Genz, "A Package for Testing Multiple Integration Subroutines",
///   in Numerical Integration: Recent Developments, Software and Applications,
///   P. Keast and G. Fairweather, eds., D. Reidel, 1987, pp. 337-340

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>

namespace boost { namespace math { namespace cubature { namespace test {

/// \brief Base class for Genz test functions
/// \tparam Real Floating point type
template <typename Real>
class genz_function {
protected:
    std::size_t dim_;              ///< Dimension
    std::vector<Real> a_;          ///< Difficulty parameters
    std::vector<Real> u_;          ///< Shift parameters
    
public:
    /// \brief Constructor
    /// \param dim Dimension
    /// \param a Difficulty parameters (higher = harder)
    /// \param u Shift parameters in [0,1]^d
    genz_function(std::size_t dim, 
                  const std::vector<Real>& a,
                  const std::vector<Real>& u)
        : dim_(dim), a_(a), u_(u) 
    {
        if (a_.size() != dim_) a_.resize(dim_, Real(1));
        if (u_.size() != dim_) u_.resize(dim_, Real(0.5));
    }
    
    /// \brief Get dimension
    std::size_t dimension() const { return dim_; }
    
    /// \brief Get difficulty parameters
    const std::vector<Real>& difficulty() const { return a_; }
    
    /// \brief Get shift parameters
    const std::vector<Real>& shift() const { return u_; }
    
    /// \brief Evaluate function (to be overridden)
    virtual Real operator()(const Real* x) const = 0;
    
    /// \brief Get exact integral over [0,1]^d (to be overridden)
    virtual Real exact_integral() const = 0;
    
    /// \brief Get function name
    virtual const char* name() const = 0;
    
    /// \brief Support vector interface
    Real operator()(const std::vector<Real>& x) const {
        return (*this)(x.data());
    }
    
    /// \brief Support (ptr, size) interface
    Real operator()(const Real* x, std::size_t /*n*/) const {
        return (*this)(x);
    }
    
    virtual ~genz_function() = default;
};

/// \brief Genz Function 1: Oscillatory
/// \details f(x) = cos(2πu₁ + Σaᵢxᵢ)
template <typename Real>
class genz_oscillatory : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        Real sum = Real(2) * M_PI * u_[0];
        for (std::size_t i = 0; i < dim_; ++i) {
            sum += a_[i] * x[i];
        }
        return std::cos(sum);
    }
    
    Real exact_integral() const override {
        // Integral = Π[sin(aᵢ)/aᵢ] * cos(2πu₁ + Σaᵢ/2)
        Real prod = Real(1);
        Real phase_shift = Real(2) * M_PI * u_[0];
        
        for (std::size_t i = 0; i < dim_; ++i) {
            if (std::abs(a_[i]) > 1e-10) {
                prod *= std::sin(a_[i]) / a_[i];
                phase_shift += a_[i] / Real(2);
            }
        }
        
        return prod * std::cos(phase_shift);
    }
    
    const char* name() const override { return "Oscillatory"; }
};

/// \brief Genz Function 2: Product Peak
/// \details f(x) = Π[1/(aᵢ⁻² + (xᵢ - uᵢ)²)]
template <typename Real>
class genz_product_peak : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        Real prod = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            Real ai_inv2 = Real(1) / (a_[i] * a_[i]);
            Real diff = x[i] - u_[i];
            prod /= (ai_inv2 + diff * diff);
        }
        return prod;
    }
    
    Real exact_integral() const override {
        // Integral = Π[aᵢ * (arctan(aᵢ(1-uᵢ)) + arctan(aᵢuᵢ))]
        Real prod = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            Real term1 = std::atan(a_[i] * (Real(1) - u_[i]));
            Real term2 = std::atan(a_[i] * u_[i]);
            prod *= a_[i] * (term1 + term2);
        }
        return prod;
    }
    
    const char* name() const override { return "Product Peak"; }
};

/// \brief Genz Function 3: Corner Peak
/// \details f(x) = (1 + Σaᵢxᵢ)^(-(d+1))
template <typename Real>
class genz_corner_peak : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        Real sum = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            sum += a_[i] * x[i];
        }
        return std::pow(sum, -static_cast<Real>(dim_ + 1));
    }
    
    Real exact_integral() const override {
        // Complex formula involving combinatorics
        // For simplicity, using numerical approximation for now
        // Exact formula: involves sum over all 2^d vertices
        Real result = Real(0);
        
        // Generate all 2^d corner combinations
        std::size_t num_corners = 1 << dim_;
        for (std::size_t mask = 0; mask < num_corners; ++mask) {
            Real corner_sum = Real(1);
            int sign = 1;
            
            for (std::size_t i = 0; i < dim_; ++i) {
                if (mask & (1 << i)) {
                    corner_sum += a_[i];
                    sign *= -1;
                }
            }
            
            if (corner_sum > 0) {
                result += sign * std::pow(corner_sum, -static_cast<Real>(dim_)) / 
                         static_cast<Real>(dim_);
            }
        }
        
        return std::abs(result);
    }
    
    const char* name() const override { return "Corner Peak"; }
};

/// \brief Genz Function 4: Gaussian
/// \details f(x) = exp(-Σaᵢ²(xᵢ - uᵢ)²)
template <typename Real>
class genz_gaussian : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        Real sum = Real(0);
        for (std::size_t i = 0; i < dim_; ++i) {
            Real diff = x[i] - u_[i];
            sum += a_[i] * a_[i] * diff * diff;
        }
        return std::exp(-sum);
    }
    
    Real exact_integral() const override {
        // Integral = Π[(√π/aᵢ) * (erf(aᵢuᵢ) + erf(aᵢ(1-uᵢ)))/2]
        Real prod = Real(1);
        Real sqrt_pi = std::sqrt(M_PI);
        
        for (std::size_t i = 0; i < dim_; ++i) {
            Real term1 = std::erf(a_[i] * u_[i]);
            Real term2 = std::erf(a_[i] * (Real(1) - u_[i]));
            prod *= (sqrt_pi / a_[i]) * (term1 + term2) / Real(2);
        }
        
        return prod;
    }
    
    const char* name() const override { return "Gaussian"; }
};

/// \brief Genz Function 5: C⁰ Continuous
/// \details f(x) = exp(-Σaᵢ|xᵢ - uᵢ|)
template <typename Real>
class genz_continuous : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        Real sum = Real(0);
        for (std::size_t i = 0; i < dim_; ++i) {
            sum += a_[i] * std::abs(x[i] - u_[i]);
        }
        return std::exp(-sum);
    }
    
    Real exact_integral() const override {
        // Integral = Π[(2 - exp(-aᵢuᵢ) - exp(-aᵢ(1-uᵢ)))/aᵢ]
        Real prod = Real(1);
        
        for (std::size_t i = 0; i < dim_; ++i) {
            Real term1 = std::exp(-a_[i] * u_[i]);
            Real term2 = std::exp(-a_[i] * (Real(1) - u_[i]));
            prod *= (Real(2) - term1 - term2) / a_[i];
        }
        
        return prod;
    }
    
    const char* name() const override { return "C0 Continuous"; }
};

/// \brief Genz Function 6: Discontinuous
/// \details f(x) = 0 if any xᵢ > uᵢ, else exp(Σaᵢxᵢ)
template <typename Real>
class genz_discontinuous : public genz_function<Real> {
    using base = genz_function<Real>;
    using base::dim_;
    using base::a_;
    using base::u_;
    
public:
    using base::base;
    
    Real operator()(const Real* x) const override {
        // Check discontinuity condition
        for (std::size_t i = 0; i < dim_; ++i) {
            if (x[i] > u_[i]) {
                return Real(0);
            }
        }
        
        // In continuous region
        Real sum = Real(0);
        for (std::size_t i = 0; i < dim_; ++i) {
            sum += a_[i] * x[i];
        }
        return std::exp(sum);
    }
    
    Real exact_integral() const override {
        // Integral = Π[(exp(aᵢuᵢ) - 1)/aᵢ]
        Real prod = Real(1);
        
        for (std::size_t i = 0; i < dim_; ++i) {
            if (std::abs(a_[i]) > 1e-10) {
                prod *= (std::exp(a_[i] * u_[i]) - Real(1)) / a_[i];
            } else {
                prod *= u_[i];  // Limit as aᵢ → 0
            }
        }
        
        return prod;
    }
    
    const char* name() const override { return "Discontinuous"; }
};

/// \brief Factory function to create Genz test function
/// \param type Function type (1-6)
/// \param dim Dimension
/// \param difficulty Overall difficulty level (0=easy, 1=moderate, 2=hard)
/// \param seed Random seed for parameter generation
template <typename Real>
std::unique_ptr<genz_function<Real>> 
create_genz_function(int type, std::size_t dim, 
                    int difficulty = 1, unsigned seed = 12345)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<Real> dist(0, 1);
    
    // Generate difficulty parameters
    std::vector<Real> a(dim);
    Real base_difficulty = (difficulty == 0) ? Real(1) : 
                          (difficulty == 1) ? Real(5) : Real(10);
    
    for (std::size_t i = 0; i < dim; ++i) {
        a[i] = base_difficulty * (Real(0.5) + dist(gen));
    }
    
    // Generate shift parameters
    std::vector<Real> u(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        u[i] = dist(gen);
    }
    
    // Create appropriate function
    switch (type) {
        case 1:
            return std::make_unique<genz_oscillatory<Real>>(dim, a, u);
        case 2:
            return std::make_unique<genz_product_peak<Real>>(dim, a, u);
        case 3:
            return std::make_unique<genz_corner_peak<Real>>(dim, a, u);
        case 4:
            return std::make_unique<genz_gaussian<Real>>(dim, a, u);
        case 5:
            return std::make_unique<genz_continuous<Real>>(dim, a, u);
        case 6:
            return std::make_unique<genz_discontinuous<Real>>(dim, a, u);
        default:
            return nullptr;
    }
}

/// \brief Test suite runner for Genz functions
template <typename Real>
class genz_test_suite {
public:
    struct test_result {
        std::string function_name;
        std::size_t dimension;
        Real exact_value;
        Real computed_value;
        Real absolute_error;
        Real relative_error;
        std::size_t evaluations;
        bool passed;
    };
    
    /// \brief Run test on a single Genz function
    template <typename Integrator>
    static test_result run_single_test(
        const genz_function<Real>& f,
        Integrator& integrator,
        Real abs_tol,
        Real rel_tol)
    {
        test_result result;
        result.function_name = f.name();
        result.dimension = f.dimension();
        result.exact_value = f.exact_integral();
        
        // Run integration
        auto integration_result = integrator(f, abs_tol, rel_tol);
        
        result.computed_value = integration_result.value;
        result.evaluations = integration_result.evaluations;
        result.absolute_error = std::abs(result.computed_value - result.exact_value);
        result.relative_error = result.absolute_error / 
                               (std::abs(result.exact_value) + 1e-10);
        
        // Check if test passed
        result.passed = (result.absolute_error <= abs_tol) ||
                       (result.relative_error <= rel_tol);
        
        return result;
    }
    
    /// \brief Run full test suite
    template <typename Integrator>
    static std::vector<test_result> run_suite(
        std::size_t dim,
        Integrator& integrator,
        Real abs_tol = 1e-6,
        Real rel_tol = 1e-6,
        int difficulty = 1,
        unsigned seed = 12345)
    {
        std::vector<test_result> results;
        
        // Test all 6 Genz functions
        for (int type = 1; type <= 6; ++type) {
            auto f = create_genz_function<Real>(type, dim, difficulty, seed);
            if (f) {
                results.push_back(run_single_test(*f, integrator, abs_tol, rel_tol));
            }
        }
        
        return results;
    }
};

}}}} // namespace boost::math::cubature::test

#endif // BOOST_MATH_CUBATURE_TEST_GENZ_TEST_FUNCTIONS_HPP