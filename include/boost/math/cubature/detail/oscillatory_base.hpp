// Copyright 2025 Boost Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_OSCILLATORY_BASE_HPP
#define BOOST_MATH_CUBATURE_DETAIL_OSCILLATORY_BASE_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/precision.hpp>
#include <complex>
#include <functional>
#include <cmath>
#include <limits>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Types of oscillatory functions supported
enum class oscillator_type {
    exponential,  // e^(iωg(x))
    cosine,       // cos(ωg(x))
    sine,         // sin(ωg(x))
    bessel_j      // J_n(ωg(x))
};

/// \brief Method to use for oscillatory integration
enum class oscillatory_method {
    automatic,     // Choose automatically based on problem
    filon,        // Filon quadrature
    levin,        // Levin collocation method
    asymptotic    // Asymptotic expansion (for very large ω)
};

/// \brief Base structure for oscillatory integrals
/// 
/// Represents integrals of the form:
/// ∫ f(x) * oscillator(ω*g(x)) dx
///
template <typename Real>
struct oscillatory_integral {
    using complex_type = std::complex<Real>;
    using amplitude_function = std::function<Real(Real)>;
    using phase_function = std::function<Real(Real)>;
    
    amplitude_function f;        // Amplitude function f(x)
    phase_function g;            // Phase function g(x)
    Real omega;                  // Frequency parameter ω
    oscillator_type type;        // Type of oscillator
    int bessel_order{0};        // Order for Bessel functions
    
    // Interval of integration
    Real a, b;
    
    // Constructor for standard oscillators
    oscillatory_integral(
        amplitude_function amplitude,
        phase_function phase,
        Real frequency,
        Real lower, Real upper,
        oscillator_type osc_type = oscillator_type::exponential)
        : f(std::move(amplitude))
        , g(std::move(phase))
        , omega(frequency)
        , type(osc_type)
        , a(lower)
        , b(upper)
    {
        validate();
    }
    
    // Evaluate the oscillatory part at x
    complex_type oscillator_value(Real x) const {
        Real phase_val = omega * g(x);
        
        switch(type) {
            case oscillator_type::exponential:
                return std::exp(complex_type(0, phase_val));
            case oscillator_type::cosine:
                return complex_type(std::cos(phase_val), 0);
            case oscillator_type::sine:
                return complex_type(std::sin(phase_val), 0);
            case oscillator_type::bessel_j:
                return complex_type(boost::math::cyl_bessel_j(bessel_order, phase_val), 0);
            default:
                return complex_type(1, 0);
        }
    }
    
    // Evaluate the full integrand
    complex_type integrand(Real x) const {
        return f(x) * oscillator_value(x);
    }
    
    // Check if frequency is high enough for oscillatory methods
    bool is_highly_oscillatory() const {
        // Estimate number of oscillations in interval
        Real path_length = std::abs(omega * (g(b) - g(a)));
        return path_length > Real(10) * boost::math::constants::pi<Real>();
    }
    
    // Estimate optimal method based on problem characteristics
    oscillatory_method suggest_method() const {
        if (!is_highly_oscillatory()) {
            return oscillatory_method::automatic; // Use standard quadrature
        }
        
        // For linear phase (g(x) = x), Filon is efficient
        // Test linearity by sampling
        Real mid = (a + b) / 2;
        Real g_mid_expected = (g(a) + g(b)) / 2;
        Real g_mid_actual = g(mid);
        Real linearity_error = std::abs(g_mid_actual - g_mid_expected) / (std::abs(g_mid_expected) + 1);
        
        if (linearity_error < Real(0.01)) {
            return oscillatory_method::filon;
        }
        
        // For very high frequencies, asymptotic methods may be best
        if (omega > Real(1000)) {
            // Check if amplitude is smooth enough
            return oscillatory_method::asymptotic;
        }
        
        // Default to Levin for general nonlinear phase
        return oscillatory_method::levin;
    }
    
private:
    void validate() const {
        if (!std::isfinite(omega) || omega <= 0) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Frequency omega must be positive and finite"));
        }
        if (!std::isfinite(a) || !std::isfinite(b)) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Integration bounds must be finite"));
        }
        if (a >= b) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Lower bound must be less than upper bound"));
        }
    }
};

/// \brief Filon moments for oscillatory integration
template <typename Real>
struct filon_moments {
    Real alpha;  // Coefficient for endpoint terms
    Real beta;   // Coefficient for even interior points
    Real gamma;  // Coefficient for odd interior points
    
    // Compute Filon moments for cos/sin oscillators
    static filon_moments compute(Real theta, oscillator_type type) {
        filon_moments m;
        
        if (std::abs(theta) < Real(0.1)) {
            // Use Taylor series for small theta
            return compute_taylor(theta, type);
        }
        
        Real theta2 = theta * theta;
        Real theta3 = theta2 * theta;
        Real sin_theta = std::sin(theta);
        Real cos_theta = std::cos(theta);
        
        if (type == oscillator_type::cosine) {
            m.alpha = (theta * sin_theta + 2 * (1 - cos_theta)) / theta3;
            m.beta = 2 * (theta * (1 + cos_theta) - 2 * sin_theta) / theta3;
            m.gamma = 4 * (sin_theta - theta * cos_theta) / theta3;
        } else if (type == oscillator_type::sine) {
            m.alpha = (theta * cos_theta + 2 * sin_theta - 2 * theta) / theta3;
            m.beta = 2 * (theta * (1 - cos_theta) - 2 * sin_theta + 2 * theta) / theta3;
            m.gamma = 4 * (theta - sin_theta) / theta3;
        } else {
            // Default to cosine for other types
            m = compute(theta, oscillator_type::cosine);
        }
        
        return m;
    }
    
private:
    static filon_moments compute_taylor(Real theta, oscillator_type type) {
        filon_moments m;
        Real theta2 = theta * theta;
        Real theta4 = theta2 * theta2;
        
        if (type == oscillator_type::cosine) {
            // Taylor series for alpha
            m.alpha = Real(2) / Real(3) - theta2 / Real(15) + theta4 * Real(2) / Real(315);
            // Taylor series for beta  
            m.beta = Real(4) / Real(3) + theta2 * Real(2) / Real(15) + theta4 * Real(2) / Real(105);
            // Taylor series for gamma
            m.gamma = Real(4) / Real(3) - theta2 * Real(2) / Real(15) + theta4 / Real(210);
        } else if (type == oscillator_type::sine) {
            // Taylor series for sine case
            Real theta3 = theta2 * theta;
            Real theta5 = theta4 * theta;
            m.alpha = theta / Real(3) - theta3 / Real(30) + theta5 / Real(840);
            m.beta = Real(2) * theta / Real(3) + theta3 / Real(10) - theta5 / Real(280);
            m.gamma = Real(4) * theta / Real(3) - theta3 * Real(2) / Real(15) + theta5 / Real(210);
        } else {
            m.alpha = m.beta = m.gamma = Real(0);
        }
        
        return m;
    }
};

/// \brief Helper to detect singularities in phase function
template <typename Real>
struct phase_analyzer {
    using phase_function = std::function<Real(Real)>;
    
    static bool has_stationary_point(
        const phase_function& g,
        Real a, Real b,
        Real& stationary_x)
    {
        // Use simple bisection to find g'(x) = 0
        // This is a simplified version - production code would use more robust methods
        
        const int max_iter = 50;
        const Real tol = std::numeric_limits<Real>::epsilon() * Real(100);
        const Real h = std::sqrt(std::numeric_limits<Real>::epsilon());
        
        auto derivative = [&](Real x) {
            // Numerical derivative
            return (g(x + h) - g(x - h)) / (Real(2) * h);
        };
        
        Real d_a = derivative(a);
        Real d_b = derivative(b);
        
        // Check if derivative changes sign
        if (d_a * d_b > 0) {
            return false; // No stationary point
        }
        
        // Bisection
        Real left = a, right = b;
        for (int iter = 0; iter < max_iter; ++iter) {
            Real mid = (left + right) / Real(2);
            Real d_mid = derivative(mid);
            
            if (std::abs(d_mid) < tol) {
                stationary_x = mid;
                return true;
            }
            
            if (d_a * d_mid < 0) {
                right = mid;
                d_b = d_mid;
            } else {
                left = mid;
                d_a = d_mid;
            }
        }
        
        stationary_x = (left + right) / Real(2);
        return true;
    }
    
    // Estimate the rate of oscillation
    static Real oscillation_rate(
        const phase_function& g,
        Real omega, Real a, Real b)
    {
        // Average number of oscillations per unit length
        Real phase_change = std::abs(omega * (g(b) - g(a)));
        Real interval_length = b - a;
        return phase_change / (interval_length * boost::math::constants::two_pi<Real>());
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_OSCILLATORY_BASE_HPP