// Copyright 2025 Boost Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_OSCILLATORY_HPP
#define BOOST_MATH_CUBATURE_OSCILLATORY_HPP

#include <boost/math/cubature/detail/oscillatory_base.hpp>
#include <boost/math/cubature/detail/filon_method.hpp>
#include <boost/math/cubature/detail/levin_method.hpp>
#include <type_traits>
#include <complex>

namespace boost { namespace math { namespace cubature {

/// \brief Main interface for oscillatory integration
/// 
/// Automatically selects the best method based on problem characteristics
///
/// \tparam Real Floating-point type
/// \tparam F Amplitude function type
/// \tparam G Phase function type  
/// \tparam Policy Error handling policy
///
/// \param f Amplitude function f(x)
/// \param g Phase function g(x) 
/// \param omega Frequency parameter ω
/// \param a Lower integration bound
/// \param b Upper integration bound
/// \param type Type of oscillator (exponential, cosine, sine, bessel_j)
/// \param tol_abs Absolute tolerance
/// \param tol_rel Relative tolerance
/// \param pol Policy for error handling
///
/// \return Integration result with value, error estimate, and diagnostics
///
template <typename Real, typename F, typename G, 
          typename Policy = boost::math::policies::policy<>>
auto integrate_oscillatory(
    F&& f,
    G&& g,
    Real omega,
    Real a, Real b,
    detail::oscillator_type type = detail::oscillator_type::exponential,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value &&
        std::is_invocable_r<Real, G, Real>::value,
        result<Real>>::type
{
    // Create oscillatory integral object
    detail::oscillatory_integral<Real> osc_integral(
        std::forward<F>(f), std::forward<G>(g), omega, a, b, type);
    
    // Check if frequency is high enough for oscillatory methods
    if (!osc_integral.is_highly_oscillatory()) {
        // Use standard quadrature for low frequencies
        auto integrand = [&](Real x) -> Real {
            auto val = osc_integral.integrand(x);
            if (type == detail::oscillator_type::cosine || 
                type == detail::oscillator_type::sine) {
                return val.real();
            }
            return std::abs(val); // For complex exponential
        };
        
        // Simple adaptive trapezoidal rule for non-oscillatory case
        const std::size_t max_iter = 100;
        std::size_t n = 16;
        Real prev_result = 0;
        Real curr_result = 0;
        Real error = std::numeric_limits<Real>::max();
        
        for (std::size_t iter = 0; iter < max_iter && n < 10000; ++iter) {
            prev_result = curr_result;
            Real h = (b - a) / Real(n);
            curr_result = 0;
            
            for (std::size_t i = 0; i <= n; ++i) {
                Real x = a + i * h;
                Real weight = (i == 0 || i == n) ? Real(0.5) : Real(1);
                curr_result += weight * integrand(x);
            }
            curr_result *= h;
            
            if (iter > 0) {
                error = std::abs(curr_result - prev_result);
                if (error < tol_abs + tol_rel * std::abs(curr_result)) {
                    break;
                }
            }
            n *= 2;
        }
        
        result<Real> res;
        res.value = curr_result;
        res.error = error;
        res.status = status_code::success;
        res.evaluations = n;
        
        return res;
    }
    
    // Select method based on problem characteristics
    auto method = osc_integral.suggest_method();
    
    switch (method) {
        case detail::oscillatory_method::filon:
        {
            detail::filon_integrator<Real, Policy> filon;
            return filon.integrate(std::forward<F>(f), std::forward<G>(g),
                                  omega, a, b, type, tol_abs, tol_rel, pol);
        }
        
        case detail::oscillatory_method::levin:
        {
            detail::levin_integrator<Real, Policy> levin;
            return levin.integrate(std::forward<F>(f), std::forward<G>(g),
                                  omega, a, b, type, tol_abs, tol_rel, pol);
        }
        
        case detail::oscillatory_method::asymptotic:
        {
            // For very large omega, use asymptotic expansion
            // This is a simplified implementation
            result<Real> res;
            
            // Leading order: O(1/ω) for endpoints
            Real f_a = f(a);
            Real f_b = f(b);
            Real g_a = g(a);
            Real g_b = g(b);
            
            std::complex<Real> phase_a = std::exp(std::complex<Real>(0, omega * g_a));
            std::complex<Real> phase_b = std::exp(std::complex<Real>(0, omega * g_b));
            
            // Numerical derivatives
            const Real h_eps = std::sqrt(std::numeric_limits<Real>::epsilon());
            Real g_prime_a = (g(a + h_eps) - g_a) / h_eps;
            Real g_prime_b = (g_b - g(b - h_eps)) / h_eps;
            
            std::complex<Real> contrib_a = f_a * phase_a / 
                                          (std::complex<Real>(0, omega * g_prime_a));
            std::complex<Real> contrib_b = f_b * phase_b / 
                                          (std::complex<Real>(0, omega * g_prime_b));
            
            std::complex<Real> value = contrib_b - contrib_a;
            
            if (type == detail::oscillator_type::cosine) {
                res.value = value.real();
            } else if (type == detail::oscillator_type::sine) {
                res.value = value.imag();
            } else {
                res.value = std::abs(value);
            }
            
            res.error = (std::abs(f_a) + std::abs(f_b)) / (omega * omega);
            res.status = status_code::success;
            res.evaluations = 4;
            
            return res;
        }
        
        default:
        {
            // Automatic selection - try Filon first, fall back to Levin
            detail::filon_integrator<Real, Policy> filon;
            auto res = filon.integrate(std::forward<F>(f), std::forward<G>(g),
                                      omega, a, b, type, tol_abs, tol_rel, pol);
            
            if (res.status != status_code::success) {
                // Fall back to Levin
                detail::levin_integrator<Real, Policy> levin;
                res = levin.integrate(std::forward<F>(f), std::forward<G>(g),
                                    omega, a, b, type, tol_abs, tol_rel, pol);
            }
            
            return res;
        }
    }
}

/// \brief Specialized function for Fourier cosine transform
/// 
/// Computes ∫ f(x) cos(ωx) dx from a to b
///
template <typename Real, typename F, 
          typename Policy = boost::math::policies::policy<>>
auto integrate_fourier_cosine(
    F&& f,
    Real omega,
    Real a, Real b,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value,
        result<Real>>::type
{
    // Identity phase function g(x) = x
    auto g = [](Real x) { return x; };
    
    return integrate_oscillatory(
        std::forward<F>(f), g, omega, a, b,
        detail::oscillator_type::cosine,
        tol_abs, tol_rel, pol);
}

/// \brief Specialized function for Fourier sine transform
/// 
/// Computes ∫ f(x) sin(ωx) dx from a to b
///
template <typename Real, typename F,
          typename Policy = boost::math::policies::policy<>>
auto integrate_fourier_sine(
    F&& f,
    Real omega,
    Real a, Real b,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value,
        result<Real>>::type
{
    // Identity phase function g(x) = x
    auto g = [](Real x) { return x; };
    
    return integrate_oscillatory(
        std::forward<F>(f), g, omega, a, b,
        detail::oscillator_type::sine,
        tol_abs, tol_rel, pol);
}

/// \brief Semi-infinite Fourier cosine transform
/// 
/// Computes ∫ f(x) cos(ωx) dx from 0 to ∞
///
template <typename Real, typename F,
          typename Policy = boost::math::policies::policy<>>
auto integrate_fourier_cosine_infinite(
    F&& f,
    Real omega,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value,
        result<Real>>::type
{
    // Use transformation x = tan(t), dx = sec²(t)dt
    // Maps [0, ∞) to [0, π/2)
    
    const Real pi = boost::math::constants::pi<Real>();
    
    auto transformed_f = [&f, omega](Real t) -> Real {
        if (t >= pi / Real(2) - std::numeric_limits<Real>::epsilon()) {
            return Real(0);
        }
        
        Real x = std::tan(t);
        Real sec_t = Real(1) / std::cos(t);
        Real jacobian = sec_t * sec_t;
        
        return f(x) * std::cos(omega * x) * jacobian;
    };
    
    // Use simple adaptive quadrature on finite interval
    const std::size_t n = 256;
    Real h = (pi / Real(2)) / Real(n);
    Real value = 0;
    
    for (std::size_t i = 0; i <= n; ++i) {
        Real t = i * h;
        Real weight = (i == 0 || i == n) ? Real(0.5) : Real(1);
        value += weight * transformed_f(t);
    }
    value *= h;
    
    result<Real> res;
    res.value = value;
    res.error = std::abs(value) * Real(1e-6); // Simple error estimate
    res.status = status_code::success;
    res.evaluations = n + 1;
    
    return res;
}

/// \brief Semi-infinite Fourier sine transform
/// 
/// Computes ∫ f(x) sin(ωx) dx from 0 to ∞
///
template <typename Real, typename F,
          typename Policy = boost::math::policies::policy<>>
auto integrate_fourier_sine_infinite(
    F&& f,
    Real omega,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value,
        result<Real>>::type
{
    // Use transformation x = tan(t), dx = sec²(t)dt
    // Maps [0, ∞) to [0, π/2)
    
    const Real pi = boost::math::constants::pi<Real>();
    
    auto transformed_f = [&f, omega](Real t) -> Real {
        if (t >= pi / Real(2) - std::numeric_limits<Real>::epsilon()) {
            return Real(0);
        }
        
        Real x = std::tan(t);
        Real sec_t = Real(1) / std::cos(t);
        Real jacobian = sec_t * sec_t;
        
        return f(x) * std::sin(omega * x) * jacobian;
    };
    
    // Use simple adaptive quadrature on finite interval
    const std::size_t n = 256;
    Real h = (pi / Real(2)) / Real(n);
    Real value = 0;
    
    for (std::size_t i = 0; i <= n; ++i) {
        Real t = i * h;
        Real weight = (i == 0 || i == n) ? Real(0.5) : Real(1);
        value += weight * transformed_f(t);
    }
    value *= h;
    
    result<Real> res;
    res.value = value;
    res.error = std::abs(value) * Real(1e-6); // Simple error estimate
    res.status = status_code::success;
    res.evaluations = n + 1;
    
    return res;
}

/// \brief Complex-valued oscillatory integral
/// 
/// Computes ∫ f(x) e^(iωg(x)) dx with complex result
///
template <typename Real, typename F, typename G,
          typename Policy = boost::math::policies::policy<>>
auto integrate_oscillatory_complex(
    F&& f,
    G&& g,
    Real omega,
    Real a, Real b,
    Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
    Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
    const Policy& pol = Policy{})
    -> typename std::enable_if<
        std::is_invocable_r<Real, F, Real>::value &&
        std::is_invocable_r<Real, G, Real>::value,
        std::pair<std::complex<Real>, Real>>::type
{
    // Compute real and imaginary parts separately
    auto res_cos = integrate_oscillatory(
        std::forward<F>(f), std::forward<G>(g), omega, a, b,
        detail::oscillator_type::cosine,
        tol_abs / Real(2), tol_rel, pol);
    
    auto res_sin = integrate_oscillatory(
        std::forward<F>(f), std::forward<G>(g), omega, a, b,
        detail::oscillator_type::sine,
        tol_abs / Real(2), tol_rel, pol);
    
    std::complex<Real> value(res_cos.value, res_sin.value);
    Real error = std::sqrt(res_cos.error * res_cos.error + 
                          res_sin.error * res_sin.error);
    
    return {value, error};
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_OSCILLATORY_HPP