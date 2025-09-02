// Copyright 2025 Boost Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_FILON_METHOD_HPP
#define BOOST_MATH_CUBATURE_DETAIL_FILON_METHOD_HPP

#include <boost/math/cubature/detail/oscillatory_base.hpp>
#include <boost/math/tools/roots.hpp>
#include <vector>
#include <algorithm>
#include <numeric>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Filon quadrature for oscillatory integrals
/// 
/// Implements classical Filon method and extended Filon-Clenshaw-Curtis
/// for integrals of the form ∫ f(x) cos(ωg(x)) dx or ∫ f(x) sin(ωg(x)) dx
///
template <typename Real, typename Policy = boost::math::policies::policy<>>
class filon_integrator {
public:
    using result_type = boost::math::cubature::result<Real>;
    using complex_type = std::complex<Real>;
    
    filon_integrator(std::size_t max_subdivisions = 100)
        : max_subdivisions_(max_subdivisions)
        , nodes_per_interval_(64) // Default to 64 nodes per subinterval
    {}
    
    /// \brief Integrate using Filon method
    template <typename F, typename G>
    result_type integrate(
        F&& f,           // Amplitude function
        G&& g,           // Phase function  
        Real omega,      // Frequency
        Real a, Real b,  // Integration bounds
        oscillator_type type,
        Real tol_abs = std::numeric_limits<Real>::epsilon() * Real(100),
        Real tol_rel = std::sqrt(std::numeric_limits<Real>::epsilon()),
        const Policy& pol = Policy{})
    {
        // Input validation
        if (!validate_input(omega, a, b)) {
            return make_error_result(status_code::invalid_input);
        }
        
        // Check if phase is linear (standard Filon case)
        bool is_linear = is_linear_phase(g, a, b);
        
        if (is_linear) {
            // Use classical Filon for linear phase
            return filon_classical(std::forward<F>(f), omega, a, b, type, tol_abs, tol_rel, pol);
        } else {
            // Use extended Filon with phase transformation
            return filon_extended(std::forward<F>(f), std::forward<G>(g), 
                                 omega, a, b, type, tol_abs, tol_rel, pol);
        }
    }
    
    /// \brief Extended Filon-Clenshaw-Curtis method
    template <typename F>
    result_type integrate_clenshaw_curtis(
        F&& f,
        Real omega,
        Real a, Real b,
        oscillator_type type,
        std::size_t n = 128,  // Number of Clenshaw-Curtis nodes
        const Policy& pol = Policy{})
    {
        // Generate Clenshaw-Curtis nodes and weights
        std::vector<Real> nodes(n + 1);
        std::vector<Real> weights(n + 1);
        
        generate_clenshaw_curtis_nodes(nodes, weights, a, b, n);
        
        // Compute Filon moments for each subinterval
        Real h = (b - a) / Real(n);
        Real theta = omega * h;
        
        auto moments = filon_moments<Real>::compute(theta, type);
        
        // Evaluate function at nodes
        std::vector<Real> f_vals(n + 1);
        for (std::size_t i = 0; i <= n; ++i) {
            f_vals[i] = f(nodes[i]);
        }
        
        // Apply Filon formula
        complex_type sum(0, 0);
        
        if (type == oscillator_type::cosine) {
            sum = apply_filon_formula_cos(f_vals, nodes, omega, h, moments);
        } else if (type == oscillator_type::sine) {
            sum = apply_filon_formula_sin(f_vals, nodes, omega, h, moments);
        }
        
        // Error estimation using Richardson extrapolation
        Real error = estimate_error_richardson(std::forward<F>(f), omega, a, b, type, sum.real());
        
        result_type res;
        res.value = sum.real();
        res.error = error;
        res.status = status_code::success;
        res.evaluations = n + 1;
        
        return res;
    }
    
private:
    std::size_t max_subdivisions_;
    std::size_t nodes_per_interval_;
    
    /// \brief Classical Filon method for linear phase g(x) = x
    template <typename F>
    result_type filon_classical(
        F&& f,
        Real omega,
        Real a, Real b,
        oscillator_type type,
        Real tol_abs, Real tol_rel,
        const Policy& pol)
    {
        // Adaptive subdivision strategy
        std::vector<Real> breakpoints = {a, b};
        
        // Add more breakpoints if interval is large
        Real interval_length = b - a;
        Real wavelength = boost::math::constants::two_pi<Real>() / omega;
        
        if (interval_length > Real(10) * wavelength) {
            // Subdivide to have roughly 10 wavelengths per subinterval
            std::size_t n_subs = static_cast<std::size_t>(interval_length / (Real(10) * wavelength));
            n_subs = std::min(n_subs, max_subdivisions_);
            
            for (std::size_t i = 1; i < n_subs; ++i) {
                breakpoints.push_back(a + i * interval_length / Real(n_subs));
            }
            std::sort(breakpoints.begin(), breakpoints.end());
        }
        
        complex_type total_integral(0, 0);
        Real total_error = 0;
        std::size_t total_evals = 0;
        
        // Integrate over each subinterval
        for (std::size_t i = 0; i < breakpoints.size() - 1; ++i) {
            Real sub_a = breakpoints[i];
            Real sub_b = breakpoints[i + 1];
            
            auto [sub_integral, sub_error, sub_evals] = 
                filon_subinterval(std::forward<F>(f), omega, sub_a, sub_b, type);
            
            total_integral += sub_integral;
            total_error += sub_error * sub_error; // Add errors in quadrature
            total_evals += sub_evals;
        }
        
        total_error = std::sqrt(total_error);
        
        result_type res;
        res.value = total_integral.real();
        res.error = total_error;
        res.evaluations = total_evals;
        res.status = (total_error <= tol_abs + tol_rel * std::abs(res.value)) 
                     ? status_code::success : status_code::maxeval_reached;
        
        return res;
    }
    
    /// \brief Extended Filon for nonlinear phase
    template <typename F, typename G>
    result_type filon_extended(
        F&& f, G&& g,
        Real omega,
        Real a, Real b,
        oscillator_type type,
        Real tol_abs, Real tol_rel,
        const Policy& pol)
    {
        // Transform to standard form using change of variables
        // u = g(x), then integral becomes ∫ f(g^{-1}(u))/g'(g^{-1}(u)) cos(ωu) du
        
        // For now, use piecewise Filon with local linearization
        // This is a simplified approach - full implementation would use
        // stationary phase analysis and contour deformation
        
        std::size_t n_intervals = std::min(std::size_t(20), max_subdivisions_);
        Real h = (b - a) / Real(n_intervals);
        
        complex_type total(0, 0);
        Real total_error = 0;
        std::size_t total_evals = 0;
        
        for (std::size_t i = 0; i < n_intervals; ++i) {
            Real x0 = a + i * h;
            Real x1 = x0 + h;
            
            // Local linearization of phase
            Real g0 = g(x0);
            Real g1 = g(x1);
            Real g_slope = (g1 - g0) / h;
            
            // Modified amplitude including Jacobian
            auto f_modified = [&f](Real x) {
                return f(x);  // Simplified - should include Jacobian correction
            };
            
            // Apply Filon to linearized problem
            Real local_omega = omega * g_slope;
            Real local_phase_shift = omega * (g0 - g_slope * x0);
            
            auto [sub_val, sub_err, sub_evals] = 
                filon_subinterval(f_modified, local_omega, x0, x1, type);
            
            // Apply phase shift
            complex_type phase_factor = std::exp(complex_type(0, local_phase_shift));
            total += sub_val * phase_factor;
            total_error += sub_err * sub_err;
            total_evals += sub_evals;
        }
        
        result_type res;
        res.value = total.real();
        res.error = std::sqrt(total_error);
        res.evaluations = total_evals;
        res.status = (res.error <= tol_abs + tol_rel * std::abs(res.value))
                     ? status_code::success : status_code::maxeval_reached;
        
        return res;
    }
    
    /// \brief Integrate over a single subinterval using Filon formula
    template <typename F>
    std::tuple<complex_type, Real, std::size_t> filon_subinterval(
        F&& f,
        Real omega,
        Real a, Real b,
        oscillator_type type)
    {
        const std::size_t n = 64; // Number of points (must be even)
        Real h = (b - a) / Real(n);
        Real theta = omega * h;
        
        // Compute Filon moments
        auto moments = filon_moments<Real>::compute(theta, type);
        
        // Evaluate function at nodes
        std::vector<Real> f_vals(n + 1);
        std::vector<Real> x_vals(n + 1);
        
        for (std::size_t i = 0; i <= n; ++i) {
            x_vals[i] = a + i * h;
            f_vals[i] = f(x_vals[i]);
        }
        
        // Apply Filon formula
        complex_type result(0, 0);
        
        if (type == oscillator_type::cosine) {
            result = apply_filon_formula_cos(f_vals, x_vals, omega, h, moments);
        } else if (type == oscillator_type::sine) {
            result = apply_filon_formula_sin(f_vals, x_vals, omega, h, moments);
        }
        
        // Simple error estimate based on difference with midpoint rule
        Real midpoint_val = f((a + b) / Real(2));
        Real error = std::abs(h * midpoint_val) / (omega * omega);
        
        return {result, error, n + 1};
    }
    
    /// \brief Apply Filon formula for cosine oscillator
    complex_type apply_filon_formula_cos(
        const std::vector<Real>& f_vals,
        const std::vector<Real>& x_vals,
        Real omega,
        Real h,
        const filon_moments<Real>& moments)
    {
        std::size_t n = f_vals.size() - 1;
        Real sum = 0;
        
        // Endpoint contributions
        sum += moments.alpha * (f_vals[0] * std::cos(omega * x_vals[0]) + 
                                f_vals[n] * std::cos(omega * x_vals[n]));
        
        // Even interior points
        for (std::size_t j = 1; j < n; j += 2) {
            if (j + 1 < n) {
                sum += moments.beta * f_vals[j + 1] * std::cos(omega * x_vals[j + 1]);
            }
        }
        
        // Odd interior points  
        for (std::size_t j = 1; j < n; j += 2) {
            sum += moments.gamma * f_vals[j] * std::cos(omega * x_vals[j]);
        }
        
        return complex_type(h * sum, 0);
    }
    
    /// \brief Apply Filon formula for sine oscillator
    complex_type apply_filon_formula_sin(
        const std::vector<Real>& f_vals,
        const std::vector<Real>& x_vals,
        Real omega,
        Real h,
        const filon_moments<Real>& moments)
    {
        std::size_t n = f_vals.size() - 1;
        Real sum = 0;
        
        // Endpoint contributions
        sum += moments.alpha * (f_vals[0] * std::sin(omega * x_vals[0]) + 
                                f_vals[n] * std::sin(omega * x_vals[n]));
        
        // Even interior points
        for (std::size_t j = 1; j < n; j += 2) {
            if (j + 1 < n) {
                sum += moments.beta * f_vals[j + 1] * std::sin(omega * x_vals[j + 1]);
            }
        }
        
        // Odd interior points
        for (std::size_t j = 1; j < n; j += 2) {
            sum += moments.gamma * f_vals[j] * std::sin(omega * x_vals[j]);
        }
        
        return complex_type(h * sum, 0);
    }
    
    /// \brief Generate Clenshaw-Curtis nodes and weights
    void generate_clenshaw_curtis_nodes(
        std::vector<Real>& nodes,
        std::vector<Real>& weights,
        Real a, Real b,
        std::size_t n)
    {
        using std::cos;
        const Real pi = boost::math::constants::pi<Real>();
        
        // Generate nodes in [-1, 1]
        for (std::size_t i = 0; i <= n; ++i) {
            Real theta = pi * Real(i) / Real(n);
            Real x_ref = cos(theta);
            // Transform to [a, b]
            nodes[i] = ((b - a) * x_ref + (b + a)) / Real(2);
        }
        
        // Compute weights (simplified - full implementation would use FFT)
        Real h = (b - a) / Real(2);
        for (std::size_t i = 0; i <= n; ++i) {
            if (i == 0 || i == n) {
                weights[i] = h / Real(n);
            } else {
                weights[i] = Real(2) * h / Real(n);
            }
        }
    }
    
    /// \brief Check if phase function is linear
    template <typename G>
    bool is_linear_phase(G&& g, Real a, Real b) {
        // Sample phase at several points and check linearity
        const std::size_t n_samples = 10;
        Real h = (b - a) / Real(n_samples);
        
        Real g0 = g(a);
        Real g1 = g(a + h);
        Real expected_slope = (g1 - g0) / h;
        
        for (std::size_t i = 2; i <= n_samples; ++i) {
            Real x = a + i * h;
            Real g_expected = g0 + expected_slope * (x - a);
            Real g_actual = g(x);
            
            Real rel_error = std::abs(g_actual - g_expected) / (std::abs(g_expected) + Real(1));
            if (rel_error > Real(0.01)) {
                return false;
            }
        }
        
        return true;
    }
    
    /// \brief Error estimation using Richardson extrapolation
    template <typename F>
    Real estimate_error_richardson(
        F&& f,
        Real omega,
        Real a, Real b,
        oscillator_type type,
        Real coarse_result)
    {
        // Compute with double the number of points
        const std::size_t n_fine = 128;
        
        auto fine_result = integrate_clenshaw_curtis(
            std::forward<F>(f), omega, a, b, type, n_fine).value;
        
        // Richardson extrapolation assuming O(h^2) error
        Real error = std::abs(fine_result - coarse_result) / Real(3);
        
        return error;
    }
    
    bool validate_input(Real omega, Real a, Real b) {
        return std::isfinite(omega) && omega > 0 &&
               std::isfinite(a) && std::isfinite(b) && a < b;
    }
    
    result_type make_error_result(status_code code) {
        result_type res;
        res.value = std::numeric_limits<Real>::quiet_NaN();
        res.error = std::numeric_limits<Real>::infinity();
        res.status = code;
        res.evaluations = 0;
        return res;
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_FILON_METHOD_HPP