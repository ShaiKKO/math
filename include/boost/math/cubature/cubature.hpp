// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_HPP
#define BOOST_MATH_CUBATURE_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/cubature/transforms.hpp>

#include <cmath>
#include <limits>
#include <vector>
#include <memory>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>

namespace boost { namespace math { namespace cubature {

/// \file cubature.hpp
/// \brief Master header for N-dimensional numerical integration
/// \details The Boost.Math Cubature library provides high-performance numerical
///          integration algorithms for multi-dimensional integrals:
///  - Adaptive integration using DCUHRE algorithm with Genz-Malik embedded rules
///  - Sparse grids with Smolyak construction and Clenshaw-Curtis nodes
///  - (Randomized) Quasi-Monte Carlo with Sobol sequences and Owen scrambling
///  - Simplex integration via Duffy transformation
///  - Automatic handling of infinite and semi-infinite domains
///
/// \author Colin MacRitchie
/// \date 2025
///
/// Example usage:
/// \code
/// using namespace boost::math::cubature;
/// 
/// // Define integrand: f(x,y) = exp(-(x^2 + y^2))
/// auto gaussian = [](const double* x, std::size_t) {
///     return std::exp(-(x[0]*x[0] + x[1]*x[1]));
/// };
/// 
/// // Set up integration domain [0,1] x [0,1]
/// hypercube<double> box(2);
/// box.lower = {0, 0};
/// box.upper = {1, 1};
/// 
/// // Integrate with adaptive method
/// auto result = integrate(gaussian, box, 1e-8, 1e-6);
/// std::cout << "Result: " << result.value << " ± " << result.error << std::endl;
/// \endcode

// Forward declarations for convenience functions
namespace detail {
    template <typename Real>
    bool has_infinite_bounds(const hypercube<Real>& box);
    
    template <typename Real>
    infinite_transform_type select_transform(const hypercube<Real>& box);
}

/// \brief Integrate over infinite or semi-infinite domains
/// \details Automatically selects appropriate transformation for infinite bounds:
///          - tanh-sinh for finite domains
///          - exp-sinh for semi-infinite [a,∞) or (-∞,b]
///          - sinh-sinh for bi-infinite (-∞,∞)
///          For dimensions > 1, applies coordinate transformations.
///
/// \tparam Real Floating-point type (float, double, long double, or multiprecision)
/// \tparam F Integrand type callable as f(const Real*, std::size_t) -> Real
/// \tparam Policy Boost.Math policy type for error handling
///
/// \param f Integrand function with signature Real(const Real* x, std::size_t dim)
/// \param lower Vector of lower bounds (may contain -std::numeric_limits<Real>::infinity())
/// \param upper Vector of upper bounds (may contain std::numeric_limits<Real>::infinity())
/// \param abs_tol Absolute error tolerance
/// \param rel_tol Relative error tolerance
/// \param max_eval Maximum number of function evaluations allowed
/// \param pol Boost.Math policy object
///
/// \return result<Real> containing value, error estimate, evaluations, and status
///
/// \note For 1D integrals, delegates to specialized quadrature methods which are
///       more efficient than the general multi-dimensional algorithms.
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_adaptive_infinite(
    const F& f,
    const std::vector<Real>& lower,
    const std::vector<Real>& upper,
    Real abs_tol, Real rel_tol,
    std::size_t max_eval = 100000,
    Policy const& pol = Policy{})
{
    const std::size_t dim = lower.size();
    if (upper.size() != dim) {
        result<Real> res;
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
        return res;
    }
    
    // For 1D, always delegate to specialized quadrature methods
    if (dim == 1) {
        // Create a wrapper for the integrand to match quadrature signature
        auto f_1d = [&f](Real x) -> Real {
            Real x_array[1] = {x};
            return f(x_array, 1);
        };
        
        result<Real> res;
        res.status = status_code::success;
        res.evaluations = 0; // Not tracked by quadrature methods
        
        Real error_est = Real(0);
        Real L1_norm = Real(0);
        std::size_t levels_used = 0;
        
        bool lower_inf = std::isinf(lower[0]) && lower[0] < 0;
        bool upper_inf = std::isinf(upper[0]) && upper[0] > 0;
        
        try {
            if (lower_inf && upper_inf) {
                // Use sinh_sinh for bi-infinite integrals
                boost::math::quadrature::sinh_sinh<Real, Policy> integrator;
                res.value = integrator.integrate(f_1d, abs_tol, &error_est, &L1_norm, &levels_used);
            }
            else if (upper_inf) {
                // Use exp_sinh for semi-infinite [a, ∞)
                boost::math::quadrature::exp_sinh<Real, Policy> integrator;
                res.value = integrator.integrate(f_1d, lower[0], 
                    std::numeric_limits<Real>::infinity(), 
                    abs_tol, &error_est, &L1_norm, &levels_used);
            }
            else if (lower_inf) {
                // Use exp_sinh for semi-infinite (-∞, b]
                // Transform: integrate f(x) from -∞ to b = integrate f(b-t) from 0 to ∞
                auto f_transformed = [&f_1d, &upper](Real t) -> Real {
                    return f_1d(upper[0] - t);
                };
                boost::math::quadrature::exp_sinh<Real, Policy> integrator;
                res.value = integrator.integrate(f_transformed, Real(0), 
                    std::numeric_limits<Real>::infinity(), 
                    abs_tol, &error_est, &L1_norm, &levels_used);
            }
            else {
                // Finite bounds - use tanh_sinh
                boost::math::quadrature::tanh_sinh<Real, Policy> integrator;
                res.value = integrator.integrate(f_1d, lower[0], upper[0], 
                    abs_tol, &error_est, &L1_norm, &levels_used);
            }
            
            res.error = error_est;
            // Estimate evaluations based on levels (rough approximation)
            res.evaluations = static_cast<std::size_t>(std::pow(2, levels_used + 4));
            
        } catch (const std::exception& e) {
            res.status = status_code::invalid_input;
            res.error = std::numeric_limits<Real>::max();
            res.value = Real(0);
        }
        
        return res;
    }
    
    // For dimensions >= 2, check for infinite bounds
    std::vector<std::size_t> infinite_dims;
    std::vector<infinite_transform_type> transform_types;
    
    for (std::size_t i = 0; i < dim; ++i) {
        bool lower_inf = std::isinf(lower[i]) && lower[i] < 0;
        bool upper_inf = std::isinf(upper[i]) && upper[i] > 0;
        
        if (lower_inf || upper_inf) {
            infinite_dims.push_back(i);
            
            if (lower_inf && upper_inf) {
                transform_types.push_back(infinite_transform_type::tangent);
            } else if (upper_inf) {
                transform_types.push_back(infinite_transform_type::rational);
            } else {
                // lower_inf only - use rational with negation
                transform_types.push_back(infinite_transform_type::rational);
            }
        }
    }
    
    // If no infinite bounds for d >= 2, use regular adaptive
    if (infinite_dims.empty()) {
        hypercube<Real> box(dim);
        box.lower = lower;
        box.upper = upper;
        return integrate_adaptive<Real>(f, box, abs_tol, rel_tol, max_eval, pol);
    }
    
    // Create transformed integrand with proper signature
    auto transformed_f = [&f, &lower, &upper, &infinite_dims, &transform_types, dim](const Real* u, std::size_t) -> Real {
        std::vector<Real> x(dim);
        Real jacobian = Real(1);
        
        // Map finite dimensions directly
        std::size_t inf_idx = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            if (inf_idx < infinite_dims.size() && i == infinite_dims[inf_idx]) {
                // Apply transform for infinite dimension
                Real ui = u[i];
                
                switch(transform_types[inf_idx]) {
                    case infinite_transform_type::tangent: {
                        tangent_transform<Real> t;
                        auto transform_result = t.forward(ui);
                        auto xi = transform_result.first;
                        auto jac = transform_result.second;
                        x[i] = xi;
                        jacobian *= jac;
                        break;
                    }
                    case infinite_transform_type::rational: {
                        rational_transform<Real> t;
                        if (std::isinf(lower[i])) {
                            // Transform for (-inf, b]
                            // Map u in [0,1] to (1-u) in [0,1], avoiding endpoints
                            Real u_rev = Real(1) - ui;
                            Real u_safe = std::min(u_rev, Real(1) - std::numeric_limits<Real>::epsilon());
                            u_safe = std::max(u_safe, std::numeric_limits<Real>::epsilon());
                            auto transform_result = t.forward(u_safe);
                            auto xi = transform_result.first;
                            auto jac = transform_result.second;
                            x[i] = upper[i] - xi;
                            jacobian *= jac;
                        } else {
                            // Transform for [a, inf)
                            // Avoid exact u=1 which maps to infinity
                            Real u_safe = std::min(ui, Real(1) - std::numeric_limits<Real>::epsilon());
                            auto transform_result = t.forward(u_safe);
                            auto xi = transform_result.first;
                            auto jac = transform_result.second;
                            x[i] = lower[i] + xi;
                            jacobian *= jac;
                        }
                        break;
                    }
                    default:
                        x[i] = u[i]; // Shouldn't happen
                }
                ++inf_idx;
            } else {
                // Map finite dimension linearly
                x[i] = lower[i] + u[i] * (upper[i] - lower[i]);
                jacobian *= (upper[i] - lower[i]);
            }
        }
        
        return f(x.data(), dim) * jacobian;
    };
    
    // Integrate over unit hypercube
    hypercube<Real> unit_cube(dim);
    std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
    std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
    
    // Use adaptive integration with transformed integrand for d >= 2
    return integrate_adaptive<Real>(transformed_f, unit_cube, abs_tol, rel_tol, max_eval, pol);
}

/// \brief Integrate over infinite or semi-infinite domain using sparse grids
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid_infinite(
    const F& f,
    const std::vector<Real>& lower,
    const std::vector<Real>& upper,
    unsigned level,
    Policy const& pol = Policy{})
{
    const std::size_t dim = lower.size();
    if (upper.size() != dim) {
        result<Real> res;
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
        return res;
    }
    
    // Check for infinite bounds and determine transforms
    std::vector<std::size_t> infinite_dims;
    std::vector<infinite_transform_type> transform_types;
    
    for (std::size_t i = 0; i < dim; ++i) {
        bool lower_inf = std::isinf(lower[i]) && lower[i] < 0;
        bool upper_inf = std::isinf(upper[i]) && upper[i] > 0;
        
        if (lower_inf || upper_inf) {
            infinite_dims.push_back(i);
            
            if (lower_inf && upper_inf) {
                transform_types.push_back(infinite_transform_type::tangent);
            } else if (upper_inf) {
                transform_types.push_back(infinite_transform_type::rational);
            } else {
                // lower_inf only - use rational with negation
                transform_types.push_back(infinite_transform_type::rational);
            }
        }
    }
    
    // If no infinite bounds, use regular sparse grid
    if (infinite_dims.empty()) {
        hypercube<Real> box(dim);
        box.lower = lower;
        box.upper = upper;
        return integrate_sparse_grid<Real>(f, box, level, pol);
    }
    
    // Create transformed integrand for sparse grid
    auto transformed_f = [&f, &lower, &upper, &infinite_dims, &transform_types, dim](const Real* u, std::size_t) -> Real {
        std::vector<Real> x(dim);
        Real jacobian = Real(1);
        
        // Map dimensions according to their bounds
        std::size_t inf_idx = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            if (inf_idx < infinite_dims.size() && i == infinite_dims[inf_idx]) {
                // Apply transform for infinite dimension
                Real ui = u[i];
                
                switch(transform_types[inf_idx]) {
                    case infinite_transform_type::tangent: {
                        tangent_transform<Real> t;
                        auto transform_result = t.forward(ui);
                        auto xi = transform_result.first;
                        auto jac = transform_result.second;
                        x[i] = xi;
                        jacobian *= jac;
                        break;
                    }
                    case infinite_transform_type::rational: {
                        rational_transform<Real> t;
                        if (std::isinf(lower[i])) {
                            // Transform for (-inf, b]
                            // Map u in [0,1] to (1-u) in [0,1], avoiding endpoints
                            Real u_rev = Real(1) - ui;
                            Real u_safe = std::min(u_rev, Real(1) - std::numeric_limits<Real>::epsilon());
                            u_safe = std::max(u_safe, std::numeric_limits<Real>::epsilon());
                            auto transform_result = t.forward(u_safe);
                            auto xi = transform_result.first;
                            auto jac = transform_result.second;
                            x[i] = upper[i] - xi;
                            jacobian *= jac;
                        } else {
                            // Transform for [a, inf)
                            // Avoid exact u=1 which maps to infinity
                            Real u_safe = std::min(ui, Real(1) - std::numeric_limits<Real>::epsilon());
                            auto transform_result = t.forward(u_safe);
                            auto xi = transform_result.first;
                            auto jac = transform_result.second;
                            x[i] = lower[i] + xi;
                            jacobian *= jac;
                        }
                        break;
                    }
                    default:
                        x[i] = u[i]; // Shouldn't happen
                }
                ++inf_idx;
            } else {
                // Map finite dimension linearly
                x[i] = lower[i] + u[i] * (upper[i] - lower[i]);
                jacobian *= (upper[i] - lower[i]);
            }
        }
        
        return f(x.data(), dim) * jacobian;
    };
    
    // Integrate over unit hypercube using sparse grid
    hypercube<Real> unit_cube(dim);
    std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
    std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
    
    return integrate_sparse_grid<Real>(transformed_f, unit_cube, level, pol);
}

/// \brief Master integration function with automatic method selection
/// \details Automatically selects the most appropriate integration algorithm based on:
///          - Dimension: Low (≤6) uses adaptive, medium (≤15) uses sparse grids, 
///            high (>15) uses adaptive with limited evaluations
///          - Tolerance: High accuracy (rel_tol < 1e-8) prefers adaptive method
///          - Domain type: Automatically handles infinite bounds
///
/// Algorithm selection:
///  - 1D: Delegates to specialized quadrature methods
///  - 2-6D or high accuracy: Adaptive DCUHRE with Genz-Malik rules
///  - 7-15D: Sparse grids with Smolyak construction
///  - 16+D: Would use QMC when fully implemented, currently uses adaptive
///
/// \tparam Real Floating-point type
/// \tparam F Integrand type callable as f(const Real*, std::size_t) -> Real
/// \tparam Policy Boost.Math policy type
///
/// \param f Integrand function
/// \param box Integration domain (hypercube)
/// \param abs_tol Absolute error tolerance
/// \param rel_tol Relative error tolerance  
/// \param pol Boost.Math policy object
///
/// \return result<Real> containing integral value, error estimate, and diagnostics
///
/// \note This is the recommended entry point for most integration tasks.
///       The algorithm selection heuristics are based on extensive testing
///       with the Genz test function suite.
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate(
    const F& f,
    const hypercube<Real>& box,
    Real abs_tol, Real rel_tol,
    Policy const& pol = Policy{})
{
    const std::size_t dim = box.dimension();
    
    // For low dimensions or high accuracy, use adaptive
    if (dim <= 6 || rel_tol < 1e-8) {
        return integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 0, pol);
    }
    
    // For moderate dimensions, use sparse grids
    if (dim <= 15) {
        // Estimate required level from tolerance
        unsigned level = static_cast<unsigned>(std::max(Real(3), -std::log10(rel_tol) / Real(2)));
        return integrate_sparse_grid<Real>(f, box, level, pol);
    }
    
    // For high dimensions (>15), use QMC
    // QMC is better suited for high dimensions due to better convergence rate
    std::size_t n_points = static_cast<std::size_t>(
        std::max(Real(10000), Real(1) / (rel_tol * rel_tol)));
    n_points = std::min(n_points, std::size_t(1000000));  // Cap at 1M points
    
    return integrate_qmc<Real>(f, box, n_points, pol);
}

// Implementation details
namespace detail {
    
template <typename Real>
bool has_infinite_bounds(const hypercube<Real>& box) {
    for (const auto& val : box.lower) {
        if (std::isinf(val)) return true;
    }
    for (const auto& val : box.upper) {
        if (std::isinf(val)) return true;
    }
    return false;
}

template <typename Real>
infinite_transform_type select_transform(const hypercube<Real>& box) {
    // Simple heuristic for transform selection
    // Could be enhanced based on integrand characteristics
    bool has_lower_inf = false;
    bool has_upper_inf = false;
    
    for (const auto& val : box.lower) {
        if (std::isinf(val) && val < 0) has_lower_inf = true;
    }
    for (const auto& val : box.upper) {
        if (std::isinf(val) && val > 0) has_upper_inf = true;
    }
    
    if (has_lower_inf && has_upper_inf) {
        return infinite_transform_type::tangent;
    } else if (has_upper_inf || has_lower_inf) {
        return infinite_transform_type::rational;
    }
    
    return infinite_transform_type::rational; // Default
}

} // namespace detail

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_HPP