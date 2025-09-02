// Copyright 2025 Boost.Math Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_INFINITE_HPP
#define BOOST_MATH_CUBATURE_INFINITE_HPP

#include <boost/math/cubature/detail/domain_transforms.hpp>
#include <boost/math/cubature/detail/algebraic_transforms.hpp>
#include <boost/math/cubature/detail/exponential_transforms.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <memory>
#include <utility>
#include <type_traits>

namespace boost { namespace math { namespace cubature {

/// \brief Result structure for infinite domain integration
template <typename Real>
struct infinite_result {
    Real value;        ///< Integral value
    Real error;        ///< Error estimate
    size_t evaluations; ///< Number of function evaluations
    bool converged;    ///< Convergence flag
    
    infinite_result() 
        : value(0), error(0), evaluations(0), converged(false) {}
    
    infinite_result(Real v, Real e, size_t n, bool c)
        : value(v), error(e), evaluations(n), converged(c) {}
};

/// \brief Transform method selection
enum class transform_method {
    automatic,      ///< Automatic selection based on function analysis
    algebraic,      ///< Algebraic transform for power-law decay
    tanh_sinh,      ///< Tanh-sinh for endpoint singularities
    exp_sinh,       ///< Exp-sinh for exponential decay
    sinh_sinh,      ///< Sinh-sinh for double exponential decay
    imt             ///< IMT for oscillatory integrals
};

namespace detail {

/// \brief Helper to select optimal transform based on domain and method
template <typename Real>
std::unique_ptr<domain_transform<Real>> 
select_transform(domain_type domain, transform_method method, 
                 Real a = 0, Real /*b*/ = 0, Real omega = 0) {
    using namespace detail;
    
    switch (method) {
        case transform_method::algebraic:
            if (domain == domain_type::infinite) {
                return std::make_unique<double_algebraic_transform<Real>>();
            } else {
                return std::make_unique<algebraic_transform<Real>>(a);
            }
            
        case transform_method::tanh_sinh:
            return std::make_unique<tanh_sinh_transform<Real>>();
            
        case transform_method::exp_sinh:
            return std::make_unique<exp_sinh_transform<Real>>(a);
            
        case transform_method::sinh_sinh:
            return std::make_unique<sinh_sinh_transform<Real>>();
            
        case transform_method::imt:
            if (omega <= 0) {
                BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                    "IMT transform requires positive frequency"));
            }
            return std::make_unique<imt_transform<Real>>(omega);
            
        default:
            // Automatic selection based on domain
            if (domain == domain_type::infinite) {
                return std::make_unique<tanh_sinh_transform<Real>>();
            } else {
                return std::make_unique<algebraic_transform<Real>>(a);
            }
    }
}

/// \brief Core integration routine with transform
template <typename Real, typename F, typename Policy>
infinite_result<Real> 
integrate_transformed(F&& f, 
                     std::unique_ptr<domain_transform<Real>>& transform,
                     Real tolerance,
                     size_t max_refinements,
                     const Policy& /*pol*/) {
    using std::abs;
    
    // Get truncation bounds
    auto [t_min, t_max] = transform->truncation_bounds(tolerance);
    
    // Transform the integrand
    auto g = transform->transform_integrand(std::forward<F>(f));
    
    // Use Gauss-Kronrod quadrature on the finite domain
    using namespace boost::math::quadrature;
    
    Real error_estimate = Real(0);
    
    // Perform adaptive integration
    Real result = gauss_kronrod<Real, 15>::integrate(g, t_min, t_max, 
                                                     max_refinements, tolerance,
                                                     &error_estimate);
    
    // Estimate evaluations based on refinement level (approximate)
    size_t evaluations = 15 * (1 << std::min(max_refinements, size_t(10)));
    
    // Check convergence
    bool converged = (error_estimate < tolerance * abs(result)) ||
                    (error_estimate < tolerance);
    
    return infinite_result<Real>(result, error_estimate, evaluations, converged);
}

} // namespace detail

/// \brief Integrate a function over a semi-infinite interval [a, ∞)
/// 
/// \param f Function to integrate
/// \param a Lower bound
/// \param tolerance Desired relative tolerance
/// \param method Transform method to use
/// \param max_refinements Maximum number of adaptive refinements
/// \return Integration result with value, error estimate, and convergence flag
template <typename F, typename Real = double,
          typename Policy = boost::math::policies::policy<>>
infinite_result<Real> 
integrate_ray(F&& f, 
              Real a = 0,
              Real tolerance = std::sqrt(std::numeric_limits<Real>::epsilon()),
              transform_method method = transform_method::automatic,
              size_t max_refinements = 15,
              const Policy& pol = Policy()) {
    
    using namespace detail;
    
    // Analyze function if automatic method selection
    auto transform = select_transform<Real>(domain_type::ray, method, a);
    
    if (method == transform_method::automatic) {
        // Sample function to determine optimal transform
        decay_type decay = transform->analyze_function(
            [&f, a](Real x) { return f(a + x); }, 10);
        
        if (decay == decay_type::exponential) {
            transform = std::make_unique<exp_sinh_transform<Real>>(a);
        } else if (decay == decay_type::power_law) {
            transform = std::make_unique<algebraic_transform<Real>>(a);
        }
        // else keep default
    }
    
    return integrate_transformed<Real>(std::forward<F>(f), transform, 
                                       tolerance, max_refinements, pol);
}

/// \brief Integrate a function over the real line (-∞, ∞)
/// 
/// \param f Function to integrate
/// \param tolerance Desired relative tolerance
/// \param method Transform method to use
/// \param max_refinements Maximum number of adaptive refinements
/// \return Integration result with value, error estimate, and convergence flag
template <typename F, typename Real = double,
          typename Policy = boost::math::policies::policy<>>
infinite_result<Real> 
integrate_real_line(F&& f,
                   Real tolerance = std::sqrt(std::numeric_limits<Real>::epsilon()),
                   transform_method method = transform_method::automatic,
                   size_t max_refinements = 15,
                   const Policy& pol = Policy()) {
    
    using namespace detail;
    
    auto transform = select_transform<Real>(domain_type::infinite, method);
    
    if (method == transform_method::automatic) {
        // Sample function to determine optimal transform
        decay_type decay = transform->analyze_function(f, 10);
        
        if (decay == decay_type::exponential) {
            // Check if it's more like exp(-x²) (use sinh-sinh) or exp(-|x|) (use tanh-sinh)
            Real f0 = abs(f(Real(0)));
            Real f1 = abs(f(Real(1)));
            Real f2 = abs(f(Real(2)));
            
            if (f0 > std::numeric_limits<Real>::epsilon()) {
                Real ratio1 = f1 / f0;
                Real ratio2 = f2 / f1;
                
                // exp(-x²) decays faster than exp(-x)
                if (ratio2 < ratio1 * ratio1) {
                    transform = std::make_unique<sinh_sinh_transform<Real>>();
                } else {
                    transform = std::make_unique<tanh_sinh_transform<Real>>();
                }
            }
        } else if (decay == decay_type::power_law) {
            transform = std::make_unique<double_algebraic_transform<Real>>();
        }
        // else keep default (tanh-sinh)
    }
    
    return integrate_transformed<Real>(std::forward<F>(f), transform,
                                       tolerance, max_refinements, pol);
}

/// \brief General infinite domain integration
/// 
/// Automatically determines the domain type and selects appropriate transform
/// 
/// \param f Function to integrate
/// \param lower Lower bound (can be -∞)
/// \param upper Upper bound (can be +∞)
/// \param tolerance Desired relative tolerance
/// \param method Transform method to use
/// \param max_refinements Maximum number of adaptive refinements
/// \return Integration result
template <typename F, typename Real = double,
          typename Policy = boost::math::policies::policy<>>
infinite_result<Real> 
integrate_infinite(F&& f,
                  Real lower = -std::numeric_limits<Real>::infinity(),
                  Real upper = std::numeric_limits<Real>::infinity(),
                  Real tolerance = std::sqrt(std::numeric_limits<Real>::epsilon()),
                  transform_method method = transform_method::automatic,
                  size_t max_refinements = 15,
                  const Policy& pol = Policy()) {
    
    using std::isinf;
    
    bool lower_inf = isinf(lower);
    bool upper_inf = isinf(upper);
    
    if (lower_inf && upper_inf) {
        // Both infinite: (-∞, ∞)
        return integrate_real_line(std::forward<F>(f), tolerance, method, 
                                  max_refinements, pol);
    } else if (upper_inf) {
        // Semi-infinite: [a, ∞)
        return integrate_ray(std::forward<F>(f), lower, tolerance, method,
                           max_refinements, pol);
    } else if (lower_inf) {
        // Semi-infinite: (-∞, b]
        // Transform to [0, ∞) by substitution x → -x
        auto g = [f, upper](Real x) { return f(upper - x); };
        return integrate_ray(g, Real(0), tolerance, method,
                           max_refinements, pol);
    } else {
        // Finite interval: use standard quadrature
        using namespace boost::math::quadrature;
        
        Real error_estimate = Real(0);
        Real result = gauss_kronrod<Real, 15>::integrate(f, lower, upper,
                                                         max_refinements, tolerance,
                                                         &error_estimate);
        
        bool converged = (error_estimate < tolerance * abs(result)) ||
                        (error_estimate < tolerance);
        
        size_t evaluations = 15 * (1 << std::min(max_refinements, size_t(10)));
        return infinite_result<Real>(result, error_estimate,
                                    evaluations, converged);
    }
}

/// \brief Integrate an oscillatory function over [0, ∞)
/// 
/// Specialized routine for integrals of the form:
/// ∫₀^∞ f(x) sin(ωx) dx or ∫₀^∞ f(x) cos(ωx) dx
/// 
/// \param f Amplitude function (without oscillatory part)
/// \param omega Oscillation frequency
/// \param is_sine True for sin(ωx), false for cos(ωx)
/// \param tolerance Desired tolerance
/// \return Integration result
template <typename F, typename Real = double,
          typename Policy = boost::math::policies::policy<>>
infinite_result<Real>
integrate_oscillatory(F&& f,
                     Real omega,
                     bool is_sine = true,
                     Real tolerance = std::sqrt(std::numeric_limits<Real>::epsilon()),
                     const Policy& pol = Policy()) {
    using std::sin;
    using std::cos;
    
    if (omega <= Real(0)) {
        BOOST_MATH_THROW_EXCEPTION(std::domain_error(
            "Oscillation frequency must be positive"));
    }
    
    // Create oscillatory integrand
    auto g = [f, omega, is_sine](Real x) {
        return f(x) * (is_sine ? sin(omega * x) : cos(omega * x));
    };
    
    // Use IMT transform for oscillatory integrals
    using namespace detail;
    std::unique_ptr<domain_transform<Real>> transform = std::make_unique<imt_transform<Real>>(omega);
    
    return integrate_transformed<Real>(g, transform, tolerance, 15, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_INFINITE_HPP