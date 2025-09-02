// Copyright 2025 Boost.Math Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_DOMAIN_TRANSFORMS_HPP
#define BOOST_MATH_CUBATURE_DETAIL_DOMAIN_TRANSFORMS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <algorithm>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Domain types for integration
enum class domain_type {
    finite,          ///< Finite interval [a, b]
    semi_infinite,   ///< Semi-infinite interval [a, ∞) or (-∞, b]
    infinite,        ///< Infinite interval (-∞, ∞)
    ray,             ///< Ray from origin [0, ∞)
    periodic         ///< Periodic domain [0, 2π]
};

/// \brief Transform decay characteristics
enum class decay_type {
    unknown,         ///< Unknown decay behavior
    exponential,     ///< Exponential decay like e^(-x)
    power_law,       ///< Power law decay like x^(-n)
    oscillatory,     ///< Oscillatory with decay
    logarithmic      ///< Logarithmic singularity
};

/// \brief Base class for domain transformations
/// 
/// Provides interface for mapping between infinite and finite domains
/// All transforms map to/from [-1, 1] or [0, 1] for standard quadrature
template <typename Real>
class domain_transform {
    static_assert(std::is_floating_point<Real>::value,
                  "Template parameter Real must be a floating-point type");

public:
    using value_type = Real;

    /// \brief Virtual destructor for polymorphic use
    virtual ~domain_transform() = default;

    /// \brief Forward transform: infinite domain → finite domain
    /// \param x Point in infinite domain
    /// \return Point in finite domain (typically [-1, 1] or [0, 1])
    virtual Real forward(Real x) const = 0;

    /// \brief Inverse transform: finite domain → infinite domain
    /// \param t Point in finite domain
    /// \return Point in infinite domain
    virtual Real inverse(Real t) const = 0;

    /// \brief Jacobian of the inverse transform dx/dt
    /// \param t Point in finite domain
    /// \return Jacobian value for change of variables
    virtual Real jacobian(Real t) const = 0;

    /// \brief Get the finite domain range
    /// \return Pair of [lower, upper] bounds for the finite domain
    virtual std::pair<Real, Real> finite_domain() const {
        return {Real(-1), Real(1)};
    }

    /// \brief Check if a point is in the valid finite domain
    bool in_finite_domain(Real t) const {
        auto [lower, upper] = finite_domain();
        return t >= lower && t <= upper;
    }

    /// \brief Get suggested truncation bounds for the transform
    /// \param tolerance Desired error tolerance
    /// \return Pair of [lower, upper] truncation bounds in finite domain
    virtual std::pair<Real, Real> truncation_bounds(Real /*tolerance*/) const {
        return finite_domain();
    }

    /// \brief Apply transform to integrand function
    /// \param f Original function on infinite domain
    /// \return Transformed function on finite domain including jacobian
    template <typename F>
    auto transform_integrand(F&& f) const {
        return [this, f = std::forward<F>(f)](Real t) -> Real {
            if (!in_finite_domain(t)) {
                return Real(0);
            }
            Real x = inverse(t);
            Real jac = jacobian(t);
            
            // Handle potential singularities
            if (!std::isfinite(x) || !std::isfinite(jac)) {
                return Real(0);
            }
            
            return f(x) * jac;
        };
    }

    /// \brief Estimate optimal transform parameters based on function behavior
    /// \param f Function to analyze
    /// \param sample_points Number of points to sample
    /// \return Estimated decay type
    template <typename F>
    decay_type analyze_function(F&& f, int sample_points = 10) const {
        using std::abs;
        using std::log;
        
        // Sample the function at several points
        std::vector<Real> x_vals, f_vals;
        x_vals.reserve(sample_points);
        f_vals.reserve(sample_points);
        
        // Use geometric spacing for sampling
        Real x = Real(1);
        Real factor = Real(2);
        
        for (int i = 0; i < sample_points; ++i) {
            Real fx = abs(f(x));
            if (std::isfinite(fx) && fx > Real(0)) {
                x_vals.push_back(x);
                f_vals.push_back(fx);
            }
            x *= factor;
        }
        
        if (x_vals.size() < 3) {
            return decay_type::unknown;
        }
        
        // Analyze decay rate
        Real exp_score = Real(0);
        Real pow_score = Real(0);
        
        for (size_t i = 1; i < x_vals.size(); ++i) {
            Real x1 = x_vals[i-1], x2 = x_vals[i];
            Real f1 = f_vals[i-1], f2 = f_vals[i];
            
            if (f1 > Real(0) && f2 > Real(0)) {
                // Check exponential decay: f(x) ~ e^(-ax)
                Real exp_rate = (log(f1) - log(f2)) / (x2 - x1);
                exp_score += exp_rate;
                
                // Check power law decay: f(x) ~ x^(-n)
                Real pow_rate = (log(f1) - log(f2)) / (log(x2) - log(x1));
                pow_score += abs(pow_rate);
            }
        }
        
        exp_score /= Real(x_vals.size() - 1);
        pow_score /= Real(x_vals.size() - 1);
        
        // Determine dominant decay type
        if (exp_score > Real(0.5)) {
            return decay_type::exponential;
        } else if (pow_score > Real(1)) {
            return decay_type::power_law;
        } else {
            return decay_type::unknown;
        }
    }
};

/// \brief Identity transform for finite domains (no transformation)
template <typename Real>
class identity_transform : public domain_transform<Real> {
public:
    Real forward(Real x) const override { 
        return x; 
    }
    
    Real inverse(Real t) const override { 
        return t; 
    }
    
    Real jacobian(Real) const override { 
        return Real(1); 
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(-1), Real(1)};
    }
};

/// \brief Linear transform for mapping [a, b] to [-1, 1]
template <typename Real>
class linear_transform : public domain_transform<Real> {
private:
    Real a_, b_;
    Real scale_, shift_;
    
public:
    linear_transform(Real a, Real b) 
        : a_(a), b_(b),
          scale_(Real(2) / (b - a)),
          shift_((b + a) / (b - a)) {
        if (b <= a) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Invalid interval: b must be greater than a"));
        }
    }
    
    Real forward(Real x) const override {
        return scale_ * x - shift_;
    }
    
    Real inverse(Real t) const override {
        return (t + shift_) / scale_;
    }
    
    Real jacobian(Real) const override {
        return Real(1) / scale_;
    }
};

/// \brief Transform selector for automatic transform selection
template <typename Real>
class transform_selector {
public:
    /// \brief Select optimal transform based on domain and function analysis
    template <typename F>
    static std::unique_ptr<domain_transform<Real>> 
    select_transform(domain_type domain, F&& /*f*/, Real a = 0, Real b = 0) {
        switch (domain) {
            case domain_type::finite:
                return std::make_unique<linear_transform<Real>>(a, b);
                
            case domain_type::semi_infinite:
            case domain_type::ray:
            case domain_type::infinite:
                // Will be implemented with specific transforms
                return std::make_unique<identity_transform<Real>>();
                
            default:
                return std::make_unique<identity_transform<Real>>();
        }
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_DOMAIN_TRANSFORMS_HPP