// Copyright 2025 Boost.Math Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_EXPONENTIAL_TRANSFORMS_HPP
#define BOOST_MATH_CUBATURE_DETAIL_EXPONENTIAL_TRANSFORMS_HPP

#include <boost/math/cubature/detail/domain_transforms.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <cmath>
#include <limits>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Tanh-sinh (double exponential) transform for (-∞, ∞)
/// 
/// Maps (-∞, ∞) → [-1, 1] using the transformation:
/// - Inverse: x = tanh(π/2 sinh(t))
/// - Jacobian: dx/dt = (π/2) cosh(t) / cosh²(π/2 sinh(t))
/// 
/// Optimal for:
/// - Functions with endpoint singularities
/// - Functions with exponential or faster decay
/// - Highly oscillatory integrands with decay
/// 
/// Features exponential clustering of points near endpoints
template <typename Real>
class tanh_sinh_transform : public domain_transform<Real> {
private:
    static constexpr Real pi_half = boost::math::constants::pi<Real>() / Real(2);
    
public:
    /// \brief Forward transform: x ∈ (-∞, ∞) to t ∈ [-1, 1]
    /// Note: This requires numerical inversion in general
    Real forward(Real x) const override {
        using std::abs;
        using std::log;
        using std::sqrt;
        
        // For |x| close to 1, use asymptotic approximation
        if (abs(x) > Real(0.999)) {
            Real sign_x = (x > Real(0)) ? Real(1) : Real(-1);
            Real one_minus_abs_x = Real(1) - abs(x);
            
            if (one_minus_abs_x < std::numeric_limits<Real>::epsilon()) {
                return sign_x * Real(4);  // Practical infinity for t
            }
            
            // Asymptotic: t ≈ log(log(2/(1-|x|)) * 2/π) for x → ±1
            Real arg = log(Real(2) / one_minus_abs_x) * Real(2) / boost::math::constants::pi<Real>();
            if (arg > Real(1)) {
                return sign_x * log(arg);
            }
        }
        
        // For moderate x, use Newton iteration
        Real t = Real(0);  // Initial guess
        const int max_iter = 20;
        const Real tol = std::numeric_limits<Real>::epsilon() * Real(100);
        
        for (int i = 0; i < max_iter; ++i) {
            Real sinh_t = sinh(t);
            Real cosh_t = cosh(t);
            Real arg = pi_half * sinh_t;
            Real tanh_arg = tanh(arg);
            
            Real f = tanh_arg - x;
            Real df = pi_half * cosh_t / (cosh(arg) * cosh(arg));
            
            Real dt = -f / df;
            t += dt;
            
            if (abs(dt) < tol) {
                break;
            }
        }
        
        return t;
    }
    
    /// \brief Inverse transform: t ∈ [-1, 1] to x ∈ (-∞, ∞)
    Real inverse(Real t) const override {
        using std::sinh;
        using std::tanh;
        using std::abs;
        
        // Handle extreme values
        if (abs(t) > Real(4)) {
            return (t > Real(0)) ? Real(1) : Real(-1);
        }
        
        Real sinh_t = sinh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (abs(arg) > Real(350)) {  // log(DBL_MAX) ≈ 700, be conservative
            return (arg > Real(0)) ? Real(1) : Real(-1);
        }
        
        return tanh(arg);
    }
    
    /// \brief Jacobian dx/dt
    Real jacobian(Real t) const override {
        using std::sinh;
        using std::cosh;
        using std::abs;
        
        // Handle extreme values
        if (abs(t) > Real(4)) {
            return Real(0);
        }
        
        Real sinh_t = sinh(t);
        Real cosh_t = cosh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (abs(arg) > Real(350)) {
            return Real(0);
        }
        
        Real cosh_arg = cosh(arg);
        return pi_half * cosh_t / (cosh_arg * cosh_arg);
    }
    
    std::pair<Real, Real> finite_domain() const override {
        // Use restricted domain for practical computation
        return {Real(-4), Real(4)};
    }
    
    /// \brief Get optimal truncation bounds based on tolerance
    std::pair<Real, Real> truncation_bounds(Real tolerance) const override {
        using std::log;
        using std::abs;
        
        // Truncation: |t| < log(log(2/tol) * 2/π)
        Real log_tol = -log(tolerance);
        if (log_tol < Real(1)) {
            return {Real(-1), Real(1)};
        }
        
        Real t_max = log(log_tol * Real(2) / boost::math::constants::pi<Real>());
        t_max = std::min(t_max, Real(4));  // Practical limit
        
        return {-t_max, t_max};
    }
};

/// \brief Exp-sinh transform for [0, ∞)
/// 
/// Maps [0, ∞) → [-1, 1] using the transformation:
/// - Inverse: x = exp(π/2 sinh(t))
/// - Jacobian: dx/dt = x * (π/2) cosh(t)
/// 
/// Optimal for:
/// - Functions with exponential decay
/// - Functions with singularities at x = 0
/// - Log-scale integration
template <typename Real>
class exp_sinh_transform : public domain_transform<Real> {
private:
    static constexpr Real pi_half = boost::math::constants::pi<Real>() / Real(2);
    Real a_;  ///< Lower bound (shift parameter)
    
public:
    /// \brief Construct exp-sinh transform for [a, ∞)
    explicit exp_sinh_transform(Real a = Real(0)) : a_(a) {}
    
    /// \brief Forward transform: x ∈ [a, ∞) to t ∈ [-1, 1]
    Real forward(Real x) const override {
        using std::log;
        using std::asinh;
        
        if (x <= a_) {
            return Real(-4);  // Practical negative infinity
        }
        
        Real y = x - a_;
        if (y < std::numeric_limits<Real>::epsilon()) {
            return Real(-4);
        }
        
        // Solve exp(π/2 sinh(t)) = y for t
        // sinh(t) = (2/π) log(y)
        // t = asinh((2/π) log(y))
        Real log_y = log(y);
        Real sinh_t = (Real(2) / boost::math::constants::pi<Real>()) * log_y;
        
        // Prevent extreme values
        if (abs(sinh_t) > Real(100)) {
            return (sinh_t > Real(0)) ? Real(4) : Real(-4);
        }
        
        return asinh(sinh_t);
    }
    
    /// \brief Inverse transform: t ∈ [-1, 1] to x ∈ [0, ∞)
    Real inverse(Real t) const override {
        using std::exp;
        using std::sinh;
        
        // Handle extreme values
        if (t < Real(-4)) {
            return a_;
        }
        
        Real sinh_t = sinh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (arg > Real(700)) {  // log(DBL_MAX)
            return std::numeric_limits<Real>::infinity();
        }
        
        // Prevent underflow
        if (arg < Real(-700)) {
            return a_;
        }
        
        return a_ + exp(arg);
    }
    
    /// \brief Jacobian dx/dt
    Real jacobian(Real t) const override {
        using std::exp;
        using std::sinh;
        using std::cosh;
        
        // Handle extreme values
        if (t < Real(-4)) {
            return Real(0);
        }
        
        Real sinh_t = sinh(t);
        Real cosh_t = cosh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (arg > Real(700)) {
            return Real(0);
        }
        
        // For very negative arg, jacobian is effectively zero
        if (arg < Real(-700)) {
            return Real(0);
        }
        
        Real x = exp(arg);
        return x * pi_half * cosh_t;
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(-4), Real(4)};
    }
    
    /// \brief Get optimal truncation bounds based on tolerance
    std::pair<Real, Real> truncation_bounds(Real tolerance) const override {
        using std::log;
        using std::asinh;
        
        // Lower bound: where x ≈ tolerance
        Real t_min = Real(-4);
        if (tolerance > std::numeric_limits<Real>::epsilon()) {
            Real log_tol = log(tolerance);
            Real sinh_t = (Real(2) / boost::math::constants::pi<Real>()) * log_tol;
            if (abs(sinh_t) < Real(100)) {
                t_min = std::max(Real(-4), asinh(sinh_t));
            }
        }
        
        // Upper bound: where jacobian becomes negligible
        Real log_tol = -log(tolerance);
        Real t_max = log(log_tol * Real(2) / boost::math::constants::pi<Real>());
        t_max = std::min(t_max, Real(4));
        
        return {t_min, t_max};
    }
};

/// \brief Sinh-sinh transform for integrals with double exponential decay
/// 
/// Maps (-∞, ∞) → [-1, 1] using x = sinh(π/2 sinh(t))
/// Optimal for integrands like exp(-x²) or exp(-|x|)
template <typename Real>
class sinh_sinh_transform : public domain_transform<Real> {
private:
    static constexpr Real pi_half = boost::math::constants::pi<Real>() / Real(2);
    
public:
    /// \brief Forward transform (requires numerical inversion)
    Real forward(Real x) const override {
        using std::asinh;
        using std::log;
        using std::abs;
        
        // For large |x|, use asymptotic approximation
        if (abs(x) > Real(100)) {
            Real sign_x = (x > Real(0)) ? Real(1) : Real(-1);
            Real log_abs_x = log(abs(x));
            Real sinh_inv = (Real(2) / boost::math::constants::pi<Real>()) * log_abs_x;
            return sign_x * asinh(sinh_inv);
        }
        
        // Newton iteration for moderate x
        Real t = asinh(x) / pi_half;  // Initial guess
        const int max_iter = 20;
        const Real tol = std::numeric_limits<Real>::epsilon() * Real(100);
        
        for (int i = 0; i < max_iter; ++i) {
            Real sinh_t = sinh(t);
            Real cosh_t = cosh(t);
            Real arg = pi_half * sinh_t;
            Real sinh_arg = sinh(arg);
            Real cosh_arg = cosh(arg);
            
            Real f = sinh_arg - x;
            Real df = pi_half * cosh_t * cosh_arg;
            
            Real dt = -f / df;
            t += dt;
            
            if (abs(dt) < tol) {
                break;
            }
        }
        
        return t;
    }
    
    /// \brief Inverse transform
    Real inverse(Real t) const override {
        using std::sinh;
        using std::abs;
        
        if (abs(t) > Real(4)) {
            return (t > Real(0)) ? 
                std::numeric_limits<Real>::infinity() :
                -std::numeric_limits<Real>::infinity();
        }
        
        Real sinh_t = sinh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (abs(arg) > Real(350)) {
            return (arg > Real(0)) ? 
                std::numeric_limits<Real>::infinity() :
                -std::numeric_limits<Real>::infinity();
        }
        
        return sinh(arg);
    }
    
    /// \brief Jacobian dx/dt
    Real jacobian(Real t) const override {
        using std::sinh;
        using std::cosh;
        using std::abs;
        
        if (abs(t) > Real(4)) {
            return Real(0);
        }
        
        Real sinh_t = sinh(t);
        Real cosh_t = cosh(t);
        Real arg = pi_half * sinh_t;
        
        // Prevent overflow
        if (abs(arg) > Real(350)) {
            return Real(0);
        }
        
        Real cosh_arg = cosh(arg);
        return pi_half * cosh_t * cosh_arg;
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(-4), Real(4)};
    }
    
    std::pair<Real, Real> truncation_bounds(Real tolerance) const override {
        using std::log;
        
        Real log_tol = -log(tolerance);
        Real t_max = log(log_tol * Real(2) / boost::math::constants::pi<Real>());
        t_max = std::min(t_max, Real(4));
        
        return {-t_max, t_max};
    }
};

/// \brief IMT (Iri-Mori-Takahasi) transform for oscillatory integrals
/// 
/// Specialized transform for integrals of the form ∫₀^∞ f(x)sin(ωx) dx
/// Maps [0, ∞) → [0, 1] with clustering near oscillation nodes
template <typename Real>
class imt_transform : public domain_transform<Real> {
private:
    Real omega_;  ///< Oscillation frequency
    Real alpha_;  ///< Decay parameter
    
public:
    /// \brief Construct IMT transform
    /// \param omega Oscillation frequency
    /// \param alpha Decay parameter (default 1)
    imt_transform(Real omega, Real alpha = Real(1)) 
        : omega_(omega), alpha_(alpha) {
        if (omega_ <= Real(0) || alpha_ <= Real(0)) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "IMT transform parameters must be positive"));
        }
    }
    
    Real forward(Real x) const override {
        using std::exp;
        if (x <= Real(0)) {
            return Real(0);
        }
        
        // t = 1 - exp(-alpha * omega * x)
        Real arg = alpha_ * omega_ * x;
        if (arg > Real(700)) {
            return Real(1);
        }
        
        return Real(1) - exp(-arg);
    }
    
    Real inverse(Real t) const override {
        using std::log;
        using boost::math::log1p;
        
        if (t <= Real(0)) {
            return Real(0);
        }
        if (t >= Real(1)) {
            return std::numeric_limits<Real>::infinity();
        }
        
        // x = -log(1 - t) / (alpha * omega)
        return -log1p(-t) / (alpha_ * omega_);
    }
    
    Real jacobian(Real t) const override {
        if (t <= Real(0) || t >= Real(1)) {
            return Real(0);
        }
        
        return Real(1) / ((Real(1) - t) * alpha_ * omega_);
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(0), Real(1)};
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_EXPONENTIAL_TRANSFORMS_HPP