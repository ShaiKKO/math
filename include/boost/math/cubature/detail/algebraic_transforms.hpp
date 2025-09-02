// Copyright 2025 Boost.Math Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_ALGEBRAIC_TRANSFORMS_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ALGEBRAIC_TRANSFORMS_HPP

#include <boost/math/cubature/detail/domain_transforms.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <cmath>
#include <limits>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Algebraic transform for semi-infinite integrals [a, ∞)
/// 
/// Maps [a, ∞) → [0, 1] using the transformation:
/// - Forward: t = (x - a)/(1 + x - a)
/// - Inverse: x = a + t/(1 - t)
/// - Jacobian: dx/dt = 1/(1 - t)²
/// 
/// Optimal for integrands with power-law decay: f(x) ~ x^(-n), n > 1
template <typename Real>
class algebraic_transform : public domain_transform<Real> {
private:
    Real a_;  ///< Lower bound of semi-infinite interval
    
public:
    /// \brief Construct algebraic transform for [a, ∞)
    explicit algebraic_transform(Real a = Real(0)) : a_(a) {}
    
    /// \brief Map x ∈ [a, ∞) to t ∈ [0, 1]
    Real forward(Real x) const override {
        if (x <= a_) {
            return Real(0);
        }
        Real y = x - a_;
        return y / (Real(1) + y);
    }
    
    /// \brief Map t ∈ [0, 1] to x ∈ [a, ∞)
    Real inverse(Real t) const override {
        if (t <= Real(0)) {
            return a_;
        }
        if (t >= Real(1)) {
            return std::numeric_limits<Real>::infinity();
        }
        return a_ + t / (Real(1) - t);
    }
    
    /// \brief Jacobian dx/dt = 1/(1-t)²
    Real jacobian(Real t) const override {
        if (t <= Real(0) || t >= Real(1)) {
            return Real(0);
        }
        Real one_minus_t = Real(1) - t;
        return Real(1) / (one_minus_t * one_minus_t);
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(0), Real(1)};
    }
    
    /// \brief Get truncation bounds for given tolerance
    std::pair<Real, Real> truncation_bounds(Real tolerance) const override {
        using std::sqrt;
        // For power-law decay, truncate at t_min where jacobian becomes large
        Real t_min = sqrt(tolerance);
        return {t_min, Real(1) - tolerance};
    }
};

/// \brief Double algebraic transform for infinite integrals (-∞, ∞)
/// 
/// Maps (-∞, ∞) → [-1, 1] using the transformation:
/// - Forward: t = x / (1 + |x|) for symmetry
/// - Inverse: x = t / (1 - |t|)
/// - Jacobian: dx/dt = 1 / (1 - |t|)²
/// 
/// Optimal for integrands with power-law decay in both directions
template <typename Real>
class double_algebraic_transform : public domain_transform<Real> {
public:
    /// \brief Map x ∈ (-∞, ∞) to t ∈ [-1, 1]
    Real forward(Real x) const override {
        using std::abs;
        Real abs_x = abs(x);
        Real t = abs_x / (Real(1) + abs_x);
        return (x >= Real(0)) ? t : -t;
    }
    
    /// \brief Map t ∈ [-1, 1] to x ∈ (-∞, ∞)
    Real inverse(Real t) const override {
        using std::abs;
        Real abs_t = abs(t);
        if (abs_t >= Real(1)) {
            return (t > Real(0)) ? 
                std::numeric_limits<Real>::infinity() :
                -std::numeric_limits<Real>::infinity();
        }
        Real x = abs_t / (Real(1) - abs_t);
        return (t >= Real(0)) ? x : -x;
    }
    
    /// \brief Jacobian dx/dt = 1/(1 - |t|)²
    Real jacobian(Real t) const override {
        using std::abs;
        Real abs_t = abs(t);
        if (abs_t >= Real(1)) {
            return Real(0);
        }
        Real one_minus_abs_t = Real(1) - abs_t;
        return Real(1) / (one_minus_abs_t * one_minus_abs_t);
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(-1), Real(1)};
    }
    
    /// \brief Get truncation bounds for given tolerance
    std::pair<Real, Real> truncation_bounds(Real tolerance) const override {
        using std::sqrt;
        // Symmetric truncation
        Real t_max = Real(1) - sqrt(tolerance);
        return {-t_max, t_max};
    }
};

/// \brief Möbius transform for rational mappings
/// 
/// General Möbius transformation: x = (at + b)/(ct + d)
/// Used for mapping between various finite and infinite domains
template <typename Real>
class mobius_transform : public domain_transform<Real> {
private:
    Real a_, b_, c_, d_;
    Real det_;  ///< Determinant ad - bc (must be non-zero)
    
public:
    /// \brief Construct Möbius transform with given coefficients
    /// \param a, b, c, d Möbius coefficients where x = (at + b)/(ct + d)
    mobius_transform(Real a, Real b, Real c, Real d) 
        : a_(a), b_(b), c_(c), d_(d), det_(a * d - b * c) {
        if (abs(det_) < std::numeric_limits<Real>::epsilon()) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Möbius transform determinant must be non-zero"));
        }
    }
    
    /// \brief Create Möbius transform mapping [t1, t2] to [x1, x2]
    static mobius_transform create_mapping(Real t1, Real t2, Real x1, Real x2) {
        // Map [t1, t2] to [x1, x2] using cross-ratio
        Real a = x2 - x1;
        Real b = x1 * t2 - x2 * t1;
        Real c = Real(1);
        Real d = -t1;
        
        Real scale = t2 - t1;
        return mobius_transform(a / scale, b / scale, c, d / scale);
    }
    
    /// \brief Forward Möbius transform (not easily invertible in general)
    Real forward(Real x) const override {
        // Solve for t: x = (at + b)/(ct + d)
        // t = (dx - b)/(a - cx)
        Real denom = a_ - c_ * x;
        if (abs(denom) < std::numeric_limits<Real>::epsilon()) {
            return (x > Real(0)) ? Real(1) : Real(-1);
        }
        return (d_ * x - b_) / denom;
    }
    
    /// \brief Inverse Möbius transform
    Real inverse(Real t) const override {
        Real denom = c_ * t + d_;
        if (abs(denom) < std::numeric_limits<Real>::epsilon()) {
            return (t > Real(0)) ? 
                std::numeric_limits<Real>::infinity() :
                -std::numeric_limits<Real>::infinity();
        }
        return (a_ * t + b_) / denom;
    }
    
    /// \brief Jacobian of Möbius transform
    Real jacobian(Real t) const override {
        Real denom = c_ * t + d_;
        if (abs(denom) < std::numeric_limits<Real>::epsilon()) {
            return Real(0);
        }
        return abs(det_) / (denom * denom);
    }
};

/// \brief Rational transform for semi-infinite integrals with specific decay
/// 
/// Maps [a, ∞) → [0, 1] using x = a + p(1-t)/t^q
/// where p and q control the mapping characteristics
template <typename Real>
class rational_transform : public domain_transform<Real> {
private:
    Real a_;     ///< Lower bound
    Real p_;     ///< Scale parameter
    Real q_;     ///< Power parameter (controls clustering near endpoints)
    
public:
    /// \brief Construct rational transform
    /// \param a Lower bound of interval
    /// \param p Scale parameter (default 1)
    /// \param q Power parameter (default 1 for standard algebraic)
    rational_transform(Real a = Real(0), Real p = Real(1), Real q = Real(1)) 
        : a_(a), p_(p), q_(q) {
        if (p_ <= Real(0) || q_ <= Real(0)) {
            BOOST_MATH_THROW_EXCEPTION(std::domain_error(
                "Rational transform parameters must be positive"));
        }
    }
    
    /// \brief Forward transform (numerical inversion)
    Real forward(Real x) const override {
        if (x <= a_) {
            return Real(0);
        }
        
        // Solve (x - a) = p(1-t)/t^q for t
        // This requires numerical root finding in general
        // For q=1, we have the standard algebraic transform
        if (abs(q_ - Real(1)) < std::numeric_limits<Real>::epsilon()) {
            Real y = (x - a_) / p_;
            return y / (Real(1) + y);
        }
        
        // For other q values, use Newton iteration
        Real t = Real(0.5);  // Initial guess
        const int max_iter = 20;
        const Real tol = std::numeric_limits<Real>::epsilon() * Real(100);
        
        for (int i = 0; i < max_iter; ++i) {
            Real tq = std::pow(t, q_);
            Real f = p_ * (Real(1) - t) / tq - (x - a_);
            Real df = -p_ * (Real(1) + (q_ - Real(1)) * t) / (tq * t);
            
            Real dt = -f / df;
            t += dt;
            
            // Clamp to valid range
            t = std::max(Real(0.01), std::min(Real(0.99), t));
            
            if (abs(dt) < tol) {
                break;
            }
        }
        
        return t;
    }
    
    /// \brief Inverse transform
    Real inverse(Real t) const override {
        if (t <= Real(0)) {
            return a_;
        }
        if (t >= Real(1)) {
            return std::numeric_limits<Real>::infinity();
        }
        
        Real tq = std::pow(t, q_);
        return a_ + p_ * (Real(1) - t) / tq;
    }
    
    /// \brief Jacobian of the transform
    Real jacobian(Real t) const override {
        if (t <= Real(0) || t >= Real(1)) {
            return Real(0);
        }
        
        Real tq = std::pow(t, q_);
        return abs(p_) * (Real(1) + (q_ - Real(1)) * t) / (tq * t);
    }
    
    std::pair<Real, Real> finite_domain() const override {
        return {Real(0), Real(1)};
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ALGEBRAIC_TRANSFORMS_HPP