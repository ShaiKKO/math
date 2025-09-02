// Copyright (c) 2025 Boost.Math Contributors
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MATH_BENCHMARK_TEST_FUNCTIONS_HPP
#define BOOST_MATH_BENCHMARK_TEST_FUNCTIONS_HPP

#include <cmath>
#include <vector>
#include <array>
#include <functional>
#include <random>
#include <numeric>
#include <algorithm>

namespace boost { namespace math { namespace benchmark {

// Base class for Genz test functions
template <typename Real>
class genz_function_base {
protected:
    std::size_t dimension_;
    std::vector<Real> c_params_;  // difficulty parameters
    std::vector<Real> w_params_;  // shift parameters
    mutable std::size_t eval_count_;
    
public:
    genz_function_base(std::size_t dim) 
        : dimension_(dim), c_params_(dim), w_params_(dim), eval_count_(0) {
        // Initialize with moderate difficulty
        std::mt19937 gen(42);  // Fixed seed for reproducibility
        std::uniform_real_distribution<Real> c_dist(0.5, 2.0);
        std::uniform_real_distribution<Real> w_dist(0.2, 0.8);
        
        for (std::size_t i = 0; i < dim; ++i) {
            c_params_[i] = c_dist(gen);
            w_params_[i] = w_dist(gen);
        }
    }
    
    void set_difficulty(Real difficulty) {
        // Higher difficulty = larger c parameters
        for (auto& c : c_params_) {
            c *= difficulty;
        }
    }
    
    void set_parameters(const std::vector<Real>& c, const std::vector<Real>& w) {
        c_params_ = c;
        w_params_ = w;
    }
    
    std::size_t get_eval_count() const { return eval_count_; }
    void reset_eval_count() { eval_count_ = 0; }
    std::size_t dimension() const { return dimension_; }
    
    // Pure virtual function for evaluation
    virtual Real operator()(const std::vector<Real>& x) const = 0;
    
    // Overload for raw pointer (for compatibility)
    Real operator()(const Real* x) const {
        std::vector<Real> vec(x, x + dimension_);
        return operator()(vec);
    }
    
    virtual ~genz_function_base() = default;
};

// 1. Oscillatory function: ∫ cos(2πw₁ + Σc_i*x_i)
template <typename Real>
class genz_oscillatory : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_oscillatory(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real sum = 2 * M_PI * w_params_[0];
        for (std::size_t i = 0; i < dimension_; ++i) {
            sum += c_params_[i] * x[i];
        }
        return std::cos(sum);
    }
    
    // Analytical solution over [0,1]^d
    Real analytical_integral() const {
        Real product = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci = c_params_[i];
            if (std::abs(ci) > 1e-10) {
                product *= std::sin(ci) / ci;
            } else {
                product *= 1;  // limit as c->0
            }
        }
        
        // Account for the phase shift
        Real phase = 2 * M_PI * w_params_[0];
        Real sum_c = std::accumulate(c_params_.begin(), c_params_.end(), Real(0));
        
        return product * std::sin(sum_c + phase);
    }
    
    static const char* name() { return "oscillatory"; }
};

// 2. Product Peak: ∫ Π(1/(c_i^-2 + (x_i-w_i)^2))
template <typename Real>
class genz_product_peak : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_product_peak(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real product = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci_inv2 = 1.0 / (c_params_[i] * c_params_[i]);
            Real diff = x[i] - w_params_[i];
            product /= (ci_inv2 + diff * diff);
        }
        return product;
    }
    
    Real analytical_integral() const {
        Real product = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci = c_params_[i];
            Real wi = w_params_[i];
            // Use atan formula for integration
            product *= ci * (std::atan(ci * (1 - wi)) + std::atan(ci * wi));
        }
        return product;
    }
    
    static const char* name() { return "product_peak"; }
};

// 3. Corner Peak: ∫ (1 + Σc_i*x_i)^(-(d+1))
template <typename Real>
class genz_corner_peak : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_corner_peak(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real sum = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            sum += c_params_[i] * x[i];
        }
        return std::pow(sum, -(Real(dimension_) + 1));
    }
    
    Real analytical_integral() const {
        // Complex analytical formula involving multinomial expansion
        // For unit hypercube [0,1]^d
        Real result = 0;
        
        // Simplified computation for specific cases
        if (dimension_ == 1) {
            Real c = c_params_[0];
            result = (1 - std::pow(1 + c, -Real(dimension_))) / (c * dimension_);
        } else if (dimension_ == 2) {
            // Use recursive formula
            Real c0 = c_params_[0];
            Real c1 = c_params_[1];
            Real sum_c = c0 + c1;
            
            result = (1.0 / (dimension_ * c0 * c1)) * 
                    (1 - std::pow(1 + c0, -Real(dimension_)) - 
                     std::pow(1 + c1, -Real(dimension_)) + 
                     std::pow(1 + sum_c, -Real(dimension_)));
        } else {
            // General case - use numerical approximation or series expansion
            // This is a placeholder - actual implementation would use
            // multinomial theorem or recursive formulas
            result = 1.0 / (dimension_ * 
                          std::accumulate(c_params_.begin(), c_params_.end(), Real(1),
                                        std::multiplies<Real>()));
        }
        
        return result;
    }
    
    static const char* name() { return "corner_peak"; }
};

// 4. Gaussian: ∫ exp(-Σc_i^2*(x_i-w_i)^2)
template <typename Real>
class genz_gaussian : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_gaussian(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real sum = 0;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real diff = x[i] - w_params_[i];
            sum += c_params_[i] * c_params_[i] * diff * diff;
        }
        return std::exp(-sum);
    }
    
    Real analytical_integral() const {
        Real product = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci = c_params_[i];
            Real wi = w_params_[i];
            // Use error function for integration
            Real sqrt_pi = std::sqrt(M_PI);
            product *= (sqrt_pi / (2 * ci)) * 
                      (std::erf(ci * (1 - wi)) + std::erf(ci * wi));
        }
        return product;
    }
    
    static const char* name() { return "gaussian"; }
};

// 5. Continuous (C⁰): ∫ exp(-Σc_i*|x_i-w_i|)
template <typename Real>
class genz_continuous : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_continuous(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real sum = 0;
        for (std::size_t i = 0; i < dimension_; ++i) {
            sum += c_params_[i] * std::abs(x[i] - w_params_[i]);
        }
        return std::exp(-sum);
    }
    
    Real analytical_integral() const {
        Real product = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci = c_params_[i];
            Real wi = w_params_[i];
            // Integration of exp(-c|x-w|) from 0 to 1
            product *= (2 - std::exp(-ci * wi) - std::exp(-ci * (1 - wi))) / ci;
        }
        return product;
    }
    
    static const char* name() { return "continuous"; }
};

// 6. Discontinuous: ∫ exp(Σc_i*x_i) if x_0>w_0 or x_1>w_1, else 0
template <typename Real>
class genz_discontinuous : public genz_function_base<Real> {
    using Base = genz_function_base<Real>;
    using Base::dimension_;
    using Base::c_params_;
    using Base::w_params_;
    using Base::eval_count_;
    
public:
    explicit genz_discontinuous(std::size_t dim) : Base(dim) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        
        // Check discontinuity condition
        if (dimension_ >= 1 && x[0] <= w_params_[0]) {
            if (dimension_ >= 2 && x[1] <= w_params_[1]) {
                return 0;
            }
        }
        
        Real sum = 0;
        for (std::size_t i = 0; i < dimension_; ++i) {
            sum += c_params_[i] * x[i];
        }
        return std::exp(sum);
    }
    
    Real analytical_integral() const {
        // Complex analytical solution due to discontinuity
        Real result = 1;
        
        // Product over all dimensions
        for (std::size_t i = 0; i < dimension_; ++i) {
            Real ci = c_params_[i];
            if (std::abs(ci) > 1e-10) {
                result *= (std::exp(ci) - 1) / ci;
            } else {
                result *= 1;  // limit as c->0
            }
        }
        
        // Adjust for discontinuity region
        if (dimension_ >= 1) {
            Real w0 = w_params_[0];
            Real correction = 1;
            
            for (std::size_t i = 0; i < dimension_; ++i) {
                Real ci = c_params_[i];
                if (i == 0) {
                    correction *= (std::exp(ci * w0) - 1) / (std::exp(ci) - 1);
                } else if (dimension_ >= 2 && i == 1) {
                    Real w1 = w_params_[1];
                    correction *= (std::exp(ci * w1) - 1) / (std::exp(ci) - 1);
                }
            }
            
            result *= (1 - correction);
        }
        
        return result;
    }
    
    static const char* name() { return "discontinuous"; }
};

// Factory for creating Genz test functions
template <typename Real>
class genz_test_suite {
public:
    enum class function_type {
        oscillatory,
        product_peak,
        corner_peak,
        gaussian,
        continuous,
        discontinuous
    };
    
    static std::unique_ptr<genz_function_base<Real>> 
    create_function(function_type type, std::size_t dimension) {
        switch (type) {
            case function_type::oscillatory:
                return std::make_unique<genz_oscillatory<Real>>(dimension);
            case function_type::product_peak:
                return std::make_unique<genz_product_peak<Real>>(dimension);
            case function_type::corner_peak:
                return std::make_unique<genz_corner_peak<Real>>(dimension);
            case function_type::gaussian:
                return std::make_unique<genz_gaussian<Real>>(dimension);
            case function_type::continuous:
                return std::make_unique<genz_continuous<Real>>(dimension);
            case function_type::discontinuous:
                return std::make_unique<genz_discontinuous<Real>>(dimension);
            default:
                throw std::invalid_argument("Unknown function type");
        }
    }
    
    static std::vector<function_type> all_functions() {
        return {
            function_type::oscillatory,
            function_type::product_peak,
            function_type::corner_peak,
            function_type::gaussian,
            function_type::continuous,
            function_type::discontinuous
        };
    }
    
    static const char* function_name(function_type type) {
        switch (type) {
            case function_type::oscillatory: return "oscillatory";
            case function_type::product_peak: return "product_peak";
            case function_type::corner_peak: return "corner_peak";
            case function_type::gaussian: return "gaussian";
            case function_type::continuous: return "continuous";
            case function_type::discontinuous: return "discontinuous";
            default: return "unknown";
        }
    }
    
    // Create test function with specific difficulty level
    static std::unique_ptr<genz_function_base<Real>>
    create_with_difficulty(function_type type, std::size_t dimension, Real difficulty) {
        auto func = create_function(type, dimension);
        func->set_difficulty(difficulty);
        return func;
    }
    
    // Get analytical integral for a specific function
    template <typename Function>
    static Real get_analytical_integral(Function& func) {
        if (auto* osc = dynamic_cast<genz_oscillatory<Real>*>(&func)) {
            return osc->analytical_integral();
        } else if (auto* peak = dynamic_cast<genz_product_peak<Real>*>(&func)) {
            return peak->analytical_integral();
        } else if (auto* corner = dynamic_cast<genz_corner_peak<Real>*>(&func)) {
            return corner->analytical_integral();
        } else if (auto* gauss = dynamic_cast<genz_gaussian<Real>*>(&func)) {
            return gauss->analytical_integral();
        } else if (auto* cont = dynamic_cast<genz_continuous<Real>*>(&func)) {
            return cont->analytical_integral();
        } else if (auto* disc = dynamic_cast<genz_discontinuous<Real>*>(&func)) {
            return disc->analytical_integral();
        }
        throw std::runtime_error("Unknown function type for analytical integral");
    }
};

// Additional test functions for specific scenarios

// Polynomial test function for exactness testing
template <typename Real>
class polynomial_function {
private:
    std::vector<std::size_t> powers_;
    mutable std::size_t eval_count_;
    
public:
    polynomial_function(const std::vector<std::size_t>& powers) 
        : powers_(powers), eval_count_(0) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real result = 1;
        for (std::size_t i = 0; i < x.size() && i < powers_.size(); ++i) {
            result *= std::pow(x[i], powers_[i]);
        }
        return result;
    }
    
    Real analytical_integral() const {
        Real result = 1;
        for (auto p : powers_) {
            result /= (p + 1);
        }
        return result;
    }
    
    std::size_t get_eval_count() const { return eval_count_; }
    void reset_eval_count() { eval_count_ = 0; }
};

// Singular integrand for stress testing
template <typename Real>
class singular_function {
private:
    Real alpha_;  // singularity strength
    std::size_t dimension_;
    mutable std::size_t eval_count_;
    
public:
    singular_function(std::size_t dim, Real alpha = 0.5) 
        : alpha_(alpha), dimension_(dim), eval_count_(0) {}
    
    Real operator()(const std::vector<Real>& x) const {
        ++eval_count_;
        Real prod = 1;
        for (std::size_t i = 0; i < dimension_; ++i) {
            prod *= std::pow(x[i], -alpha_);
        }
        return prod;
    }
    
    Real analytical_integral() const {
        if (alpha_ >= 1) {
            return std::numeric_limits<Real>::infinity();
        }
        return std::pow(1.0 / (1 - alpha_), dimension_);
    }
    
    std::size_t get_eval_count() const { return eval_count_; }
    void reset_eval_count() { eval_count_ = 0; }
};

}}} // namespace boost::math::benchmark

#endif // BOOST_MATH_BENCHMARK_TEST_FUNCTIONS_HPP