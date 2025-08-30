// Test sparse grid diagnostic functionality
#include <boost/math/cubature/detail/sparse_grid_impl.hpp>
#include <boost/math/cubature/regions.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

// Sharp Gaussian that causes problems
template <typename Real>
class sharp_gaussian {
    std::vector<Real> a_;
    std::vector<Real> u_;
public:
    sharp_gaussian() : a_{5.0, 5.0, 5.0}, u_{0.5, 0.5, 0.5} {}
    
    Real operator()(const Real* x, std::size_t) const {
        Real sum = 0;
        for (int i = 0; i < 3; ++i) {
            Real diff = x[i] - u_[i];
            sum += a_[i] * a_[i] * diff * diff;
        }
        return std::exp(-sum);
    }
};

// Smooth function that should work well
template <typename Real>
class smooth_function {
public:
    Real operator()(const Real* x, std::size_t) const {
        return std::cos(x[0]) * std::sin(x[1]) * std::exp(-x[2]);
    }
};

int main() {
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    std::cout << "=== Sparse Grid Diagnostic Test ===\n\n";
    
    // Test 1: Sharp Gaussian (problematic)
    std::cout << "Test 1: Sharp Gaussian (exp(-25*sum((x_i - 0.5)^2)))\n";
    std::cout << "----------------------------------------------\n";
    {
        sharp_gaussian<double> f;
        
        for (unsigned level = 3; level <= 5; ++level) {
            std::cout << "\nLevel " << level << ":\n";
            
            // Create grid with diagnostics enabled
            smolyak_grid<double> grid(3, level, true);
            auto result = grid.evaluate(f, box);
            auto diag = grid.get_diagnostics();
            
            std::cout << "  Result: " << std::scientific << std::setprecision(6) 
                      << result.value << "\n";
            std::cout << "  Status: " << (result.status == status_code::success ? "SUCCESS" : 
                                         result.status == status_code::cancelled ? "WARNING (numerical issues)" : "ERROR") << "\n";
            std::cout << "  Nodes: " << result.evaluations << "\n";
            
            std::cout << "\n  Weight Diagnostics:\n";
            std::cout << "    Positive weights: " << diag.num_positive_weights 
                      << " (sum: " << diag.sum_positive_weights << ")\n";
            std::cout << "    Negative weights: " << diag.num_negative_weights 
                      << " (sum: -" << diag.sum_negative_weights << ")\n";
            std::cout << "    Cancellation ratio: " << std::fixed << std::setprecision(2) 
                      << diag.weight_cancellation_ratio << "\n";
            std::cout << "    Has issues: " << (diag.has_weight_issues ? "YES" : "NO") << "\n";
            
            if (result.value < 0) {
                std::cout << "  *** WARNING: Negative result for positive integrand! ***\n";
            }
        }
    }
    
    std::cout << "\n\nTest 2: Smooth Function (cos(x)*sin(y)*exp(-z))\n";
    std::cout << "----------------------------------------------\n";
    {
        smooth_function<double> f;
        
        for (unsigned level = 3; level <= 5; ++level) {
            std::cout << "\nLevel " << level << ":\n";
            
            // Create grid with diagnostics enabled
            smolyak_grid<double> grid(3, level, true);
            auto result = grid.evaluate(f, box);
            auto diag = grid.get_diagnostics();
            
            std::cout << "  Result: " << std::scientific << std::setprecision(6) 
                      << result.value << "\n";
            std::cout << "  Status: " << (result.status == status_code::success ? "SUCCESS" : "WARNING") << "\n";
            std::cout << "  Nodes: " << result.evaluations << "\n";
            
            std::cout << "\n  Weight Diagnostics:\n";
            std::cout << "    Positive weights: " << diag.num_positive_weights 
                      << " (sum: " << diag.sum_positive_weights << ")\n";
            std::cout << "    Negative weights: " << diag.num_negative_weights 
                      << " (sum: -" << diag.sum_negative_weights << ")\n";
            std::cout << "    Cancellation ratio: " << std::fixed << std::setprecision(2) 
                      << diag.weight_cancellation_ratio << "\n";
            std::cout << "    Has issues: " << (diag.has_weight_issues ? "YES" : "NO") << "\n";
        }
    }
    
    std::cout << "\n=== Diagnostic Summary ===\n";
    std::cout << "- Sharp peaked functions show high weight cancellation ratios\n";
    std::cout << "- Negative results for positive integrands indicate grid inadequacy\n";
    std::cout << "- Smooth functions show lower cancellation ratios and better results\n";
    
    return 0;
}