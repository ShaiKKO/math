// Copyright 2025 Colin MacRitchie/Ripple Group
// Test for Gauss-Hermite quadrature rules

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>

// Simple reimplementation to test
namespace test {

template <typename Real>
std::pair<std::vector<Real>, std::vector<Real>> compute_gauss_hermite(std::size_t n) {
    std::vector<Real> nodes(n);
    std::vector<Real> weights(n);
    
    if (n == 1) {
        nodes[0] = Real(0);
        weights[0] = std::sqrt(Real(3.14159265358979323846));
        return {nodes, weights};
    }
    
    // Compute nodes
    const Real pi = Real(3.14159265358979323846);
    
    for (std::size_t i = 0; i < (n + 1) / 2; ++i) {
        // Initial guess
        Real x;
        if (i == 0) {
            x = std::sqrt(Real(2 * n + 1)) - Real(1.85575) * std::pow(Real(2 * n + 1), Real(-1.0/6.0));
        } else {
            int k = n / 2 - i;
            Real theta = pi * (Real(4 * k + 3)) / Real(4 * n + 2);
            x = std::sqrt(Real(2 * n + 1)) * std::cos(theta);
            x *= Real(1) - Real(0.165) / Real(n) - Real(0.28) / (Real(n) * Real(n));
        }
        
        // Newton refinement
        const Real tol = Real(1e-14);
        for (int iter = 0; iter < 10; ++iter) {
            Real x_old = x;
            
            // Evaluate H_n(x) using recurrence
            Real H_prev = Real(1);  // H_0
            Real H = Real(2) * x;    // H_1
            
            for (std::size_t k = 2; k <= n; ++k) {
                Real H_next = Real(2) * x * H - Real(2 * (k - 1)) * H_prev;
                H_prev = H;
                H = H_next;
            }
            
            // Derivative H'_n = 2n * H_{n-1}
            Real H_deriv = Real(2 * n) * H_prev;
            
            x = x_old - H / H_deriv;
            
            if (std::abs(x - x_old) < tol) break;
        }
        
        nodes[i] = -x;
        nodes[n - 1 - i] = x;
    }
    
    // Sort nodes
    std::sort(nodes.begin(), nodes.end());
    
    // Compute weights: w_i = 2^(n+1) * n! * sqrt(pi) / [H'_n(x_i)]^2
    const Real sqrt_pi = std::sqrt(pi);
    
    for (std::size_t i = 0; i < n; ++i) {
        Real x = nodes[i];
        
        // Evaluate H_{n-1}(x)
        Real H_prev = Real(1);  // H_0
        Real H = Real(2) * x;    // H_1
        
        for (std::size_t k = 2; k < n; ++k) {
            Real H_next = Real(2) * x * H - Real(2 * (k - 1)) * H_prev;
            H_prev = H;
            H = H_next;
        }
        
        // For n >= 2, H_{n-1} is in H
        if (n == 1) {
            H = Real(1);  // H_0 for n=1
        }
        
        // Weight formula
        Real factorial_n = Real(1);
        for (std::size_t k = 2; k <= n; ++k) {
            factorial_n *= Real(k);
        }
        
        weights[i] = (std::pow(Real(2), n) * factorial_n * sqrt_pi) / 
                     (Real(n * n) * H * H);
    }
    
    return {nodes, weights};
}

}  // namespace test

int main() {
    std::cout << "Testing Gauss-Hermite Quadrature Rules\n";
    std::cout << "======================================\n\n";
    
    // Test for various n values
    std::vector<std::size_t> test_n = {1, 2, 3, 4, 5};
    
    for (auto n : test_n) {
        std::cout << "n = " << n << " points:\n";
        auto result = test::compute_gauss_hermite<double>(n);
        auto nodes = result.first;
        auto weights = result.second;
        
        // Display nodes and weights
        std::cout << std::fixed << std::setprecision(10);
        for (std::size_t i = 0; i < n; ++i) {
            std::cout << "  x[" << i << "] = " << std::setw(15) << nodes[i]
                      << ", w[" << i << "] = " << std::setw(15) << weights[i] << "\n";
        }
        
        // Check weight sum (should be sqrt(pi))
        double weight_sum = 0;
        for (auto w : weights) weight_sum += w;
        double expected_sum = std::sqrt(3.14159265358979323846);
        std::cout << "  Weight sum: " << weight_sum 
                  << " (expected: " << expected_sum << ")\n";
        std::cout << "  Error: " << std::scientific << std::abs(weight_sum - expected_sum) << "\n";
        
        // Test polynomial exactness
        // Integrate x^(2n-2) * exp(-x^2) which should be exact
        if (n >= 2) {
            int degree = 2 * n - 2;
            double integral = 0;
            for (std::size_t i = 0; i < n; ++i) {
                integral += std::pow(nodes[i], degree) * weights[i];
            }
            
            // Exact value for x^(2k) * exp(-x^2) over (-inf,inf)
            // = (2k-1)!! * sqrt(pi) / 2^k
            double exact = 1.0;
            for (int k = 1; k <= degree/2; ++k) {
                exact *= (2.0 * k - 1.0) / 2.0;
            }
            exact *= std::sqrt(3.14159265358979323846);
            
            std::cout << "  Polynomial test (degree " << degree << "):\n";
            std::cout << "    Computed: " << std::fixed << integral << "\n";
            std::cout << "    Exact:    " << exact << "\n";
            std::cout << "    Error:    " << std::scientific 
                      << std::abs(integral - exact) << "\n";
        }
        
        std::cout << "\n";
    }
    
    // Special test for n=3 with known values
    std::cout << "Special test n=3 (known values):\n";
    auto result3 = test::compute_gauss_hermite<double>(3);
    auto n3_nodes = result3.first;
    auto n3_weights = result3.second;
    
    double expected_nodes[3] = {-std::sqrt(1.5), 0.0, std::sqrt(1.5)};
    std::cout << "Nodes comparison:\n";
    for (int i = 0; i < 3; ++i) {
        std::cout << "  Computed: " << std::fixed << std::setprecision(10) << n3_nodes[i]
                  << ", Expected: " << expected_nodes[i]
                  << ", Error: " << std::scientific << std::abs(n3_nodes[i] - expected_nodes[i]) << "\n";
    }
    
    return 0;
}