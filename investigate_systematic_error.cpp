#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>
#include <cmath>

using namespace boost::math::cubature::detail::gm;

int main() {
    std::cout << std::fixed << std::setprecision(15);
    
    // Test both degree-7 and degree-9 rules for comparison
    std::cout << "=== 2D Degree-7 Rule ===\n";
    {
        auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
        auto weights = rule_fam<7, 2, family_9_7>::weights<double>();
        
        double integral_x0_sq = 0.0;
        double integral_x1_sq = 0.0;
        double integral_x0_x1 = 0.0;
        double integral_const = 0.0;
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            double x0 = nodes[i][0];
            double x1 = nodes[i][1];
            double w = weights[i];
            
            integral_const += w;
            integral_x0_sq += w * (x0 * x0);
            integral_x1_sq += w * (x1 * x1);
            integral_x0_x1 += w * (x0 * x1);
        }
        
        std::cout << "∫ 1 dx: " << integral_const << " (expected: 1.0, error: " << std::abs(integral_const - 1.0) << ")\n";
        std::cout << "∫ x0² dx: " << integral_x0_sq << " (expected: 0.333333, error: " << std::abs(integral_x0_sq - 1.0/3.0) << ")\n";
        std::cout << "∫ x1² dx: " << integral_x1_sq << " (expected: 0.333333, error: " << std::abs(integral_x1_sq - 1.0/3.0) << ")\n";
        std::cout << "∫ x0*x1 dx: " << integral_x0_x1 << " (expected: 0.25, error: " << std::abs(integral_x0_x1 - 0.25) << ")\n";
        
        // Check the ratio of actual to expected for x0²
        double ratio = integral_x0_sq / (1.0/3.0);
        std::cout << "Ratio (actual/expected) for x0²: " << ratio << "\n";
        std::cout << "Error as fraction: " << (1.0/3.0 - integral_x0_sq) << " = " << (1.0/3.0 - integral_x0_sq) * 16 << "/16\n\n";
    }
    
    std::cout << "=== 2D Degree-9 Rule ===\n";
    {
        auto nodes = rule_fam<9, 2, family_9_7>::nodes<double>();
        auto weights = rule_fam<9, 2, family_9_7>::weights<double>();
        
        double integral_x0_sq = 0.0;
        double integral_x1_sq = 0.0;
        double integral_x0_x1 = 0.0;
        double integral_const = 0.0;
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            double x0 = nodes[i][0];
            double x1 = nodes[i][1];
            double w = weights[i];
            
            integral_const += w;
            integral_x0_sq += w * (x0 * x0);
            integral_x1_sq += w * (x1 * x1);
            integral_x0_x1 += w * (x0 * x1);
        }
        
        std::cout << "∫ 1 dx: " << integral_const << " (expected: 1.0, error: " << std::abs(integral_const - 1.0) << ")\n";
        std::cout << "∫ x0² dx: " << integral_x0_sq << " (expected: 0.333333, error: " << std::abs(integral_x0_sq - 1.0/3.0) << ")\n";
        std::cout << "∫ x1² dx: " << integral_x1_sq << " (expected: 0.333333, error: " << std::abs(integral_x1_sq - 1.0/3.0) << ")\n";
        std::cout << "∫ x0*x1 dx: " << integral_x0_x1 << " (expected: 0.25, error: " << std::abs(integral_x0_x1 - 0.25) << ")\n";
        
        // Check the ratio of actual to expected for x0²
        double ratio = integral_x0_sq / (1.0/3.0);
        std::cout << "Ratio (actual/expected) for x0²: " << ratio << "\n";
        std::cout << "Error as fraction: " << (1.0/3.0 - integral_x0_sq) << " = " << (1.0/3.0 - integral_x0_sq) * 16 << "/16\n\n";
    }
    
    // Test manual calculation on [-1,1]² domain to verify transformation
    std::cout << "=== Manual verification on [-1,1]² domain ===\n";
    {
        // Expected integral of x² on [-1,1]² should be:
        // ∫₋₁¹ ∫₋₁¹ x² dx dy = ∫₋₁¹ [x³/3]₋₁¹ dy = ∫₋₁¹ (1/3 - (-1/3)) dy = ∫₋₁¹ (2/3) dy = (2/3) * 2 = 4/3
        double expected_on_minus1_1 = 4.0/3.0;
        std::cout << "Expected ∫ x² on [-1,1]²: " << expected_on_minus1_1 << "\n";
        
        // Using degree-7 weights directly (before scaling)
        const long double w_center = 2.914823000350649604L;
        const long double w_axis_l2 = 0.123593980320432341L;
        const long double w_axis_l1 = 0.049205508173196015L;
        const long double w_pair_l1_l1 = 0.009708933337374202L;
        const long double w_full_diag = 0.088785828081335042L;
        
        const long double l1 = 0.955907315804538915L;
        const long double l2 = 0.406057174738239546L;
        const long double l0 = 0.686075797561756295L;
        
        double integral_raw = 0.0;
        
        // Center: (0,0)
        integral_raw += w_center * (0.0 * 0.0);
        
        // axis_l2: (±l2, 0), (0, ±l2)
        integral_raw += w_axis_l2 * (l2 * l2);   // (+l2, 0)
        integral_raw += w_axis_l2 * ((-l2) * (-l2)); // (-l2, 0)
        integral_raw += w_axis_l2 * (0.0 * 0.0); // (0, +l2)
        integral_raw += w_axis_l2 * (0.0 * 0.0); // (0, -l2)
        
        // axis_l1: (±l1, 0), (0, ±l1)
        integral_raw += w_axis_l1 * (l1 * l1);   // (+l1, 0)
        integral_raw += w_axis_l1 * ((-l1) * (-l1)); // (-l1, 0)
        integral_raw += w_axis_l1 * (0.0 * 0.0); // (0, +l1)
        integral_raw += w_axis_l1 * (0.0 * 0.0); // (0, -l1)
        
        // pair_l1_l1: (±l1, ±l1)
        integral_raw += w_pair_l1_l1 * (l1 * l1);     // (+l1, +l1)
        integral_raw += w_pair_l1_l1 * (l1 * l1);     // (+l1, -l1)
        integral_raw += w_pair_l1_l1 * ((-l1) * (-l1)); // (-l1, +l1)
        integral_raw += w_pair_l1_l1 * ((-l1) * (-l1)); // (-l1, -l1)
        
        // full_diag_l0: (±l0, ±l0)
        integral_raw += w_full_diag * (l0 * l0);     // (+l0, +l0)
        integral_raw += w_full_diag * (l0 * l0);     // (+l0, -l0)
        integral_raw += w_full_diag * ((-l0) * (-l0)); // (-l0, +l0)
        integral_raw += w_full_diag * ((-l0) * (-l0)); // (-l0, -l0)
        
        std::cout << "Raw integral on [-1,1]²: " << integral_raw << "\n";
        std::cout << "Error on [-1,1]²: " << std::abs(integral_raw - expected_on_minus1_1) << "\n";
        std::cout << "Scaled to [0,1]² (÷4): " << integral_raw / 4.0 << "\n";
    }
    
    return 0;
}
