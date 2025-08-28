#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

int main() {
    // Test 2D degree-7 rule
    auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
    auto weights = rule_fam<7, 2, family_9_7>::weights<double>();

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "2D degree-7 rule has " << nodes.size() << " nodes\n";
    
    // Check weight sum
    double weight_sum = 0.0;
    for (auto w : weights) weight_sum += w;
    std::cout << "Weight sum: " << weight_sum << " (should be 1.0)\n";
    
    // Test various polynomials
    double integral_constant = 0.0;
    double integral_x0 = 0.0;
    double integral_x1 = 0.0;
    double integral_x0_squared = 0.0;
    double integral_x1_squared = 0.0;
    double integral_x0_x1 = 0.0;
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        double x1 = nodes[i][1];
        double w = weights[i];
        
        integral_constant += w * 1.0;
        integral_x0 += w * x0;
        integral_x1 += w * x1;
        integral_x0_squared += w * (x0 * x0);
        integral_x1_squared += w * (x1 * x1);
        integral_x0_x1 += w * (x0 * x1);
    }
    
    std::cout << "\nResults for 2D degree-7:\n";
    std::cout << "∫ 1 dx: " << integral_constant << " (expected: 1.0, error: " << std::abs(integral_constant - 1.0) << ")\n";
    std::cout << "∫ x0 dx: " << integral_x0 << " (expected: 0.5, error: " << std::abs(integral_x0 - 0.5) << ")\n";
    std::cout << "∫ x1 dx: " << integral_x1 << " (expected: 0.5, error: " << std::abs(integral_x1 - 0.5) << ")\n";
    std::cout << "∫ x0² dx: " << integral_x0_squared << " (expected: 0.333333, error: " << std::abs(integral_x0_squared - 1.0/3.0) << ")\n";
    std::cout << "∫ x1² dx: " << integral_x1_squared << " (expected: 0.333333, error: " << std::abs(integral_x1_squared - 1.0/3.0) << ")\n";
    std::cout << "∫ x0*x1 dx: " << integral_x0_x1 << " (expected: 0.25, error: " << std::abs(integral_x0_x1 - 0.25) << ")\n";
    
    return 0;
}
