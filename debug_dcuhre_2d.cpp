#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>
#include <vector>
#include <cmath>

using namespace boost::math::cubature::detail::gm;

int main() {
    // Test the corrected 2D degree-7 rule
    auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
    auto weights = rule_fam<7, 2, family_9_7>::weights<double>();
    
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "2D degree-7 rule (DCUHRE corrected) has " << nodes.size() << " nodes\n\n";
    
    // Expected structure based on DCUHRE specification:
    // 1 center + 4 axis_l2 + 4 axis_l1 + 4 pair_l1_l1 + 4 full_diag_l0 = 17 nodes
    
    const long double l1 = 0.955907315804538915L;
    const long double l2 = 0.406057174738239546L;
    const long double l0 = 0.686075797561756295L;
    auto map = [](double x) { return (x + 1.0) / 2.0; };
    
    std::cout << "Expected node structure:\n";
    std::cout << "Center: (" << 0.5 << ", " << 0.5 << ")\n";
    std::cout << "axis_l2: (" << map(-l2) << ", 0.5), (" << map(l2) << ", 0.5), (0.5, " << map(-l2) << "), (0.5, " << map(l2) << ")\n";
    std::cout << "axis_l1: (" << map(-l1) << ", 0.5), (" << map(l1) << ", 0.5), (0.5, " << map(-l1) << "), (0.5, " << map(l1) << ")\n";
    std::cout << "pair_l1_l1: (" << map(-l1) << ", " << map(-l1) << "), (" << map(-l1) << ", " << map(l1) << "), (" << map(l1) << ", " << map(-l1) << "), (" << map(l1) << ", " << map(l1) << ")\n";
    std::cout << "full_diag_l0: (" << map(-l0) << ", " << map(-l0) << "), (" << map(-l0) << ", " << map(l0) << "), (" << map(l0) << ", " << map(-l0) << "), (" << map(l0) << ", " << map(l0) << ")\n\n";
    
    std::cout << "Actual generated nodes:\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::cout << i << ": (" << nodes[i][0] << ", " << nodes[i][1] << ") w=" << weights[i] << "\n";
    }
    
    // Check weight sum
    double weight_sum = 0.0;
    for (auto w : weights) weight_sum += w;
    std::cout << "\nWeight sum: " << weight_sum << " (expected: 1.0)\n";
    
    // Test integral of x0^2
    double integral_x0_sq = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        integral_x0_sq += weights[i] * (x0 * x0);
    }
    
    std::cout << "∫ x0² dx: " << integral_x0_sq << " (expected: 0.333333333333333)\n";
    std::cout << "Error: " << std::abs(integral_x0_sq - 1.0/3.0) << "\n";
    
    // Let's also check if the issue is with a specific orbit
    std::cout << "\nContributions by orbit:\n";
    size_t k = 0;
    
    // Center
    double center_contrib = weights[k] * (nodes[k][0] * nodes[k][0]);
    std::cout << "Center: " << center_contrib << "\n";
    k++;
    
    // axis_l2
    double axis_l2_contrib = 0.0;
    for (int i = 0; i < 4; ++i) {
        axis_l2_contrib += weights[k] * (nodes[k][0] * nodes[k][0]);
        k++;
    }
    std::cout << "axis_l2: " << axis_l2_contrib << "\n";
    
    // axis_l1
    double axis_l1_contrib = 0.0;
    for (int i = 0; i < 4; ++i) {
        axis_l1_contrib += weights[k] * (nodes[k][0] * nodes[k][0]);
        k++;
    }
    std::cout << "axis_l1: " << axis_l1_contrib << "\n";
    
    // pair_l1_l1
    double pair_contrib = 0.0;
    for (int i = 0; i < 4; ++i) {
        pair_contrib += weights[k] * (nodes[k][0] * nodes[k][0]);
        k++;
    }
    std::cout << "pair_l1_l1: " << pair_contrib << "\n";
    
    // full_diag_l0
    double diag_contrib = 0.0;
    for (int i = 0; i < 4; ++i) {
        diag_contrib += weights[k] * (nodes[k][0] * nodes[k][0]);
        k++;
    }
    std::cout << "full_diag_l0: " << diag_contrib << "\n";
    
    std::cout << "Total: " << (center_contrib + axis_l2_contrib + axis_l1_contrib + pair_contrib + diag_contrib) << "\n";
    
    return 0;
}
