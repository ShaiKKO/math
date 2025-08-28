#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

int main() {
    std::cout << std::fixed << std::setprecision(15);
    
    // Test hypothesis: What if the weights are already for [0,1]² domain?
    // In that case, we shouldn't scale by 1/4
    
    std::cout << "=== Testing if weights are already for [0,1]² domain ===\n";
    
    // Get the current implementation (which scales by 1/4)
    auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
    auto weights = rule_fam<7, 2, family_9_7>::weights<double>();
    
    // Test with weights scaled by 4 (i.e., remove the 1/4 scaling)
    std::cout << "Current implementation (scaled by 1/4):\n";
    double integral_current = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        integral_current += weights[i] * (x0 * x0);
    }
    std::cout << "∫ x0² dx: " << integral_current << " (error: " << std::abs(integral_current - 1.0/3.0) << ")\n";
    
    // Test with weights multiplied by 4 (no scaling)
    std::cout << "\nWith weights multiplied by 4 (no domain scaling):\n";
    double integral_no_scale = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        integral_no_scale += (weights[i] * 4.0) * (x0 * x0);
    }
    std::cout << "∫ x0² dx: " << integral_no_scale << " (error: " << std::abs(integral_no_scale - 1.0/3.0) << ")\n";
    
    // Check weight sums
    double weight_sum_current = 0.0;
    double weight_sum_no_scale = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        weight_sum_current += weights[i];
        weight_sum_no_scale += weights[i] * 4.0;
    }
    std::cout << "\nWeight sums:\n";
    std::cout << "Current (scaled): " << weight_sum_current << "\n";
    std::cout << "No scaling: " << weight_sum_no_scale << "\n";
    
    // Test the raw DCUHRE weights directly on the mapped nodes
    std::cout << "\n=== Testing raw DCUHRE weights on [0,1]² nodes ===\n";
    
    const long double w_center = 2.914823000350649604L;
    const long double w_axis_l2 = 0.123593980320432341L;
    const long double w_axis_l1 = 0.049205508173196015L;
    const long double w_pair_l1_l1 = 0.009708933337374202L;
    const long double w_full_diag = 0.088785828081335042L;
    
    double integral_raw_on_01 = 0.0;
    double weight_sum_raw = 0.0;
    
    // Apply raw weights to [0,1]² nodes
    size_t k = 0;
    
    // Center
    integral_raw_on_01 += w_center * (nodes[k][0] * nodes[k][0]);
    weight_sum_raw += w_center;
    k++;
    
    // axis_l2 (4 nodes)
    for (int i = 0; i < 4; ++i) {
        integral_raw_on_01 += w_axis_l2 * (nodes[k][0] * nodes[k][0]);
        weight_sum_raw += w_axis_l2;
        k++;
    }
    
    // axis_l1 (4 nodes)
    for (int i = 0; i < 4; ++i) {
        integral_raw_on_01 += w_axis_l1 * (nodes[k][0] * nodes[k][0]);
        weight_sum_raw += w_axis_l1;
        k++;
    }
    
    // pair_l1_l1 (4 nodes)
    for (int i = 0; i < 4; ++i) {
        integral_raw_on_01 += w_pair_l1_l1 * (nodes[k][0] * nodes[k][0]);
        weight_sum_raw += w_pair_l1_l1;
        k++;
    }
    
    // full_diag_l0 (4 nodes)
    for (int i = 0; i < 4; ++i) {
        integral_raw_on_01 += w_full_diag * (nodes[k][0] * nodes[k][0]);
        weight_sum_raw += w_full_diag;
        k++;
    }
    
    std::cout << "Raw DCUHRE weights on [0,1]² nodes:\n";
    std::cout << "∫ x0² dx: " << integral_raw_on_01 << " (error: " << std::abs(integral_raw_on_01 - 1.0/3.0) << ")\n";
    std::cout << "Weight sum: " << weight_sum_raw << "\n";
    
    return 0;
}
