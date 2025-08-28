#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>
#include <vector>
#include <cmath>

using namespace boost::math::cubature::detail::gm;

int main() {
    // Get current implementation
    auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
    auto weights = rule_fam<7, 2, family_9_7>::weights<double>();
    
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Current implementation has " << nodes.size() << " nodes\n\n";
    
    // Expected CSV nodes (degree-7 only, on [-1,1]^2)
    struct CSVNode { double x, y, weight; std::string type; };
    std::vector<CSVNode> csv_nodes = {
        {0.000000000000000000, 0.000000000000000000, 2.914823000350649451, "center"},
        {-0.406057174738239546, 0.000000000000000000, 0.123593980320432464, "axis_l2"},
        {0.000000000000000000, -0.406057174738239546, 0.123593980320432464, "axis_l2"},
        {0.000000000000000000, 0.406057174738239546, 0.123593980320432464, "axis_l2"},
        {0.406057174738239546, 0.000000000000000000, 0.123593980320432464, "axis_l2"},
        {-0.955907315804538915, 0.000000000000000000, 0.049205508173196022, "axis_l1"},
        {0.000000000000000000, -0.955907315804538915, 0.049205508173196022, "axis_l1"},
        {0.000000000000000000, 0.955907315804538915, 0.049205508173196022, "axis_l1"},
        {0.955907315804538915, 0.000000000000000000, 0.049205508173196022, "axis_l1"},
        {-0.955907315804538915, -0.955907315804538915, 0.088785828081335019, "pair_l1_l1"},
        {-0.955907315804538915, 0.955907315804538915, 0.088785828081335019, "pair_l1_l1"},
        {0.955907315804538915, -0.955907315804538915, 0.088785828081335019, "pair_l1_l1"},
        {0.955907315804538915, 0.955907315804538915, 0.088785828081335019, "pair_l1_l1"},
        {-0.686075797561756295, -0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {-0.686075797561756295, 0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {0.686075797561756295, -0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {0.686075797561756295, 0.686075797561756295, 0.088785828081335019, "full_diag_l0"}
    };
    
    std::cout << "CSV has " << csv_nodes.size() << " nodes\n\n";
    
    // Map CSV nodes to [0,1]^2 and scale weights
    auto map = [](double x) { return (x + 1.0) / 2.0; };
    double scale = 1.0 / 4.0;
    
    std::cout << "CSV nodes mapped to [0,1]^2:\n";
    for (size_t i = 0; i < csv_nodes.size(); ++i) {
        double x0 = map(csv_nodes[i].x);
        double x1 = map(csv_nodes[i].y);
        double w = csv_nodes[i].weight * scale;
        std::cout << i << ": (" << x0 << ", " << x1 << ") w=" << w << " [" << csv_nodes[i].type << "]\n";
    }
    
    std::cout << "\nCurrent implementation nodes:\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::cout << i << ": (" << nodes[i][0] << ", " << nodes[i][1] << ") w=" << weights[i] << "\n";
    }
    
    // Test integral of x0^2
    double csv_integral = 0.0;
    for (const auto& node : csv_nodes) {
        double x0 = map(node.x);
        double w = node.weight * scale;
        csv_integral += w * (x0 * x0);
    }
    
    double impl_integral = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        impl_integral += weights[i] * (x0 * x0);
    }
    
    std::cout << "\nIntegral of x0^2:\n";
    std::cout << "CSV: " << csv_integral << "\n";
    std::cout << "Implementation: " << impl_integral << "\n";
    std::cout << "Expected: " << (1.0/3.0) << "\n";
    
    return 0;
}
