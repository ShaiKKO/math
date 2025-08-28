#include <iostream>
#include <iomanip>
#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

int main() {
    // Get 2D degree-7 rule
    auto nodes = rule_fam<7, 2, family_9_7>::nodes<double>();
    auto weights = rule_fam<7, 2, family_9_7>::weights<double>();

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "2D degree-7 rule has " << nodes.size() << " nodes\n";

    // Check weight sum
    double weight_sum = 0.0;
    for (auto w : weights) weight_sum += w;
    std::cout << "Weight sum: " << weight_sum << " (should be 1.0)\n\n";

    // Show the raw constants
    const long double l1 = 0.955907315804538915L;
    const long double l2 = 0.406057174738239546L;
    const long double l0 = 0.686075797561756295L;
    std::cout << "Constants: l1=" << l1 << " l2=" << l2 << " l0=" << l0 << "\n\n";

    // Map function
    auto map = [](long double x) { return (x + 1.0L) / 2.0L; };

    // Show expected CSV nodes (mapped to [0,1]^2) and their types
    std::cout << "Expected nodes from CSV (mapped to [0,1]^2):\n";
    std::cout << "Center: (" << 0.5 << ", " << 0.5 << ")\n";
    std::cout << "axis_l2: (" << map(-l2) << ", " << 0.5 << "), (" << 0.5 << ", " << map(-l2) << "), (" << 0.5 << ", " << map(l2) << "), (" << map(l2) << ", " << 0.5 << ")\n";
    std::cout << "axis_l1: (" << map(-l1) << ", " << 0.5 << "), (" << 0.5 << ", " << map(-l1) << "), (" << 0.5 << ", " << map(l1) << "), (" << map(l1) << ", " << 0.5 << ")\n";
    std::cout << "pair_l1_l1: (" << map(-l1) << ", " << map(-l1) << "), (" << map(-l1) << ", " << map(l1) << "), (" << map(l1) << ", " << map(-l1) << "), (" << map(l1) << ", " << map(l1) << ")\n";
    std::cout << "full_diag_l0: (" << map(-l0) << ", " << map(-l0) << "), (" << map(-l0) << ", " << map(l0) << "), (" << map(l0) << ", " << map(-l0) << "), (" << map(l0) << ", " << map(l0) << ")\n\n";

    // Show actual generated nodes and compute integral
    std::cout << "Actually generated nodes:\n";
    double integral = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double x0 = nodes[i][0];
        double x1 = nodes[i][1];
        double f_val = x0 * x0;  // x0^2
        double contrib = weights[i] * f_val;
        integral += contrib;
        std::cout << "Node " << i << ": (" << x0 << ", " << x1 << ") weight=" << weights[i] << " f(x0^2)=" << f_val << " contrib=" << contrib << "\n";
    }

    std::cout << "\nIntegral of x0^2: " << integral << " (expected: 0.333333333333333)\n";
    std::cout << "Error: " << std::abs(integral - 1.0/3.0) << "\n";

    // Let's also manually verify a few key contributions
    std::cout << "\nManual verification:\n";
    std::cout << "Center contrib: " << weights[0] << " * " << (0.5*0.5) << " = " << weights[0] * 0.25 << "\n";

    return 0;
}
