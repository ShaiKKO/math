#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

struct Node {
    double x, y, weight;
    std::string type;
};

int main() {
    // CSV data for 2D degree-7 (manually extracted)
    std::vector<Node> csv_nodes = {
        {0.000000000000000000, 0.000000000000000000, 2.914823000350649451, "center"},
        {-0.406057174738239546, 0.000000000000000000, 0.123593980320432464, "axis_l2"},
        {0.000000000000000000, -0.406057174738239546, 0.123593980320432464, "axis_l2"},
        {0.000000000000000000, 0.406057174738239546, 0.123593980320432464, "axis_l2"},
        {0.406057174738239546, 0.000000000000000000, 0.123593980320432464, "axis_l2"},
        {-0.955907315804538915, 0.000000000000000000, 0.049205508173196022, "axis_l1"},
        {0.000000000000000000, -0.955907315804538915, 0.049205508173196022, "axis_l1"},
        {0.000000000000000000, 0.955907315804538915, 0.049205508173196022, "axis_l1"},
        {0.955907315804538915, 0.000000000000000000, 0.049205508173196022, "axis_l1"},
        {-0.955907315804538915, -0.955907315804538915, 0.009708933337374199, "pair_l1_l1"},
        {-0.955907315804538915, 0.955907315804538915, 0.009708933337374199, "pair_l1_l1"},
        {0.955907315804538915, -0.955907315804538915, 0.009708933337374199, "pair_l1_l1"},
        {0.955907315804538915, 0.955907315804538915, 0.009708933337374199, "pair_l1_l1"},
        {-0.686075797561756295, -0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {-0.686075797561756295, 0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {0.686075797561756295, -0.686075797561756295, 0.088785828081335019, "full_diag_l0"},
        {0.686075797561756295, 0.686075797561756295, 0.088785828081335019, "full_diag_l0"}
    };
    
    std::cout << std::fixed << std::setprecision(15);
    
    // Check weight sum on [-1,1]^2 (should be 4)
    double weight_sum = 0.0;
    double integral_x_squared_on_minus1_1 = 0.0;
    for (const auto& node : csv_nodes) {
        weight_sum += node.weight;
        integral_x_squared_on_minus1_1 += node.weight * (node.x * node.x);
    }
    std::cout << "Weight sum on [-1,1]^2: " << weight_sum << " (expected: 4.0)\n";
    std::cout << "∫ x² on [-1,1]²: " << integral_x_squared_on_minus1_1 << " (expected: 4/3 = " << 4.0/3.0 << ")\n";
    std::cout << "Error on [-1,1]²: " << std::abs(integral_x_squared_on_minus1_1 - 4.0/3.0) << "\n";
    
    // Map to [0,1]^2 and scale weights
    auto map = [](double x) { return (x + 1.0) / 2.0; };
    double scale = 1.0;  // Try no scaling - maybe weights are already for [0,1]^2
    
    double weight_sum_01 = 0.0;
    double integral_constant = 0.0;
    double integral_x0 = 0.0;
    double integral_x1 = 0.0;
    double integral_x0_squared = 0.0;
    double integral_x1_squared = 0.0;
    double integral_x0_x1 = 0.0;

    std::cout << "\nTesting multiple polynomials:\n";
    for (const auto& node : csv_nodes) {
        double x0_01 = map(node.x);
        double x1_01 = map(node.y);
        double weight_01 = node.weight * scale;

        weight_sum_01 += weight_01;
        integral_constant += weight_01 * 1.0;
        integral_x0 += weight_01 * x0_01;
        integral_x1 += weight_01 * x1_01;
        integral_x0_squared += weight_01 * (x0_01 * x0_01);
        integral_x1_squared += weight_01 * (x1_01 * x1_01);
        integral_x0_x1 += weight_01 * (x0_01 * x1_01);
    }

    std::cout << "\nResults:\n";
    std::cout << "Weight sum on [0,1]^2: " << weight_sum_01 << " (expected: 1.0)\n";
    std::cout << "∫ 1 dx: " << integral_constant << " (expected: 1.0, error: " << std::abs(integral_constant - 1.0) << ")\n";
    std::cout << "∫ x0 dx: " << integral_x0 << " (expected: 0.5, error: " << std::abs(integral_x0 - 0.5) << ")\n";
    std::cout << "∫ x1 dx: " << integral_x1 << " (expected: 0.5, error: " << std::abs(integral_x1 - 0.5) << ")\n";
    std::cout << "∫ x0² dx: " << integral_x0_squared << " (expected: 0.333333, error: " << std::abs(integral_x0_squared - 1.0/3.0) << ")\n";
    std::cout << "∫ x1² dx: " << integral_x1_squared << " (expected: 0.333333, error: " << std::abs(integral_x1_squared - 1.0/3.0) << ")\n";
    std::cout << "∫ x0*x1 dx: " << integral_x0_x1 << " (expected: 0.25, error: " << std::abs(integral_x0_x1 - 0.25) << ")\n";
    
    return 0;
}
