#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/cubature/cubature.hpp>

using namespace boost::math::cubature;

int main() {
    using Real = double;
    
    // Simple test: exp(-x) on [0, inf)
    auto exp_func = [](const Real* x, std::size_t) -> Real {
        std::cout << "    Evaluating at x=" << x[0] << ", exp(-x)=" << std::exp(-x[0]) << std::endl;
        return std::exp(-x[0]);
    };
    
    std::vector<Real> lower = {0};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    std::cout << "Testing 1D [0,inf) with sparse grid level 3..." << std::endl;
    auto result = integrate_sparse_grid_infinite<Real>(exp_func, lower, upper, 3);
    
    std::cout << "Result: " << result.value << std::endl;
    std::cout << "Expected: 1" << std::endl;
    std::cout << "Status: " << static_cast<int>(result.status) << std::endl;
    
    // Test 2D case
    std::cout << "\nTesting 2D [0,1] x [0,inf)..." << std::endl;
    auto exp_2d = [](const Real* x, std::size_t) -> Real {
        std::cout << "    Evaluating at x=" << x[0] << ", y=" << x[1] 
                  << ", exp(-y)=" << std::exp(-x[1]) << std::endl;
        if (std::isnan(x[1]) || std::isinf(x[1])) {
            std::cout << "    WARNING: Invalid y value!" << std::endl;
        }
        return std::exp(-x[1]);
    };
    
    std::vector<Real> lower2 = {0, 0};
    std::vector<Real> upper2 = {1, std::numeric_limits<Real>::infinity()};
    
    auto result2 = integrate_sparse_grid_infinite<Real>(exp_2d, lower2, upper2, 2);
    std::cout << "Result: " << result2.value << std::endl;
    std::cout << "Expected: 1" << std::endl;
    
    return 0;
}