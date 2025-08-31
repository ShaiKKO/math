#include <boost/math/cubature/sparse_grid.hpp>
#include <iostream>

int main() {
    using boost::math::cubature::hypercube;
    using boost::math::cubature::integrate_sparse_grid;
    
    auto f = [](const double* x, std::size_t) -> double {
        return x[0] + x[1];
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    auto res = integrate_sparse_grid<double>(f, box, 0);
    
    std::cout << "Status: " << static_cast<int>(res.status) << "\n";
    std::cout << "Value: " << res.value << "\n";
    std::cout << "Evaluations: " << res.evaluations << "\n";
    
    return 0;
}
