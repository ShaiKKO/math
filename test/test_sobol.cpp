#include <boost/random/sobol.hpp>
#include <iostream>

int main() {
    // Test 2D Sobol
    boost::random::sobol_engine<std::uint64_t, 2> sobol;
    
    for (int i = 0; i < 5; ++i) {
        auto point = sobol();
        std::cout << "Point " << i << ": ";
        for (int d = 0; d < 2; ++d) {
            std::cout << point[d] << " ";
        }
        std::cout << "\n";
    }
    
    return 0;
}
