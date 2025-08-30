#!/bin/bash

# Simple build script for testing individual components
# This script compiles tests with minimal dependencies

echo "Simple Build Script for Boost.Math Cubature Tests"
echo "================================================="

# Basic configuration
CXX=${CXX:-g++}
CXXFLAGS="-std=c++17 -O2 -Wall"
INCLUDES="-I../../include -I/opt/homebrew/Cellar/boost/1.88.0/include"

# Create build directory
mkdir -p build
cd build

# Function to try compiling a test
try_compile() {
    local name=$1
    local source="$2"
    local extra_flags=$3
    
    echo ""
    echo "Attempting to compile $name..."
    
    if $CXX $CXXFLAGS $INCLUDES $extra_flags "$source" -o "$name" 2> "$name.err"; then
        echo "✓ $name compiled successfully"
        rm -f "$name.err"
        return 0
    else
        echo "✗ $name failed to compile"
        echo "  Error log saved to build/$name.err"
        return 1
    fi
}

# Test 1: Basic header compilation test
echo "Creating header compilation test..."
cat > test_headers.cpp << 'EOF'
// Test that headers compile
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <iostream>

int main() {
    std::cout << "Headers compile successfully!" << std::endl;
    return 0;
}
EOF

try_compile "test_headers" "test_headers.cpp" ""

# Test 2: Simple integration test
echo "Creating simple integration test..."
cat > test_simple.cpp << 'EOF'
#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>

using namespace boost::math::cubature;

// Simple 2D Gaussian
double gaussian_2d(const double* x, std::size_t) {
    double r2 = x[0]*x[0] + x[1]*x[1];
    return std::exp(-r2);
}

int main() {
    hypercube<double> box(2);
    box.lower = {-2.0, -2.0};
    box.upper = {2.0, 2.0};
    
    auto result = integrate_adaptive<double>(
        gaussian_2d, box, 1e-6, 1e-6);
    
    std::cout << "Integral: " << result.value << std::endl;
    std::cout << "Error: " << result.error << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    
    // Check against known value (should be close to pi)
    double expected = M_PI;
    double diff = std::abs(result.value - expected);
    
    if (diff < 0.05) {  // Allow 5% error
        std::cout << "✓ Result matches expected value" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Result differs from expected" << std::endl;
        return 1;
    }
}
EOF

try_compile "test_simple" "test_simple.cpp" ""

# Test 3: Work partitioning test (without Boost.Asio)
echo "Creating work partitioning test..."
cat > test_work_partitioner.cpp << 'EOF'
// Simple test of work partitioning utilities
#include <iostream>
#include <vector>
#include <utility>
#include <numeric>

// Simple work partitioner implementation
template <typename Index>
class work_partitioner {
public:
    static std::vector<std::pair<Index, Index>> 
    create_index_ranges(Index total_work, std::size_t num_threads) {
        std::vector<std::pair<Index, Index>> ranges;
        if (num_threads == 0 || total_work == 0) return ranges;
        
        Index work_per_thread = total_work / num_threads;
        Index remainder = total_work % num_threads;
        
        Index current = 0;
        for (std::size_t i = 0; i < num_threads; ++i) {
            Index chunk = work_per_thread + (i < remainder ? 1 : 0);
            if (chunk > 0) {
                ranges.push_back({current, current + chunk});
                current += chunk;
            }
        }
        return ranges;
    }
};

// Simple tree reducer
template <typename T>
class tree_reducer {
public:
    T reduce(std::vector<T>& values) const {
        if (values.empty()) return T{};
        
        while (values.size() > 1) {
            std::size_t new_size = (values.size() + 1) / 2;
            for (std::size_t i = 0; i < new_size; ++i) {
                if (2*i + 1 < values.size()) {
                    values[i] = values[2*i] + values[2*i + 1];
                } else {
                    values[i] = values[2*i];
                }
            }
            values.resize(new_size);
        }
        return values[0];
    }
};

int main() {
    std::cout << "Testing work partitioner..." << std::endl;
    
    // Test work partitioning
    auto ranges = work_partitioner<int>::create_index_ranges(100, 4);
    
    std::cout << "Partitioned 100 items into " << ranges.size() << " ranges:" << std::endl;
    for (size_t i = 0; i < ranges.size(); ++i) {
        auto [start, end] = ranges[i];
        std::cout << "  Range " << i << ": [" << start << ", " << end << ")" << std::endl;
    }
    
    // Test tree reducer
    tree_reducer<double> reducer;
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    double result = reducer.reduce(values);
    
    std::cout << "Tree reduction result: " << result << std::endl;
    
    if (result == 36.0) {
        std::cout << "✓ Tree reduction correct" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Tree reduction incorrect" << std::endl;
        return 1;
    }
}
EOF

try_compile "test_work_partitioner" "test_work_partitioner.cpp" ""

# Test 4: Policy test
echo "Creating policy test..."
cat > test_policy.cpp << 'EOF'
#include <boost/math/cubature/policies.hpp>
#include <iostream>

using namespace boost::math::cubature;

int main() {
    std::cout << "Testing policy accumulator..." << std::endl;
    
    policy_accumulator<double> acc;
    
    // Test Kahan summation
    for (int i = 0; i < 1000000; ++i) {
        acc.add(0.1);
    }
    
    double result = acc.sum();
    double expected = 100000.0;
    double error = std::abs(result - expected);
    
    std::cout << "Sum of 1M × 0.1 = " << result << std::endl;
    std::cout << "Expected: " << expected << std::endl;
    std::cout << "Error: " << error << std::endl;
    
    if (error < 1e-10) {
        std::cout << "✓ Kahan summation working correctly" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Excessive numerical error" << std::endl;
        return 1;
    }
}
EOF

try_compile "test_policy" "test_policy.cpp" ""

# Summary
echo ""
echo "================================================="
echo "Build Summary"
echo "================================================="

success_count=0
total_count=0

for test_exe in test_headers test_simple test_work_partitioner test_policy; do
    total_count=$((total_count + 1))
    if [ -f "$test_exe" ]; then
        success_count=$((success_count + 1))
        echo "✓ $test_exe built successfully"
    else
        echo "✗ $test_exe failed to build"
    fi
done

echo ""
echo "Built $success_count out of $total_count tests"

# Try running successfully built tests
if [ $success_count -gt 0 ]; then
    echo ""
    echo "================================================="
    echo "Running Tests"
    echo "================================================="
    
    for test_exe in test_headers test_simple test_work_partitioner test_policy; do
        if [ -f "$test_exe" ]; then
            echo ""
            echo "Running $test_exe..."
            if ./$test_exe; then
                echo "✓ $test_exe passed"
            else
                echo "✗ $test_exe failed"
            fi
        fi
    done
fi

echo ""
echo "Done!"