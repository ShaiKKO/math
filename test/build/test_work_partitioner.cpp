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
