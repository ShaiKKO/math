# Copyright (c) 2025 Boost.Math Contributors
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.14)
project(cubature_benchmarks VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Boost
find_package(Boost 1.83 REQUIRED)

# Add include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${Boost_INCLUDE_DIRS}
)

# Compiler flags for optimization
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-O3 -march=native -mtune=native)
    add_compile_options(-Wall -Wextra -Wpedantic)
    # Enable vectorization reports (optional)
    # add_compile_options(-fopt-info-vec-optimized)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/O2 /arch:AVX2)
    add_compile_options(/W4)
endif()

# Boost.Math benchmarks executable
add_executable(boost_math_benchmarks
    boost_math_benchmarks.cpp
)

target_link_libraries(boost_math_benchmarks
    Boost::boost
)

# Comparison benchmark executable
add_executable(comparison_benchmark
    comparison_benchmark.cpp
)

target_link_libraries(comparison_benchmark
    Boost::boost
)

# Results analyzer executable
add_executable(results_analyzer
    results_analyzer.cpp
)

target_link_libraries(results_analyzer
    Boost::boost
)

# Optional: Add test executable for Genz functions
add_executable(test_genz_functions
    test_genz_functions.cpp
)

target_link_libraries(test_genz_functions
    Boost::boost
)

# Install targets
install(TARGETS 
    boost_math_benchmarks
    comparison_benchmark
    results_analyzer
    RUNTIME DESTINATION bin
)

# Install headers
install(FILES
    benchmark_base.hpp
    test_functions.hpp
    DESTINATION include/boost/math/benchmark
)

# Add custom target to run all benchmarks
add_custom_target(run_benchmarks
    COMMAND boost_math_benchmarks boost_math_results
    COMMAND comparison_benchmark comparison_results
    COMMAND results_analyzer boost_math_results.csv boost_math_report.md
    COMMAND results_analyzer comparison_results_all.csv comparison_report.md
    DEPENDS boost_math_benchmarks comparison_benchmark results_analyzer
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running all benchmarks and generating reports"
)

# Add custom target for quick benchmark (subset)
add_custom_target(run_quick_benchmark
    COMMAND boost_math_benchmarks quick_results --quick
    COMMAND results_analyzer quick_results.csv quick_report.md
    DEPENDS boost_math_benchmarks results_analyzer
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running quick benchmark subset"
)

# Print configuration
message(STATUS "================================")
message(STATUS "Cubature Benchmarks Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Boost Version: ${Boost_VERSION}")
message(STATUS "================================")