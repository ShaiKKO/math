# Copyright 2025 Colin MacRitchie/Ripple Group
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.14)
project(BoostMathCubatureTests)

# Find Boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../../include)
include_directories(${Boost_INCLUDE_DIRS})

# Enable testing
enable_testing()

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    else()
        add_compile_options(-O3)
    endif()
elseif(MSVC)
    add_compile_options(/W4)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(/Od)
    else()
        add_compile_options(/O2)
    endif()
endif()

# Test executables
set(CUBATURE_TESTS
    test_adaptive
    test_sparse_grid
    test_transforms
    test_qmc
    test_genz_malik
    test_simplex
)

# Add test executables
foreach(test ${CUBATURE_TESTS})
    add_executable(${test} ${test}.cpp)
    
    # Link with Boost.Test
    target_link_libraries(${test} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    
    # Add to CTest
    add_test(NAME cubature_${test} COMMAND ${test})
    
    # Set test properties
    set_tests_properties(cubature_${test} PROPERTIES
        TIMEOUT 60
        LABELS "cubature"
    )
endforeach()

# Add standalone tests (without Boost.Test dependency)
set(STANDALONE_TESTS
    test_adaptive_standalone
    test_qmc_standalone
)

foreach(test ${STANDALONE_TESTS})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}.cpp)
        add_executable(${test} ${test}.cpp)
        # No Boost.Test linking needed for standalone tests
        add_test(NAME cubature_${test} COMMAND ${test})
        set_tests_properties(cubature_${test} PROPERTIES
            TIMEOUT 60
            LABELS "cubature standalone"
        )
    endif()
endforeach()

# Add a custom target to run all cubature tests
add_custom_target(cubature_test
    COMMAND ${CMAKE_CTEST_COMMAND} -L cubature --output-on-failure
    DEPENDS ${CUBATURE_TESTS}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running cubature tests"
)