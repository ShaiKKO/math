#!/bin/bash

# Build and test script for Boost.Math Cubature library
# This script compiles and runs all tests with various configurations

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "Boost.Math Cubature - Build and Test Script"
echo "=========================================="

# Configuration
BUILD_DIR="build"
SOURCE_DIR="$(dirname "$0")"
INCLUDE_DIR="$SOURCE_DIR/include"
TEST_DIR="$SOURCE_DIR/test"

# Compiler selection (can be overridden by environment)
CXX=${CXX:-g++}
echo "Using compiler: $CXX"

# Check for required tools
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${RED}Error: $1 is not installed${NC}"
        exit 1
    fi
}

check_command $CXX
check_command make

# Create build directory
echo -e "\n${YELLOW}Creating build directory...${NC}"
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Compiler flags
CXXFLAGS="-std=c++17 -O2 -Wall -Wextra -Wpedantic"
INCLUDES="-I$INCLUDE_DIR -I/usr/local/include"
LIBS="-lboost_system -lboost_thread -pthread"

# For macOS, might need different paths
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Check for Homebrew boost
    if [ -d "/opt/homebrew/opt/boost" ]; then
        INCLUDES="$INCLUDES -I/opt/homebrew/opt/boost/include"
        LIBS="-L/opt/homebrew/opt/boost/lib $LIBS"
    elif [ -d "/usr/local/opt/boost" ]; then
        INCLUDES="$INCLUDES -I/usr/local/opt/boost/include"
        LIBS="-L/usr/local/opt/boost/lib $LIBS"
    fi
fi

echo -e "\n${YELLOW}Compiler flags:${NC} $CXXFLAGS"
echo -e "${YELLOW}Includes:${NC} $INCLUDES"
echo -e "${YELLOW}Libraries:${NC} $LIBS"

# Function to compile a test
compile_test() {
    local test_name=$1
    local source_file="$TEST_DIR/${test_name}.cpp"
    local output_file="./${test_name}"
    
    if [ ! -f "$source_file" ]; then
        echo -e "${YELLOW}Skipping $test_name (source not found)${NC}"
        return 1
    fi
    
    echo -e "\n${YELLOW}Compiling $test_name...${NC}"
    
    # Special handling for certain tests
    local extra_flags=""
    if [[ "$test_name" == "multiprecision_test" ]]; then
        extra_flags="-DBOOST_MP_USE_QUAD"
    elif [[ "$test_name" == "parallel_determinism_test" ]]; then
        extra_flags="-DBOOST_ASIO_NO_DEPRECATED"
    fi
    
    if $CXX $CXXFLAGS $extra_flags $INCLUDES "$source_file" -o "$output_file" $LIBS 2>&1; then
        echo -e "${GREEN}✓ $test_name compiled successfully${NC}"
        return 0
    else
        echo -e "${RED}✗ Failed to compile $test_name${NC}"
        return 1
    fi
}

# Function to run a test
run_test() {
    local test_name=$1
    local output_file="./${test_name}"
    
    if [ ! -f "$output_file" ]; then
        return 1
    fi
    
    echo -e "\n${YELLOW}Running $test_name...${NC}"
    
    if timeout 60 "$output_file"; then
        echo -e "${GREEN}✓ $test_name passed${NC}"
        return 0
    else
        echo -e "${RED}✗ $test_name failed${NC}"
        return 1
    fi
}

# List of tests to compile and run
TESTS=(
    "genz_integration_test"
    "multiprecision_test"
    "parallel_determinism_test"
    "gm_rules_exactness_test"
    "gm97_exactness_test"
    "gm97_raw_weight_sums_test"
    "gm97_weight_sums_test"
)

# Compile all tests
echo -e "\n${YELLOW}=== Compilation Phase ===${NC}"
COMPILED_TESTS=()
FAILED_COMPILATIONS=()

for test in "${TESTS[@]}"; do
    if compile_test "$test"; then
        COMPILED_TESTS+=("$test")
    else
        FAILED_COMPILATIONS+=("$test")
    fi
done

# Report compilation results
echo -e "\n${YELLOW}=== Compilation Summary ===${NC}"
echo -e "Compiled successfully: ${#COMPILED_TESTS[@]}/${#TESTS[@]}"
if [ ${#FAILED_COMPILATIONS[@]} -gt 0 ]; then
    echo -e "${RED}Failed to compile:${NC}"
    for test in "${FAILED_COMPILATIONS[@]}"; do
        echo "  - $test"
    done
fi

# Run compiled tests
if [ ${#COMPILED_TESTS[@]} -gt 0 ]; then
    echo -e "\n${YELLOW}=== Test Execution Phase ===${NC}"
    PASSED_TESTS=()
    FAILED_TESTS=()
    
    for test in "${COMPILED_TESTS[@]}"; do
        if run_test "$test"; then
            PASSED_TESTS+=("$test")
        else
            FAILED_TESTS+=("$test")
        fi
    done
    
    # Final summary
    echo -e "\n${YELLOW}========================================${NC}"
    echo -e "${YELLOW}=== Final Test Summary ===${NC}"
    echo -e "${YELLOW}========================================${NC}"
    
    echo -e "\n${GREEN}Passed: ${#PASSED_TESTS[@]}/${#COMPILED_TESTS[@]}${NC}"
    if [ ${#PASSED_TESTS[@]} -gt 0 ]; then
        for test in "${PASSED_TESTS[@]}"; do
            echo -e "  ${GREEN}✓${NC} $test"
        done
    fi
    
    if [ ${#FAILED_TESTS[@]} -gt 0 ]; then
        echo -e "\n${RED}Failed: ${#FAILED_TESTS[@]}/${#COMPILED_TESTS[@]}${NC}"
        for test in "${FAILED_TESTS[@]}"; do
            echo -e "  ${RED}✗${NC} $test"
        done
    fi
    
    # Exit with appropriate code
    if [ ${#FAILED_TESTS[@]} -eq 0 ] && [ ${#FAILED_COMPILATIONS[@]} -eq 0 ]; then
        echo -e "\n${GREEN}All tests passed successfully!${NC}"
        exit 0
    else
        echo -e "\n${RED}Some tests failed. Please review the output above.${NC}"
        exit 1
    fi
else
    echo -e "\n${RED}No tests were compiled successfully${NC}"
    exit 1
fi