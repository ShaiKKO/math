cmake_minimum_required(VERSION 3.14)
project(BoostMathTests)

# Copyright 2021 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

include(BoostTestJamfile OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(HAVE_BOOST_TEST)

    if (BOOST_MATH_ENABLE_CUDA)

        message(STATUS "Building boost.math with CUDA")

        find_package(CUDA REQUIRED)
        enable_language(CUDA)
        set(CMAKE_CUDA_EXTENSIONS OFF)

        enable_testing()

        boost_test_jamfile(FILE cuda_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception Boost::unit_test_framework ${CUDA_LIBRARIES} INCLUDE_DIRECTORIES ${CUDA_INCLUDE_DIRS}  )

    elseif (BOOST_MATH_ENABLE_NVRTC)

        message(STATUS "Building boost.math with NVRTC")

        find_package(CUDA REQUIRED)

        enable_testing()

        set(CUDA_nvrtc_LIBRARY /usr/local/cuda/lib64/libnvrtc.so)

        if (BOOST_MATH_NVRTC_CI_RUN)
            boost_test_jamfile(FILE nvrtc_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception ${CUDA_nvrtc_LIBRARY} ${CUDA_LIBRARIES} ${CUDA_CUDA_LIBRARY} COMPILE_DEFINITIONS BOOST_MATH_NVRTC_CI_RUN=1 INCLUDE_DIRECTORIES ${CUDA_INCLUDE_DIRS} )
        else ()
            boost_test_jamfile(FILE nvrtc_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception ${CUDA_nvrtc_LIBRARY} ${CUDA_LIBRARIES} ${CUDA_CUDA_LIBRARY} INCLUDE_DIRECTORIES ${CUDA_INCLUDE_DIRS} )
        endif()

    elseif (BOOST_MATH_ENABLE_SYCL)

        message(STATUS "Building boost.math with SYCL")

        set(CMAKE_CXX_COMPILER "icpx")
        set(CMAKE_C_COMPILER "icx")

        enable_testing()

        boost_test_jamfile(FILE sycl_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception Boost::unit_test_framework sycl COMPILE_OPTIONS -fsycl )
    else()

        boost_test(SOURCES check_cmake_version.cpp ARGUMENTS ${PROJECT_VERSION} LINK_LIBRARIES Boost::core Boost::config)
        add_subdirectory(compile_test)


        # Cubature header-only smoke test
        add_executable(cubature_smoke_test cubature_smoke_test.cpp)
        target_include_directories(cubature_smoke_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
          target_compile_options(cubature_smoke_test PRIVATE -Wall -Wextra -Werror)
        elseif (MSVC)
          target_compile_options(cubature_smoke_test PRIVATE /W4 /WX)
        endif()
        enable_testing()
        add_test(NAME cubature_smoke_test COMMAND cubature_smoke_test)

    endif()

endif()

# Fallback when BoostTestJamfile is unavailable
if(NOT HAVE_BOOST_TEST)
  add_executable(cubature_smoke_test cubature_smoke_test.cpp)
  target_include_directories(cubature_smoke_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_features(cubature_smoke_test PRIVATE cxx_std_17)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(cubature_smoke_test PRIVATE -Wall -Wextra -Werror)
  elseif (MSVC)
    target_compile_options(cubature_smoke_test PRIVATE /W4 /WX)
  endif()
  enable_testing()
  add_test(NAME cubature_smoke_test COMMAND cubature_smoke_test)
endif()

# Genz–Malik 2D exactness test (enable once 2D tables land)
add_executable(gm_rules_exactness_test gm_rules_exactness_test.cpp)
target_include_directories(gm_rules_exactness_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(gm_rules_exactness_test PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gm_rules_exactness_test PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(gm_rules_exactness_test PRIVATE /W4 /WX)
endif()
add_test(NAME gm_rules_exactness_test COMMAND gm_rules_exactness_test)

add_executable(gm97_weight_sums_test gm97_weight_sums_test.cpp)
target_include_directories(gm97_weight_sums_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(gm97_weight_sums_test PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gm97_weight_sums_test PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(gm97_weight_sums_test PRIVATE /W4 /WX)
endif()
add_test(NAME gm97_weight_sums_test COMMAND gm97_weight_sums_test)


add_executable(gm97_exactness_test gm97_exactness_test.cpp)
target_include_directories(gm97_exactness_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(gm97_exactness_test PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gm97_exactness_test PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(gm97_exactness_test PRIVATE /W4 /WX)
endif()
add_test(NAME gm97_exactness_test COMMAND gm97_exactness_test)


add_executable(gm97_raw_weight_sums_test gm97_raw_weight_sums_test.cpp)
target_include_directories(gm97_raw_weight_sums_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(gm97_raw_weight_sums_test PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gm97_raw_weight_sums_test PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(gm97_raw_weight_sums_test PRIVATE /W4 /WX)
endif()
add_test(NAME gm97_raw_weight_sums_test COMMAND gm97_raw_weight_sums_test)

add_executable(test_sparse_grid_vector test_sparse_grid_vector.cpp)
target_include_directories(test_sparse_grid_vector PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(test_sparse_grid_vector PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(test_sparse_grid_vector PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(test_sparse_grid_vector PRIVATE /W4 /WX)
endif()
add_test(NAME test_sparse_grid_vector COMMAND test_sparse_grid_vector)

add_executable(adaptive_integration_test adaptive_integration_test.cpp)
target_include_directories(adaptive_integration_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_features(adaptive_integration_test PRIVATE cxx_std_14)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(adaptive_integration_test PRIVATE -Wall -Wextra -Werror)
elseif (MSVC)
  target_compile_options(adaptive_integration_test PRIVATE /W4 /WX)
endif()
add_test(NAME adaptive_integration_test COMMAND adaptive_integration_test)



