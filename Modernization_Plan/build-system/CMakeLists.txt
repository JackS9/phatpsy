# PHATPSY Modern Platform
# CMake configuration for C++/Python conversion

cmake_minimum_required(VERSION 3.16)
project(PHATPSY VERSION 2.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_GPU_SUPPORT "Build GPU acceleration support" OFF)
option(BUILD_WEBMO_PLUGIN "Build WebMO integration plugin" ON)
option(BUILD_TESTS "Build test suite" ON)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# Optional dependencies
if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(pybind11 REQUIRED)
endif()

if(BUILD_GPU_SUPPORT)
    find_package(CUDA QUIET)
    if(NOT CUDA_FOUND)
        find_package(hip QUIET)
    endif()
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Core library subdirectories
add_subdirectory(src/core)
add_subdirectory(src/integrals)
add_subdirectory(src/scf)
add_subdirectory(src/analysis)

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

# WebMO plugin
if(BUILD_WEBMO_PLUGIN)
    add_subdirectory(webmo)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS phatpsy-core phatpsy-integrals phatpsy-scf phatpsy-analysis
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Export targets for find_package
install(EXPORT phatpsy-targets
        FILE phatpsy-targets.cmake
        NAMESPACE phatpsy::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phatpsy)

# Package configuration
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/phatpsy-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/phatpsy-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phatpsy)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/phatpsy-config-version.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/phatpsy-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/phatpsy-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phatpsy)
