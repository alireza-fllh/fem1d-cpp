cmake_minimum_required(VERSION 3.25)

project(fem1d
    VERSION 0.1.0
    DESCRIPTION "Linear finite elements for the 1D elastic bar -(EA u')' = f"
    LANGUAGES CXX)

option(FEM1D_BUILD_TESTS "Build the unit tests" ON)
option(FEM1D_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(FEM1D_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/ProjectOptions.cmake)

add_library(fem1d
    src/mesh.cpp
    src/element.cpp
    src/linear_system.cpp
    src/assembler.cpp
    src/boundary_conditions.cpp
    src/linear_solver.cpp
    src/bar_solver.cpp)
add_library(fem1d::fem1d ALIAS fem1d)

target_include_directories(fem1d PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(fem1d PUBLIC cxx_std_17)
set_target_properties(fem1d PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(fem1d PRIVATE fem1d_project_options)

if(FEM1D_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
