include(FetchContent)

FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
    GIT_SHALLOW    TRUE
    SYSTEM)
FetchContent_MakeAvailable(Catch2)

add_executable(fem1d_tests
    test_mesh.cpp
    test_element.cpp
    test_linear_system.cpp
    test_assembler.cpp
    test_linear_solver.cpp
    test_boundary_conditions.cpp
    test_bar_scenarios.cpp
    test_convergence.cpp
    test_io.cpp
    test_golden.cpp
    test_cli.cpp)

target_link_libraries(fem1d_tests PRIVATE
    fem1d::fem1d
    fem1d_cli
    Catch2::Catch2WithMain
    fem1d_project_options)

target_compile_definitions(fem1d_tests PRIVATE
    FEM1D_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(fem1d_tests)
