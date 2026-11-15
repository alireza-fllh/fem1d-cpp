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
    test_convergence.cpp)

target_link_libraries(fem1d_tests PRIVATE
    fem1d::fem1d
    Catch2::Catch2WithMain
    fem1d_project_options)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(fem1d_tests)
