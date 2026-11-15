add_executable(fem1d_bench bench_solvers.cpp)
target_link_libraries(fem1d_bench PRIVATE fem1d_cli fem1d_project_options)
target_compile_definitions(fem1d_bench PRIVATE FEM1D_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

if(FEM1D_BUILD_TESTS)
    # Only checks that the benchmark runs; timings from CI machines mean little.
    add_test(NAME benchmark_smoke COMMAND fem1d_bench --quick)
endif()
