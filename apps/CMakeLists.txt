# The command-line front end. Its logic is a library so the tests can call it.
add_library(fem1d_cli
    cli/cases.cpp
    cli/cli.cpp)
target_include_directories(fem1d_cli PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fem1d_cli
    PUBLIC fem1d::fem1d
    PRIVATE fem1d_project_options)

add_executable(fem1d_app main.cpp)
set_target_properties(fem1d_app PROPERTIES OUTPUT_NAME fem1d)
target_link_libraries(fem1d_app PRIVATE fem1d_cli fem1d_project_options)
