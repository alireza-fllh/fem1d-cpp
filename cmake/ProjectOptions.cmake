# Warnings and sanitizers shared by every target in this project.
# Link it PRIVATE so the flags never leak to consumers of the library.

add_library(fem1d_project_options INTERFACE)

set(_gnu_like "$<CXX_COMPILER_ID:GNU,Clang,AppleClang>")
set(_msvc "$<CXX_COMPILER_ID:MSVC>")

target_compile_options(fem1d_project_options INTERFACE
    $<${_gnu_like}:-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
                   -Wold-style-cast -Wnon-virtual-dtor -Woverloaded-virtual>
    $<${_msvc}:/W4 /permissive->)

if(FEM1D_WARNINGS_AS_ERRORS)
    target_compile_options(fem1d_project_options INTERFACE
        $<${_gnu_like}:-Werror>
        $<${_msvc}:/WX>)
endif()

if(FEM1D_ENABLE_SANITIZERS)
    if(MSVC)
        message(FATAL_ERROR "FEM1D_ENABLE_SANITIZERS is only supported with GCC and Clang")
    endif()
    set(_sanitize_flags -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_compile_options(fem1d_project_options INTERFACE ${_sanitize_flags})
    target_link_options(fem1d_project_options INTERFACE ${_sanitize_flags})
endif()
