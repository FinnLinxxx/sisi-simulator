cmake_minimum_required(VERSION 3.7)
project(SISI-Simulator)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose type of build." FORCE)
    message(STATUS "Build-Type set ---> Release")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# LTO Check
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.9)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_supported OUTPUT lto_error)
    if(lto_supported)
        message(STATUS "LTO (Link Time Optimization) ENABLED.")
    else()
        message(STATUS "LTO DISABLED: ${lto_error}")
    endif()
else()
    message(STATUS "CMake Version < 3.9: LTO check skipped for compatibility.")
    set(lto_supported FALSE)
endif()

# Suppress Windows/MSVC crying
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Not everything from yaml-cpp is needed:
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

add_subdirectory(external/yaml-cpp)
include_directories(external/yaml-cpp/include)

add_executable(sisi sisi.cpp)
target_link_libraries(sisi yaml-cpp)

if(lto_supported)
    set_target_properties(sisi PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Fast Math Option
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sisi PRIVATE
        -O3 
        -march=native 
        -ffast-math
    )
elseif(MSVC)
    target_compile_options(sisi PRIVATE
        /O2       # Max Speed
        /Ob2      
        /Oi       
        /Ot       
        /fp:fast  
    )
endif()

set(SISI_WORKSPACE_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml")
set(SISI_WORKSPACE_CONFIG_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/config.example.yaml")

if(EXISTS "${SISI_WORKSPACE_CONFIG}")
    message(STATUS "Using existing config.yaml from workspace root.")
elseif(EXISTS "${SISI_WORKSPACE_CONFIG_EXAMPLE}")
    message(WARNING
        "config.yaml was not found in the workspace root. "
        "Creating it from config.example.yaml. "
        "Please review and adjust config.yaml before use."
    )
    configure_file(
        "${SISI_WORKSPACE_CONFIG_EXAMPLE}"
        "${SISI_WORKSPACE_CONFIG}"
        COPYONLY
    )
else()
    message(WARNING
        "Neither config.yaml nor config.example.yaml was found in the workspace root."
    )
endif()

# -lm 
if(UNIX)
    target_link_libraries(sisi m)
endif()
