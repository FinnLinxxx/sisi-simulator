cmake_minimum_required(VERSION 3.7)
project(SISI-Simulator)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
