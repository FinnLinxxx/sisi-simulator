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

# Only for Windows/Visual Studio: Copys config to .exe path
if(MSVC)
    add_custom_command(TARGET sisi POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
        "$<TARGET_FILE_DIR:sisi>/config.yaml"
        COMMENT "Windows: Copy the configuration to the Build folder..."
    )
endif()

# -lm 
if(UNIX)
    target_link_libraries(sisi m)
endif()
