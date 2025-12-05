cmake_minimum_required(VERSION 3.10)
project(C_CASOPT VERSION 1.0.0 LANGUAGES CXX)

# --- 1. Default to Debug mode if not specified
# This ensures -g flags (debug symbols) are passed to the compiler
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TARGET_EXE "c_casopt")

set(SOURCES
    src/main.cpp
    src/IO.cpp
)

add_executable(${TARGET_EXE} ${SOURCES})

target_include_directories(${TARGET_EXE} PRIVATE ${PROJECT_SOURCE_DIR})

# Add debug flags explicitly if you want extra info, though CMAKE_BUILD_TYPE=Debug handles most
target_compile_options(${TARGET_EXE} PRIVATE 
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Debug>:-O0 -g3> # Optimization off, max debug info in Debug mode
)

set_target_properties(${TARGET_EXE} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})


find_program(GDB_EXECUTABLE gdb)

if(GDB_EXECUTABLE)
    add_custom_target(debug
        COMMAND ${GDB_EXECUTABLE} -ex run --args $<TARGET_FILE:${TARGET_EXE}>
        DEPENDS ${TARGET_EXE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building and starting GDB for ${TARGET_EXE}..."
        USES_TERMINAL
    )
else()
    message(STATUS "GDB not found. The 'make debug' target will not be available.")
endif()

# --- Add tests (maybe)
