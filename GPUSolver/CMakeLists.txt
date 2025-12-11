cmake_minimum_required(VERSION 3.18)

# 1. Initialize project
project(APU_FEM VERSION 1.0.0 LANGUAGES CXX CUDA)

# --- 2. AGGRESSIVE RPATH CONFIGURATION ---
# We force CMake to write the installation RPATH into the build tree executable.
# We also explicitly initialize the RPATH list.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "") # Initialize empty

# --- Standard Settings ---
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)

set(TARGET_EXE "poissfem")
set(SOURCES src/main.cpp src/IO.cpp src/ComputeModel.cpp)

add_executable(${TARGET_EXE} ${SOURCES})

# 3. Find Dependencies
find_package(CUDAToolkit REQUIRED)

# Find cuDSS (using Spack environment)
find_library(CUDSS_LIB NAMES cudss)
find_path(CUDSS_INCLUDE_DIR NAMES cudss.h)

if(NOT CUDSS_LIB OR NOT CUDSS_INCLUDE_DIR)
    message(FATAL_ERROR "cuDSS not found. Check 'module load cudss'.\nLib: ${CUDSS_LIB}")
endif()

# 4. MANUAL RPATH INJECTION (The Fix)
# We ask CMake where it found cublas and cudss, then add those folders to RPATH.

# Get path to cublas
get_target_property(CUBLAS_LIB_PATH CUDA::cublas LOCATION)
get_filename_component(CUBLAS_LIB_DIR ${CUBLAS_LIB_PATH} DIRECTORY)

# Get path to cudss
get_filename_component(CUDSS_LIB_DIR ${CUDSS_LIB} DIRECTORY)

# Add both to the RPATH of the executable
list(APPEND CMAKE_INSTALL_RPATH "${CUBLAS_LIB_DIR}")
list(APPEND CMAKE_INSTALL_RPATH "${CUDSS_LIB_DIR}")

# Print to verify (check this in your output!)
message(STATUS "injecting RPATH for Runtime: ${CUBLAS_LIB_DIR}")
message(STATUS "injecting RPATH for Runtime: ${CUDSS_LIB_DIR}")

set_target_properties(${TARGET_EXE} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_RPATH}")

# 5. Link Libraries
target_link_libraries(${TARGET_EXE} PRIVATE 
    ${CUDSS_LIB} 
    CUDA::cublas 
    CUDA::cudart
)

# 6. Include Directories
target_include_directories(${TARGET_EXE} PRIVATE 
    ${PROJECT_SOURCE_DIR} 
    ${CUDSS_INCLUDE_DIR}
)

# --- Compile Options (ignoring warnings for now to keep log clean) ---
target_compile_options(${TARGET_EXE} PRIVATE 
    $<$<CONFIG:Debug>:-O0 -g3> 
)

set_target_properties(${TARGET_EXE} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})