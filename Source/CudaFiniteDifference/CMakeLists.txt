cmake_minimum_required(VERSION 3.12)
project(cuda_finite_difference CUDA)

# Create CUDA library
add_library(${PROJECT_NAME}
    FiniteDifference2D.cu
    FiniteDifference2D.cuh
)

# Set include directories
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set CUDA specific properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

# Export the target
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_ARCHITECTURES native
)