cmake_minimum_required(VERSION 3.12)
project(cuda_finite_difference CUDA CXX)

# Create CUDA library
add_library(${PROJECT_NAME}
    FiniteDifference2D.cu
    FiniteDifference2D.cuh
    DeviceMemory.cuh 
    MemoryManagerFiniteDifference2D.cu
    MemoryManagerFiniteDifference2D.cuh 
    RunFiniteDifference2D.h 
    RunFiniteDifference2D.cu 
)

# Set include directories
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
          ${CUDA_INCLUDE_DIRS}

)


target_link_libraries(${PROJECT_NAME} PRIVATE
    
    ${CUDA_LIBRARIES}
    cuda
    cudart
    ${NVML_LIBRARY}
)


# Set CUDA specific properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_RESOLVE_DEVICE_SYMBOLS ON 
    CUDA_ARCHITECTURES "70;75;86" 
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

# Export the target
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_ARCHITECTURES native
)