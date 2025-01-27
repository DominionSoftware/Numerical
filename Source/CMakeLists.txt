cmake_minimum_required(VERSION 3.12)

set(CMAKE_CUDA_ARCHITECTURES 70;75;86)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
if(UNIX AND NOT APPLE)
    set(CMAKE_C_COMPILER /usr/bin/gcc-11)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-11)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fconcepts-ts)
endif()

set(PROJECT_NAME Numerical)
project(${PROJECT_NAME} LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 20)


if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(FATAL_ERROR "GCC version must be at least 11 for C++20 support")
endif()

if(WIN32)
    option(USE_SCIPLOT "Enable SciPlot integration" ON)
endif(WIN32)



# Look for CUDA 12.6
find_package(CUDA 12.6...12.8)

if(CUDA_FOUND)
    # Enable CUDA support for the project
    enable_language(CUDA)
    message(STATUS "CUDA 12.6 found. Adding CUDA support.")
    
    # Set CUDA specific flags
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
 
    # Add the CUDA subproject
    add_subdirectory(CudaFiniteDifference)
else()
    message(STATUS "CUDA 12.6 not found. Skipping CUDA components.")
endif()

# Add sciplot dependency
if(USE_SCIPLOT)
    if(WIN32)
        set(SCIPLOT_DIR "D:/Projects/sciplot.ins/share/sciplot")
    endif(WIN32)
    
    find_package(SCIPLOT)
    if(SCIPLOT_FOUND)
        set(SCIPLOT_SOURCE "${SCIPLOT_DIR}/../../include")
        message(STATUS "SciPlot found and enabled")
    else()
        message(WARNING "SciPlot not found. Disabling SciPlot integration.")
        set(USE_SCIPLOT OFF)
    endif()
endif()
# Add Eigen dependency
if(WIN32)
    set(Eigen3_DIR "D:/Projects/eigen.ins/share/eigen3/cmake")
endif(WIN32)
find_package(Eigen3 REQUIRED NO_MODULE)

# Add Google Test as a subdirectory
if(WIN32)
    set(CMAKE_PREFIX_PATH "D:/Projects/googletest.ins" ${CMAKE_PREFIX_PATH})
endif(WIN32)
find_package(GTest REQUIRED)
if(WIN32)
    set(GTEST_INCLUDES "D:/Projects/googletest.ins/include")
endif(WIN32)
# Source files
set(SOURCE
    main.cpp
    Linspace.h
    FiniteDifference.h
    Integral.h
    Roots.h 
    Abscissae.h 
    AbscissaeFromData.h 
    WriteVector.h 
    valarray2vector.h
    GenerateLegendreCoeffs.h 
    GaussKronrodConstants.h
    GaussKronrodRules.h 
    ClenshawCurtis.h 
)

# Add the main executable
add_executable(${PROJECT_NAME} ${SOURCE})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCIPLOT_SOURCE} ${EIGEN3_INCLUDE_DIRS})


if(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

if(CUDA_FOUND)
    # Link CUDA components if available
    target_link_libraries(${PROJECT_NAME} PRIVATE cuda_finite_difference)
endif()

# Create a separate target for tests
add_executable(${PROJECT_NAME}_tests
    TestMain.cpp
    TestFiniteDifference.cpp
    TestIntegration.cpp
    TestRootFinding.cpp
    TestFiniteDifference2D.cpp
)
if(USE_SCIPLOT)
    target_compile_definitions(${PROJECT_NAME}_tests PRIVATE USE_SCIPLOT)
endif()
# Link the test executable to Google Test and other libraries
target_link_libraries(${PROJECT_NAME}_tests PRIVATE GTest::GTest GTest::Main ${EIGEN3_LIBRARIES})

if(CUDA_FOUND)
    # Link CUDA components to tests if available
    target_link_libraries(${PROJECT_NAME}_tests PRIVATE cuda_finite_difference)
endif()

# Link the sciplot library if needed in tests
target_include_directories(${PROJECT_NAME}_tests PRIVATE ${SCIPLOT_SOURCE} ${GTEST_INCLUDES} ${EIGEN3_INCLUDE_DIR})

# Enable testing
enable_testing()
add_test(NAME ${PROJECT_NAME}_tests COMMAND ${PROJECT_NAME}_tests)

# Create a separate target for generation of data
add_executable(${PROJECT_NAME}_data_generation GenerateData.cpp)

# Link the test executable to Google Test and other libraries
target_link_libraries(${PROJECT_NAME}_data_generation PRIVATE ${EIGEN3_LIBRARIES})

# Link the sciplot library if needed in tests
target_include_directories(${PROJECT_NAME}_data_generation PRIVATE ${EIGEN3_INCLUDE_DIR})


# Dump all variables for debugging
get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
