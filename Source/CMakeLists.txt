cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME Numerical)
project(${PROJECT_NAME} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

# Look for CUDA 12.6
find_package(CUDA 12.6)

if(CUDA_FOUND)
    # Enable CUDA support for the project
    enable_language(CUDA)
    message(STATUS "CUDA 12.6 found. Adding CUDA support.")
    
    # Set CUDA specific flags
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES "70;75;86")
    # Add the CUDA subproject
    add_subdirectory(CUDAFiniteDifference)
else()
    message(STATUS "CUDA 12.6 not found. Skipping CUDA components.")
endif()

# Add sciplot dependency
set(CMAKE_PREFIX_PATH "Y:/sciplot.ins" ${CMAKE_PREFIX_PATH})
set(SCIPLOT_DIR "Y:/Projects/sciplot.ins")
find_package(SCIPLOT REQUIRED)
set(SCIPLOT_SOURCE "${SCIPLOT_DIR}/../../include")

# Add Eigen dependency
set(CMAKE_PREFIX_PATH "Y:/eigen340.ins" ${CMAKE_PREFIX_PATH})
set(Eigen3_DIR "Y:/eigen.ins/share/eigen3/cmake")
find_package(Eigen3 REQUIRED NO_MODULE)

# Add Google Test as a subdirectory
set(CMAKE_PREFIX_PATH "Y:/googletest.ins" ${CMAKE_PREFIX_PATH})
find_package(GTest REQUIRED)
set(GTEST_INCLUDES "Y:/googletest.ins/include")

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
)

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
