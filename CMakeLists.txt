cmake_minimum_required(VERSION 3.5)
project(poseRegularizationToolbox)

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DSE3PROD)

include_directories(include)

set(SOURCE_FILES
        src/manifoldGeneric.cpp
        src/manifoldSE3prod.cpp)

add_library(poseRegularizationToolbox ${SOURCE_FILES})

if(WITH_PYTHON)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)
    include(UseCython)
    
    add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
        ${CMAKE_CURRENT_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/wrapper)
        
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/wrapper/pose_denoising.pyx
        PROPERTIES CYTHON_IS_CXX TRUE)
    
    cython_add_module(pose_denoising python/wrapper/pose_denoising.pyx src/manifoldGeneric.cpp src/manifoldSE3prod.cpp)

endif(WITH_PYTHON)

# test

option(DO_TESTS "Perform unit tests" OFF)

if(DO_TESTS)
    enable_testing()
    add_subdirectory(test/ext/gtest)

    set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)

    file(GLOB TEST_SRC_FILES test/*.cpp)
    add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
    add_dependencies(${PROJECT_TEST_NAME} googletest)

    include_directories(${GTEST_INCLUDE_DIRS})

    add_test(test1 ${PROJECT_TEST_NAME})
endif(DO_TESTS)

