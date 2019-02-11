cmake_minimum_required(VERSION 3.5)
project(poseRegularizationToolbox)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 11)

if (WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories(include)

add_library(poseRegularizationToolbox src/pose_regularization.cpp)

# test

option(BUILD_TESTS "Perform unit tests" OFF)

if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(test/ext/gtest)

    set(PROJECT_TEST_NAME run_test)

    # for MSVC ¯\_(ツ)_/¯ - https://github.com/xmfbit/gtest-VisualStudio-CMake-example/blob/master/CMakeLists.txt
    if (MSVC)
        set(variables
                CMAKE_CXX_FLAGS_DEBUG
                CMAKE_CXX_FLAGS_RELEASE
                CMAKE_CXX_FLAGS_RELWITHDEBINFO
                CMAKE_CXX_FLAGS_MINSIZEREL
                )
        foreach (variable ${variables})
            if (${variable} MATCHES "/MD")
                string(REGEX REPLACE "/MD" "/MT" ${variable} "${${variable}}")
            endif ()
        endforeach ()
    endif ()

    find_package(MATIO REQUIRED)

    file(GLOB TEST_SRC_FILES test/*.cpp)
    add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
    add_dependencies(${PROJECT_TEST_NAME} googletest)

    include_directories(BEFORE ${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})

    add_test(test1 ${PROJECT_TEST_NAME})
    target_include_directories(${PROJECT_TEST_NAME} PRIVATE ${MATIO_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_TEST_NAME} poseRegularizationToolbox ${GTEST_LIBRARIES} ${MATIO_LIBRARIES})
    set_property(TARGET ${PROJECT_TEST_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test")
    if (NOT WIN32)
        target_link_libraries(${PROJECT_TEST_NAME} pthread)
    endif (NOT WIN32)

endif (BUILD_TESTS)

