cmake_minimum_required(VERSION 3.5)
project(poseRegularizationToolbox)

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DSE3PROD)

include_directories(include)

set(SOURCE_FILES
        src/manifoldGeneric.cpp
        src/manifoldSE3prod.cpp)

add_library(poseRegularizationToolbox ${SOURCE_FILES})