cmake_minimum_required(VERSION 2.8.8)
project(gtest_builder C CXX)
include(ExternalProject)

set(GTEST_DISABLE_PTHREADS OFF)

if (WIN32)
    set(GTEST_DISABLE_PTHREADS ON)
endif (WIN32)

ExternalProject_Add(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        CMAKE_ARGS
        -Dgtest_disable_pthreads=${GTEST_DISABLE_PTHREADS}
        -DBUILD_GTEST=ON
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
        # Disable install step
        INSTALL_COMMAND ""
        )

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include PARENT_SCOPE)
set(GMOCK_INCLUDE_DIRS ${source_dir}/googlemock/include PARENT_SCOPE)

# Specify MainTest's link libraries
ExternalProject_Get_Property(googletest binary_dir)

if (WIN32)
    set(GTEST_LIBRARY ${binary_dir}/lib/${CMAKE_SHARED_MODULE_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(GTEST_MAIN_LIBRARY ${binary_dir}/lib/${CMAKE_SHARED_MODULE_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX})
else (WIN32)
    set(GTEST_LIBRARY ${binary_dir}/lib/${CMAKE_SHARED_MODULE_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(GTEST_MAIN_LIBRARY ${binary_dir}/lib/${CMAKE_SHARED_MODULE_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX})
endif (WIN32)


set(GTEST_LIBRARIES ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} PARENT_SCOPE)
