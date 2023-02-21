cmake_minimum_required(VERSION 3.10)
project(fast-marching-method-test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing(true)

add_subdirectory(gtest-1.7.0)

# Default to release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CompilerFlags
          CMAKE_CXX_FLAGS
          CMAKE_CXX_FLAGS_DEBUG
          CMAKE_CXX_FLAGS_RELEASE
          CMAKE_C_FLAGS
          CMAKE_C_FLAGS_DEBUG
          CMAKE_C_FLAGS_RELEASE)
  foreach(CompilerFlag ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
  endforeach()
  message(STATUS "CXX flags (release): ${CMAKE_CXX_FLAGS_RELEASE}")
  message(STATUS "CXX flags (debug): ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

# Project Compiler Flags
# add_definitions(-Wall)

include_directories(
  ${gtest_SOURCE_DIR}
  ${gtest_SOURCE_DIR}/include)

add_executable(fast-marching-method-test
  main.cpp
  eikonal_solvers_test.cpp
  signed_arrival_time_test.cpp)

target_link_libraries(fast-marching-method-test gtest gtest_main)

add_test(NAME fast-marching-method-test COMMAND fast-marching-method-test)
