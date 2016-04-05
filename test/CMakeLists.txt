CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(fast-marching-method)

# Default to release build
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Project Compiler Flags
# ADD_DEFINITIONS(-Wall)

INCLUDE_DIRECTORIES(
  ../../ppm-io/include
  ../include
  include)

ADD_EXECUTABLE(fmm
  main.cpp)
