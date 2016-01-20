CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(fast-marching-method)

# Default to release build
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Project Compiler Flags
#ADD_DEFINITIONS(-Wall -Wextra)

SET(fmm_HEADERS
)

SET(fmm_SOURCES
  main.cpp
)

INCLUDE_DIRECTORIES(
  ../include
  include)

ADD_EXECUTABLE(fmm
  ${fmm_SOURCES}
  ${fmm_HEADERS})
