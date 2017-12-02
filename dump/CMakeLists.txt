CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(fast-marching-method-dump)

# Default to release build
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

IF(MSVC)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  SET(CompilerFlags
          CMAKE_CXX_FLAGS
          CMAKE_CXX_FLAGS_DEBUG
          CMAKE_CXX_FLAGS_RELEASE
          CMAKE_C_FLAGS
          CMAKE_C_FLAGS_DEBUG
          CMAKE_C_FLAGS_RELEASE)
  FOREACH(CompilerFlag ${CompilerFlags})
    STRING(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
  ENDFOREACH()
  MESSAGE(STATUS "CXX flags (release): ${CMAKE_CXX_FLAGS_RELEASE}")
  MESSAGE(STATUS "CXX flags (debug): ${CMAKE_CXX_FLAGS_DEBUG}")
ENDIF()

# Project Compiler Flags
# ADD_DEFINITIONS(-Wall)

INCLUDE_DIRECTORIES(
  ../../ppm-io/include)

ADD_EXECUTABLE(fast-marching-method-dump
  main.cpp)
