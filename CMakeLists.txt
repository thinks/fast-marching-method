cmake_minimum_required(VERSION 3.10)
project(fast-marching-method)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Default to release build
if(NOT CMAKE_BUILD_TYPE)
  #set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
  #set(CMAKE_BUILD_TYPE Debug)
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

# Warnings Compiler Flags
set(MyWarningFlags
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Wno-gnu-zero-variadic-macro-arguments>
)
# message(STATUS "Warning flags: ${MyWarningFlags}")

if(MSVC)
  string( REPLACE "/DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
else()
  string( REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()


enable_testing()
add_subdirectory(tests)
add_subdirectory(examples/cpp)
add_subdirectory(bindings/python)