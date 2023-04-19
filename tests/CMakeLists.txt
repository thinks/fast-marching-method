include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        750d67d809700ae8fca6d610f7b41b71aa161808
  SYSTEM
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")
set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
set_target_properties(gmock PROPERTIES CXX_CLANG_TIDY "")
set_target_properties(gmock_main PROPERTIES CXX_CLANG_TIDY "")

add_executable(fast-marching-method-test
  main.cpp
  eikonal_solvers_test.cpp
  signed_arrival_time_test.cpp)

target_link_libraries(fast-marching-method-test PRIVATE GTest::gtest GTest::gtest_main)

# Warnings Compiler Flags
target_compile_options(fast-marching-method-test PRIVATE ${MyWarningFlags})

add_test(NAME fast-marching-method-test COMMAND fast-marching-method-test)

# Python bindings test
add_test(NAME py-bindings-test
    COMMAND ${PYTHON_EXECUTABLE} py-bindings-test.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    
set_tests_properties(py-bindings-test PROPERTIES
    ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${CMAKE_BINARY_DIR}/bindings/python")