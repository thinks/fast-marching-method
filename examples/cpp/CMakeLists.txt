add_executable(fast-marching-minimal-example
  minimal_example.cpp)

# Warnings Compiler Flags
target_compile_options(fast-marching-minimal-example PRIVATE ${MyWarningFlags})