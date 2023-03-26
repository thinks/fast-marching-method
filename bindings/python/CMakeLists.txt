include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.10.4
)
FetchContent_MakeAvailable(pybind11)

# A linked issue is generated at the moment
# See https://github.com/pybind/pybind11/pull/4301
pybind11_add_module(py_fast_marching_method py-fast-marching-method.cpp)