# The Fast Marching Method
The *Fast Marching Method* (FMM), in its simplest form, can be used to compute the arrival times at grid cells for a monotonously expanding interface. One application of this method is to set the speed of the interface to one, which enables computation of distance fields, where the closest distance to the interface is assigned to every cell in a grid. This repository contains an implementation of the FMM in arbitrary dimensions (actually two or more), although typical usage is limited to 2D and 3D. The code is designed to be simple to incorporate into existing projects, and robustness has been prioritized over speed optimizations. All code in this repository is released under the [MIT license](https://en.wikipedia.org/wiki/MIT_License). If you have any comments or suggestions please feel free to make a pull request.

This note is divided into two major sections. First, we provide examples on how to use the code along with other practical details such as running the accompanying tests. Thereafter, we described the technical choices that were made in the implementation with references to relevant literature.

## Usage
This section describes how to use the FMM implementation provided in this repository. First, we provide some examples on how to call these methods, together with some discussion related to valid inputs. Thereafter, we give instructions on how to run the accompanying tests. Note that there are also a lot of examples within the test code itself, that can be found in the [test folder](https://github.com/thinks/fast-marching-method/tree/master/test). 

The FMM implementation in this repository is contained in a single [header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fast_marching_method/fast_marching_method.hpp). This makes it very easy to add as a dependency to existing projects. Further, the code has no external dependencies other than the standard C++ libraries. All interfaces use standard types, such as `std::array` and `std::vector`. The code contains a fairly large number of `assert` statements, making it easier to debug when the `NDEBUG` preprocessor variable is not defined. However, the code runs very slowly in debug mode so data set sizes may need to be adjusted accordingly.

### Methods
Most of the functions in the single [header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fast_marching_method/fast_marching_method.hpp) (which is all that needs to be included) are in `namespace detail` and should not be called directly. Instead, the is a single entry point provided by the `ArrivalTime` function should be used. As the name suggests this function computes arrival times at grid cells. A conceptual example illustrates what is meant by this.

![alt text](https://github.com/thinks/fast-marching-method/blob/master/img/fmm_readme_concept.png "Conceptual example")

In the figure above, the green circle (_left_) was used as input to compute arrival times on a grid (_right_). Locations outside the circle have positive arrival times (or distances depending on interpretation), here shown in red. Similarly, locations inside the circle have negative arrival times, here shown in blue. The intensity of the colors gives the distance to the interface (i.e. circle boundary). This is why cells close to the interface appear black, since the red or blue component is small. Next, we give an example demonstrating how to write code to generate an image similar to the one shown above.

The input to the `ArrivalTime` is given as grid cells with known distances (or arrival times depending on interpretation). The following code snippet computes a low resolution version of the image shown above.

```cpp
using namespace std;
namespace fmm = thinks::fast_marching_method;

auto circle_boundary_indices = vector<array<int32_t, 2>>{
  {{5, 3}}, {{6, 3}}, {{7, 3}}, {{8, 3}}, {{9, 3}}, {{10, 3}}, {{4, 4}},
  {{5, 4}}, {{10, 4}}, {{11, 4}}, {{3, 5}}, {{4, 5}}, {{11, 5}}, {{12, 5}},
  {{3, 6}}, {{12, 6}}, {{3, 7}}, {{12, 7}}, {{3, 8}}, {{12, 8}}, {{3, 9}},
  {{12, 9}}, {{3, 10}}, {{4, 10}}, {{11, 10}}, {{12, 10}}, {{4, 11}},
  {{5, 11}}, {{10, 11}}, {{11, 11}}, {{5, 12}}, {{6, 12}}, {{7, 12}},
  {{8, 12}}, {{9, 12}}, {{10, 12}},
};
auto circle_boundary_distances = vector<float>{
    0.0417385f, 0.0164635f, 0.0029808f, 0.0029808f, 0.0164635f, 0.0417385f,
    0.0293592f, -0.0111773f, -0.0111773f, 0.0293592f, 0.0417385f, -0.0111773f,
    -0.0111773f, 0.0417385f, 0.0164635f, 0.0164635f, 0.0029808f, 0.0029808f,
    0.0029808f, 0.0029808f, 0.0164635f, 0.0164635f, 0.0417385f, -0.0111773f,
    -0.0111773f, 0.0417385f, 0.0293592f, -0.0111773f, -0.0111773f, 0.0293592f,
    0.0417385f, 0.0164635f, 0.0029808f, 0.0029808f, 0.0164635f, 0.0417385f
};

auto grid_size = array<size_t, 2>{{16, 16}};
auto grid_spacing = array<float, 2>{{1.f/16, 1.f/16}};
auto uniform_speed = 1.f;

auto arrival_times = fmm::SignedArrivalTime(
  grid_size,
  circle_boundary_indices,
  circle_boundary_times,
  fmm::UniformSpeedEikonalSolver<float, 2>(grid_spacing, uniform_speed));
```

First, we define our input, the cell coordinates for which we have known distances. These are stored in two separate lists, one for the coordinates of the cells (`circle_boundary_indices`) and one for the corresponding distances (`circle_boundary_distances`). Normally, these values would of course not be hard-coded like this, but rather generated by some function. Thereafter, we specify the size (`grid_size`) and dimensions (`grid_spacing`) of the grid. Here the grid dimensions are set so that the domain is [0, 1] in each dimension. In order to be able to interpret the arrival times as euclidean distance, a uniform speed of one is set for the entire grid (`uniform_speed`). The speed is passed to an Eikonal solver, which will de discussed in further detail in the following section. The resulting image is shown below.

![alt text](https://github.com/thinks/fast-marching-method/blob/master/img/fmm_readme_input_values.png "Code example")

Cells with known distances are shaded darker grey in the left image. Known input values may be interpreted as radii that intersect the input shape. We note that negative distances are given for cell where the center is inside the circle. In the next section we discuss the use of Eikonal solvers, which allow easy customization of the algorithm while re-using the basic ideas.

### Eikonal Solvers
The basic idea of the FMM algorithm is to propagate given information at known locations to other locations in a numerically reasonable way. 

Eikonal solvers:
* `UniformSpeedEikonalSolver`
* `HighAccuracyUniformSpeedEikonalSolver`
* `VaryingSpeedEikonalSolver`
* `HighAccuracyVaryingSpeedEikonalSolver`


### Input Validation


### Tests
In order to run the tests you need to have [CMake](https://cmake.org/) installed. The tests are implemented in the [Google Test](https://github.com/google/googletest) framework, which is included as part of this repository. 

Running the tests is simple. In a terminal do the following:

```bash
$ cd d:
$ git clone git@github.com:/thinks/fast-marching-method.git D:/fmm
$ mkdir fmm-build
$ cd fmm-build
$ cmake ../fmm/test -DCMAKE_BUILD_TYPE=Release
$ cmake --build . 
$ ctest
```

In order, the following is being done:
* Clone the source code to a directory `D:/fmm`.
* Create an out-of-source build directory `fmm-build`.
* Create the default project files for your machine in the build directory (change `Release` to `Debug` for a debug build).
* Builds the tests.
* Runs the tests. 

If the tests pass you should see something like:

```
Test project D:/fmm-build
    Start 1: fast-marching-method-test
1/1 Test #1: fast-marching-method-test ........     Passed      0.33 sec

100% tests passed, 0 tests failed out of 1

Total test time (real) =    0.35 sec

```

For more detailed test output you can run the test executable directly:

```
$ D:/fmm-build/fast-marching-method-test.exe
```

## Implementation
This section describes the implementation of the FMM implementation from a more technical perspective. 

### Code Design




google coding guidelines

references for further reading

### Future Work
* Termination criteria for narrow band marching.
* Comparison with fast sweeping method and vector distance.


### References
**[1]** J.A. Sethian. A fast marching level set method for monotonically advancing fronts. *Proceeding of the National Academy of Sciences of the USA - Paper Edition*, 93(4):1591-1595, 1996.

**[2]** J. Rickett and S. Fomel. Short note: A second-order fast marching eikonal solver. *Technical Report, Stanford Exploration Project*, 2000.

**[3]** J.A. Baerentzen. On the implementation of fast marching methods for 3D lattices. *Technical Report. IMM-TR-2001-13*, 2001.

**[4]** M.W. Jones, J.A. Baerentzen, and M. Sramek. 3D Distance Fields: A Survey of Techniques and Applications. *IEEE Transactions on Visualization and Computer Graphics*, 12(4):581-599, July/August 2006.

Bridson book















