# The Fast Marching Method
The *Fast Marching Method* (FMM), in its simplest form, can be used to compute the arrival times at grid points for a monotonously expanding interface. One application of this method is to compute distance fields, where the closest distance to an interface is computed at every cell in a grid. This repository contains an implementation of the FMM in arbitrary dimensions (actually two or more), although typical usage is limited to 2D and 3D. The code is designed to be simple to incorporate into existing projects, and robustness has been prioritized over speed optimizations. All code in this repository is released under the [MIT license](https://en.wikipedia.org/wiki/MIT_License). If you have any comments or suggestions please feel free to make a pull request.

This note is divided into two major sections. First, we provide examples on how to use the code, along with other practical details such as running the accompanying tests. Thereafter, we described the technical choices that were made in the implementation with references to relevant literature.

## Usage
This section describes how to use the FMM functions provided in this repository. First, we provide some examples on how to call these methods, together with some discussion related to valid inputs. Thereafter, we give instructions on how to run the accompanying tests. Note that there are also a lot of examples within the test code itself, which can be found in the [test folder](https://github.com/thinks/fast-marching-method/tree/master/test). 

The FMM implementation in this repository is contained in a single [header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fast_marching_method/fast_marching_method.hpp). This makes it very easy to add as a dependency to existing projects. Further, the code has no external dependencies other than the standard C++ libraries. All interfaces use standard types, such as `std::array` and `std::vector`. The code contains a fairly large number of `assert` statements, making it easier to debug when the `NDEBUG` preprocessor variable is not defined. However, the code runs very slowly in debug mode so data set sizes may need to be adjusted accordingly.

### Methods
The single [header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fast_marching_method/fast_marching_method.hpp) (which is all that needs to be included) contains only two "public" functions. The other functions are in a detail `namespace` and should not be called. The two functions are:
* `SignedArrivalTime`
* `UnsignedArrivalTime`

They both compute arrival times, but in the first case arrival times for cells inside a connected component will be negated. The interfaces of these two functions are very similar, and one could also argue the same for the output. 

![alt text](https://github.com/thinks/fast-marching-method/blob/master/img/input.png "Signed vs unsigned")

In the figure above, the green circle (_left_) is given as input to the two methods. Positive arrival times (or distances depending on interpretation) are shown in red, negative arrival times are shown in blue. In the case of `UnsignedArrivalTime` (_middle_) all arrival times are positive, regardless of being inside or outside the circle. For the `SignedArrivalTime` function (_right_), locations inside the circle have negative distances. Note that the magnitudes of the arrival times are identical for both functions, the only difference is the sign for locations inside the circle. Next, we give an example showing the code used to generate the images discussed in this paragraph.

First of all, the input to the FMM functions is given as grid cells with known distances (or arrival times depending on interpretation). From this boundary condition distances at other grid cells are computed.  

From a more technical point of view, the code required to generate the arrival times for locations 
There are, however, some subtle differences in how they should be called. An example will illustrate the difference.

### Input Validation

### Speed Function

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



### References
[Sethian96] J.A. Sethian A fast marching level set method for monotonically advancing fronts. *Proceeding of the National Academy of Sciences of the USA - Paper Edition*. 93(4):1591-1595, 1996


















