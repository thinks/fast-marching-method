# The Fast Marching Method
The *Fast Marching Method* (FMM) is a method for computing distance fields. The FMM was introduced to the level set community by Sethian in 1996. An in-depth theoretical discussion of the mathematical foundations of the method is beyond the scope of this note, those interested in such details are referred to the vast body of literature on the subject. This repository contains an implementation of the FMM in arbitrary dimensions, although typical usage is most likely limited to 2D and 3D. The code is designed to be simple to incorporate into existing projects and robustness has been prioritized over speed. A number of tests have been performed and the results are given in both visual and numeric versions below. 

The basic idea of the FMM is that distance is propagated from known values to more distant locations. Computations are done on a grid providing simple topology for computing gradients. Fast marching methods have also been adapted to more general manifold structures but such methods are not implemented in this repository. Grid cells are classified to be in one of three states: (1) *frozen*; (2) *narrow band*; or *far*. In this context, *frozen* means that the distance value in that cell has its final value. A list of frozen cells is provided as input to the FMM functions. The distances in these cells are the boundary conditions from which the rest of the grid is solved. As a first step, the *face neighbors* (4-connected in 2D) of all frozen cells are marked as *narrow band* (if they are not already *frozen*!) and tentative distances are estimated for those cells. Since distance information should propagate outwards from frozen cells, the *narrow band* cell with the smallest distance is then marked as *frozen* and its neighbors (if not already *frozen*!) are added to the list of *narrow band* cells. Finding the *narrow band* cell with the smallest distance can be time consuming, so a min-heap is used to store information about the *narrow band* cells. This allows *log(N)* access to the cell with the smallest distance value, where *N* is the number of cells currently in the *narrow band*. The method continues until all cells have been frozen, at which point no further updates are possible.

what is the fast marching method? 
what does it solve? 
why is it useful?
references for further reading

## Design
The FMM implementation in this repository is in a [single header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fastMarchingMethod.hpp). This makes it very easy to add as a dependency to existing projects. Further, the code has no external dependencies other than the standard C++ libraries. All interfaces use standard types, such as `std::array` and `std::vector`. As mentioned earlier robustness and correctness have been prioritized over speed. The code contains a fairly large number of `assert` statements, making it easier to debug when the `NDEBUG` preprocessor variable is not defined. However, the code runs very slowly in debug mode so test data set sizes might need to be adjusted accordingly.

Besides the FMM implmentation, a some tests are provided in this repository in the [test folder](https://github.com/thinks/fast-marching-method/tree/master/test). The test code is of course not required to use the implementation, but demonstrates the correctness of the implementation. Also, if changes are made in a cloned repository extending the tests to cover new features or bug fixes becomes possible. If you do find any improvements feel free to make a pull request! An example file running the tests is included [here](https://github.com/thinks/fast-marching-method/blob/master/test/main.cpp) together with a very simple [CMake project](https://github.com/thinks/fast-marching-method/blob/master/test/CMakeLists.txt). The results of the testing that has currently been performed are discussed below.

## Usage
The easiest way to see how to use the FMM distance functions in this repository is to have a look at the accompanying tests in [this](https://github.com/thinks/fast-marching-method/blob/master/test/include/thinks/testFastMarchingMethod.hpp) header file. Those test functions set up valid inputs to the FMM functions and analyze the output.



## Tests
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


## Future Work
* Higher order spatial gradients.
* Non-uniform speed function.
* Specialize Eikonal solver for 2d/3d.

