# The Fast Marching Method
The *Fast Marching Method* (FMM) is a method for computing distance fields. The FMM was introduced to the level set community by Sethian in 1996. We will not give a complete theoretical discussion of the mathematical foundations of the method here, those interested in such details are referred to the vast body of literature on the subject. This repository contains an implementation of the FMM in arbitrary dimensions, although typical usage is most likely limited to 2D and 3D. The code was designed to be simple to incorporate into existing projects and robustness has been prioritized over speed. A number of tests have been performed and the results are given in both visual and numeric versions below. 

what is the fast marching method? 
what does it solve? 
why is it useful?
references for further reading

## Design
The FMM implementation in this repository is in a [single header file](https://github.com/thinks/fast-marching-method/blob/master/include/thinks/fastMarchingMethod.hpp). This makes it very easy to add as a dependency to existing projects. Further, the code has no external dependencies other than the standard C++ libraries. All interfaces use standard types, such as `std::array` and `std::vector`. As mentioned earlier robustness and correctness have been prioritized over speed. The code contains a fairly large number of `assert` statements, making it easier to debug when the `NDEBUG` preprocessor variable is not defined. However, the code runs very slowly in debug mode so test data set sizes might need to be adjusted accordingly.

Besides the FMM implmentation, a some tests are provided in this repository in the [test folder](https://github.com/thinks/fast-marching-method/tree/master/test). The test code is of course not required to use the implementation, but demonstrates the correctness of the implementation. Also, if changes are made in a cloned repository extending the tests to cover new features or bug fixes becomes possible. If you do find any improvements feel free to make a pull request! An example file running the tests is included [here](https://github.com/thinks/fast-marching-method/blob/master/test/main.cpp) together with a very simple [CMake project](https://github.com/thinks/fast-marching-method/blob/master/test/CMakeLists.txt). The results of the testing that has currently been performed are discussed below.

## Usage
The easiest way to see how to use the FMM distance functions in this repository is to have a look at the accompanying tests in [this](https://github.com/thinks/fast-marching-method/blob/master/test/include/thinks/testFastMarchingMethod.hpp) header file.



## Testing and Verification
![hello world](https://github.com/thinks/fast-marching-method/blob/master/test/img/unsigned_grad_mag_float.png?raw=true "grad mag")


## Future Work
* Higher order spatial gradients.
* Specialize Eikonal solver for 2d/3d.
