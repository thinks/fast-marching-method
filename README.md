# The Fast Marching Method
The *Fast Marching Method* (FMM) is a method for computing distance fields. The FMM was introduced to the level set community by Sethian in 1996. We will not give a complete theoretical discussion of the mathematical foundations of the method here, those interested in such details are referred to the vast body of literature on the subject. This repository contains an implementation of the FMM in arbitrary dimensions, although typical usage is most likely limited to 2D and 3D. The code was designed to be simple to incorporate into existing projects and robustness has been prioritized over speed. A number of tests have been performed and the results are given in both visual and numeric versions below. 

what is the fast marching method? 
what does it solve? 
why is it useful?
references for further reading

## Design


single header file
prefer robustness and correctness over speed.
numerical and visual testing

## Testing and Verification
![hello world](https://github.com/thinks/fast-marching-method/blob/master/test/img/unsigned_grad_mag_float.png?raw=true "grad mag")


## Future Work
* Higher order spatial gradients.
* Specialize Eikonal solver for 2d/3d.
