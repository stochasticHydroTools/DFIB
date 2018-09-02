# Divergence-Free Immersed Boundary (DFIB) Method

This repository contains MATLAB demo codes for the  DFIB method, which substantially improves the volume (mass) conservation of the conventional immersed boundary method. Currently, the method is limited to periodic boundary conditions.

For details of the method, see the paper:

**An Immersed Boundary Method with Divergence-Free Velocity Interpolation and Force Spreading**, Y. Bao, A. Donev, B.E. Griffith, D.M. McQueen, and C.S. Peskin, Journal of Computational Physics, Vol 347, 183-206 (2017). [JCP](https://www.sciencedirect.com/science/article/pii/S0021999117304953) [arXiv](https://arxiv.org/abs/1701.07169)

## Table of Contents
* `DFIB_SpreadInterp2D_MEX/`: DFIB spreading and interpolation operators in 2D.
* `DFIB_SpreadInterp3D_MEX/`: DFIB spreading and interpolation operators in 3D.
* `DFIBsolver2D/`: contains codes for solving IB equations in 2D.
* `DFIBsolver3D/`: contains codes for solving IB equations in 3D. 
* `Kernels_MEX/`: contains MEX C codes for spreading/interpolating various IB kernels (Peskin kernels).
* `EXAMPLES/`: contains _main_ programs for demos. 

## Instructions
* Before running any tests in `EXAMPLES/`, first build the MEX files. In MATLAB command line, type:

```
cd Kernels_MEX
mex KernelGrid2D.c Kernels.c
mex KernelGrid3D.c Kernels.c
```
* To run the 3D surface tension test, you also need to build the MEX files:

```
cd EXAMPLES
mex ForceSurfTension.c
mex TetrahedronSphereVolume.c
```

