# Mixed precision wave propagation to leverage AI hardware in seismic modelling

Hardware is increasingly being designed to support AI/ML applications that make use of reduced precision arithmetic with FP16 (or FP8) as opposed to the standard FP32/FP64.
Seismic modelling and imaging is a major HPC application and requires efficient leveraging of hardware, so we aim to scale problems such that they remain in the FP16 range.
Scaling problems so they remain in the desired FP16 range, without producing any over/underflow will allow for effective use of new hardware systems and a gain in both speed and memory usage.
This work is done primarily using the Devito DSL, a Python package used for solving PDEs with finite difference methods. Users outline a problem in Python and Devito generates C code to execute the simulation.

This repository stores all of the work done in the process of finding a generalised method for scaling equations to FP16 range, thus allowing them to be solved on new AI hardware.

## Contents
- [`fd_coefficients`](./fd_coefficients.py) is a Python script that, given a derivative order, expansion point and set of values, calculates the coefficients for the corresponding finite difference stencil. The algorithm is based on work in the PhD thesis of Ed Caunt. In Devito, the stencil coefficients are generated automatically, but this script will be used when using the mpmath package to hard-code a half-precision implementation
- [`factor_scaling_acousitc.ipynb`](./factor_scaling_acousitc.ipynb) is a Jupyter notebook that details an attempt to apply the scaling method outlined in the 2020 work 'Seismic modeling and inversion using half-precision floating-point numbers' by Gabriel Fabien-Ouellet to the 2d acoustic wave equation. This was unsuccessful, with the conclusion being that the method used in Fabien-Ouellet 2020 requires coupling which the 2d acoustic equation does not have
- [`unit_scaling_acoustic.ipynb`](./unit_scaling_acoustic.ipynb) is a Jupyter notebook that details a successful attempt to scale the 2d acoustic wave equation by defining units for the physical parameters of the problem such that they are between 0 and 1. This keeps every number in the simulation within the FP16 range. The method was done initially with a constant velocity model, before generating an arbitrary 4-layer velocity model and applying maximum, mean and weighted average velocity scaling schemes, all of which were successful
- [`scaling-fletcher.ipynb`](./scaling-fletcher.ipynb) is a Jupyter notebook that details the process of applying the unit scaling method to the coupled TTI equations defined in Fletcher et al 2009. We first drop the equations from 3d to 2d before applying the same method as defined in unit_scaling_acoustic.ipynb
- [`kahan_summation.py`](./kahan_summation.py) is a Python script that will implement the Kahan summation, a method of preserving accuracy when performing floating-point operations. This is a stretch goal for the project so, as yet, is unfinished
- [`Figures`](./Figures) is a folder that contains various figures used in notebooks in this repository 
