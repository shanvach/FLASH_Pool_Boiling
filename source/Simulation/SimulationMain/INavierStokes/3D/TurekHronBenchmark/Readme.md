# Turek and Hron Benchmark Problem

## Overview
This problem is defined by [Turek and Hron (2009)][Turek2009] as a flexible beam mounted to a rigid cylinder.  The problem is originally 2D, however it is being tested in FLASH using the 3D solvers and only a few cells in the span direction.  Homogeneous boundary conditions are also applied in the span to ensure that the flow and structure maintain planar uniformity.

## System Parameters
The parameters for the system are taken from [Turek and Hron (2009)][Turek2009] Table 12.  Recapitulated here:

|                    |      Parameter     | FSI-1 | FSI-2 | FSI-3 |
|--------------------|:------------------:|:-----:|:-----:|:-----:|
| Density Ratio      |  $\rho_s / \rho_f$ |     1 |    10 |     1 |
| Reynolds Number    |   $Re=U d / \nu$   |    20 |   100 |   200 |
| Relative Stiffness | $E/( \rho_f U^2 )$ | 3.5e4 | 1.4e3 | 1.4e3 |
| Poisson's Ratio    |       $\nu_s$      |   0.4 |   0.4 |   0.4 |

These are similar to the polyethylene in glycerine configuration discussed early in the text.  However, for convenience these "nicer looking" parameters are used instead.

The geometry in the FLASH simulation has been scaled relative to the diameter of the cylinder.  These new values are given in the next Table.

| Non-dimensional Geometry quantity (relative to diameter) | Symbol |  Value |
|----------------------------------------------------------|:------:|-------:|
| Cylinder diameter                                        |   $d$  |      1 |
| Channel length                                           |   $L$  |     25 |
| Channel height                                           |   $H$  |    4.1 |
| Center of cylinder                                       |   $C$  | (2, 2) |
| Cylinder radius                                          |   $r$  |    0.5 |
| Beam length                                              |   $l$  |    3.2 |
| Beam thickness                                           |   $h$  |    0.2 |
| Reference point                                          |   $A$  | (6, 2) |


## The Fluid Grid
The problem is currently setup to use the uniform grid solver.  The boundary conditions are a parabola for the inflow, no-slip channel walls on the top and bottom, outflow on the exit, and periodic(?) on the normal.

The starting guess of the cell size is 1/50 the diameter of the cylinder.


## The structural body
The body is drawn in Gmsh and then processed in Matlab.  There are three body files in this directory, the prefix FSI-%d of each corresponds to the parameters from the data above.  




[Turek2009]: http://dx.doi.org/10.1007/3-540-34596-5_15

