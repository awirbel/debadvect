# debadvect
This repository contains a development version of my FEM model for englacial debris transport that has been published in 2018:

 A. Wirbel, A. H. Jarosch, and L. Nicholson, “Modelling debris transport within glaciers by advection in a full-Stokes ice flow model,” The Cryosphere, Vol. 12, Iss. 1, p. 189–204, 2018, doi: [10.5194/tc-12-189-2018](http://dx.doi.org/10.5194/tc-12-189-2018).

A running test example of the benchmark tests in Section 4 can be found in 2D_btest. The model version used in the paper is release v1.0.0.

This model makes use of the FEM ice flow model [icetools](https://github.com/alexjarosch/icetools).

## Requirements
* FEniCS (container coming soon)
* Gmsh
* numpy

A singularity container which fulfills all requirements will be provided soon.
