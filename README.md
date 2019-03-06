# debadvect
This repository contains a development version of my FEM model for englacial debris transport that has been published in 2018:

 A. Wirbel, A. H. Jarosch, and L. Nicholson, “Modelling debris transport within glaciers by advection in a full-Stokes ice flow model,” The Cryosphere, Vol. 12, Iss. 1, p. 189–204, 2018, doi: [10.5194/tc-12-189-2018](http://dx.doi.org/10.5194/tc-12-189-2018).

A running test example of the benchmark tests in Section 4 can be found in `2D_btest`. The model version used in the paper is release v1.0.0.

This model makes use of the FEM ice flow model [icetools](https://github.com/alexjarosch/icetools).

## Requirements
* [singularity v2.6](https://www.sylabs.io/guides/2.6/user-guide/installation.html)
* [Gmsh](http://gmsh.info/)

A singularity container with FEniCS v2016.2 installed is required to run the example in `2D_btest`.

After you have cloned the repo to your machine with
```shell
git clone https://github.com/awirbel/debadvect.git
```
you can download the singularity container form icetools release v1.2.0 called [fenics_icetools.simg](https://github.com/alexjarosch/icetools/releases/download/v1.2.0/fenics_icetools.simg) (approx 1.5 GB).
Place this container in the `2D_btest` directory.

Given that your singularity and gmsh installs are working, you should be able to run the `2D_btest` case as described in its readme file.
