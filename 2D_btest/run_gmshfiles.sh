#!/bin/bash

## name of geometries
NAME="box"
DIM=2

gmsh ${NAME}.geo -${DIM}
gmsh ${NAME}_fine.geo -${DIM}
gmsh ${NAME}_coarse.geo -${DIM}

singularity exec fenics_icetools.simg dolfin-convert ${NAME}.msh ${NAME}.xml
singularity exec fenics_icetools.simg dolfin-convert ${NAME}_fine.msh ${NAME}_fine.xml
singularity exec fenics_icetools.simg dolfin-convert ${NAME}_coarse.msh ${NAME}_coarse.xml
