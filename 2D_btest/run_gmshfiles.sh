#!/bin/bash

## name of geometries
NAME="box"
DIM=2

gmsh ${NAME}.geo -${DIM}
gmsh ${NAME}_fine.geo -${DIM}
gmsh ${NAME}_coarse.geo -${DIM}

dolfin-convert ${NAME}.msh ${NAME}.xml
dolfin-convert ${NAME}_fine.msh ${NAME}_fine.xml
dolfin-convert ${NAME}_coarse.msh ${NAME}_coarse.xml
