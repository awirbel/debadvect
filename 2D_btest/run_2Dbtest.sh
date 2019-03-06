#!/bin/bash
set -e

singularity exec fenics_icetools.simg python cfg.py

singularity exec fenics_icetools.simg python makeFiles.py

singularity exec fenics_icetools.simg python read_inputfiles.py

singularity exec fenics_icetools.simg python refine_gl2D.py
singularity exec fenics_icetools.simg python advect_gl.py

for i in {1..20}
do
    singularity exec fenics_icetools.simg python refine_gl2D.py
    singularity exec fenics_icetools.simg python advect_gl.py

done
