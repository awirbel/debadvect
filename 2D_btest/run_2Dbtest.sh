#!/bin/bash
set -e


python cfg.py

python makeFiles.py

python read_inputfiles.py

python refine_gl2D.py
python advect_gl.py

for i in {1..20}
do
    python refine_gl2D.py
    python advect_gl.py

done

