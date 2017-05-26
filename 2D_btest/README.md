# 2D_btest

2D benchmark test example following numerical example 4 in de Frutos et al. (2014). To run the test example, follow steps 1-3:

1) Change the path in cfg\_datadir.py to your current working directory
2) Create the required mesh geometries by running run\_gmshfiles.sh
3) Start the model by running run\_2Dbtest.sh

The default settings can be found and also changed in cfg.py and cfg\_datadir.py. The run script performs 21 advection time steps, for the default settings of dt = 0.1*pi, 20 advection time steps give one full rotation of 2pi. Results are saved to concentration\_t%d.pvd (t%d - respective advection time step) for visualisation (first advection time step is t1 = 0).
