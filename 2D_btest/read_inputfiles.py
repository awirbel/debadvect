"""

    This function just reads all the required input files for the 2D three body rotation test

    Copyright 2017 Anna Wirbel
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

"""

from dolfin import *
from cfg_datadir import data_dir
import os
import pickle

# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
data_dir1 = os.path.join(data_dir, 'input_2Dbtest')

# Load mesh from file
mesh_vel = Mesh()
mname = "%s" % (fname)
hdf = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir, mname+".hdf5"), "r")
hdf.read(mesh_vel, mname, False)
hdf.close()
mesh_coarse = Mesh()
m1name = "%s_coarse" % (fname)
hdf1 = HDF5File(mesh_coarse.mpi_comm(), os.path.join(data_dir, m1name+".hdf5"), "r")
hdf1.read(mesh_coarse, m1name, False)
hdf1.close()
mesh_a = Mesh()
m2name = "%s_fine" % (fname)
hdf2 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, m2name+".hdf5"), "r")
hdf2.read(mesh_a, m2name, False)
hdf2.close()

# Define Function Spaces
Q = FunctionSpace(mesh_a, "CG", 1)
V2 = VectorFunctionSpace(mesh_vel, "CG", 2)
V_coarse = VectorFunctionSpace(mesh_coarse, "CG", 2)


# Read concentration and velocities from file
c0 = Function(Q)
c0inname = "concentration_refine_%s" % (fname)
hdf3 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir1, c0inname+".hdf5"), "r")
hdf3.read(c0, c0inname)
hdf3.close()
velocity = Function(V2)
velinname = "velocity_%s" % (fname)
hdf4 = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir1, velinname+".hdf5"), "r")
hdf4.read(velocity, velinname)
hdf4.close()
velocity_coarse = Function(V_coarse)
velin2name = "velocity_coarse_%s" % (fname)
hdf5 = HDF5File(mesh_coarse.mpi_comm(), os.path.join(data_dir1, velin2name+".hdf5"), "r")
hdf5.read(velocity_coarse, velin2name)
hdf5.close()

# Save input parameters to file
params['mass_initial'] = 934.887940398
params['mass_previous'] = 934.887940398
params['velmax'] = 50.0

with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)

# Save to file
mfilename = "refined%s" % (fname)
fm = HDF5File(mesh_a.mpi_comm(),os.path.join(data_dir, mfilename+".hdf5"), 'w')
fm.write(mesh_a, mfilename)

cfilename = "concentration_refine_%s" % (fname)
f1 = HDF5File(mesh_a.mpi_comm(),os.path.join(data_dir, cfilename+".hdf5"), 'w')
f1.write(c0, cfilename)


v2filename = "velocity_coarse_%s" % (fname)
f2 = HDF5File(mesh_coarse.mpi_comm(),os.path.join(data_dir, v2filename+".hdf5"), 'w')
f2.write(velocity_coarse, v2filename)

v3filename = "velocity_%s" % (fname)
f3 = HDF5File(mesh_vel.mpi_comm(),os.path.join(data_dir, v3filename+".hdf5"), 'w')
f3.write(velocity, v3filename)
