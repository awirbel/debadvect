"""

    This function reads the refined mesh and adapts functions onto the refined mesh

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
import numpy as np
import pickle
import os
from cfg_datadir import data_dir, cmaxF


def refine_ui(velocity, u0):

    """
    Adapt velocity and initial concentration to new refined mesh

    :param c0: concentration
    :param velocity: velocity
    :return: mesh_new - refined mesh; conc_init, vel_new, bound_new - initial concentration, velocity, boundaries on
            refined mesh

    """

    # Load refined mesh from file
    mesh_new = Mesh("test.xml")

    # Adapt functions to refined mesh
    adapt(u0, mesh_new)
    conc_init = u0.child()
    adapt(velocity, mesh_new)
    vel_new = velocity.child()

    # Mark exterior boundary
    bound_new = FacetFunction("size_t", mesh_new)
    for f in facets(mesh_new):
        # All exterior facets
        if f.exterior():
                bound_new[f] = 1


    return[mesh_new, conc_init, vel_new, bound_new]


# +++++++++++++++++++++++++++++++ Main program ++++++++++++++++++++++++++++++++
# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters['allow_extrapolation'] = True
velmax = params['velmax']
dt_ADV = params['dt_ADV']
t1 = params['t1']

# Load mesh from file
mesh_vel = Mesh()
mname = "%s" % (fname)
hdf = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir, mname+".hdf5"), "r")
hdf.read(mesh_vel, mname, False)
hdf.close()
mesh_a = Mesh()
maname = "refined%s" % (fname)
hdf2 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, maname+".hdf5"), "r")
hdf2.read(mesh_a, maname, False)
hdf2.close()
mesh = mesh_vel

# Define Function Spaces
Q = FunctionSpace(mesh_a, "CG", 1)
V2 = VectorFunctionSpace(mesh_vel, "CG", 2)

# Read concentration and velocities from file
c0 = Function(Q)
uinname = "concentration_refine_%s" % (fname)
hdf3 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, uinname+".hdf5"), "r")
hdf3.read(c0, uinname)
hdf3.close()
velocity = Function(V2)
velinname = "velocity_%s" % (fname)
hdf4 = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir, velinname+".hdf5"), "r")
hdf4.read(velocity, velinname)
hdf4.close()

# call refine_ui to adapt concentration and velocity to the new mesh
[mesh_new, conc_init, vel_new, bound_new] = refine_ui(velocity, c0)


# check if mass conservation due to refinement and interpolation on newly refined mesh
m_previous = assemble(c0 * dx)
m_new = assemble(conc_init * dx)
print "Masses", m_previous, m_new
print "due to interpolation the new mass is ", (100. / m_previous) * m_new, "percent"
if m_new < (0.95 * m_previous) or m_new > (1.05 * m_previous):
    print "NEW MESH REFINEMENT has led to mass loss/gain > 5%"

# Define time step for subsequent advection computations
dt1 = ((cmaxF * mesh_new.hmin()) / (velmax))
if t1 == 0:
    Q_1 =FunctionSpace(mesh_new, "CG", 1)
    mass_now = assemble(conc_init*dx)
    params['mass_start'] = mass_now

# Update parameters and save to file
params['velmax'] = velmax
params['dt'] = dt1
with open(os.path.join(data_dir, "data.pkl"), "w") as f:
    pickle.dump(params, f)


# Save updated mesh and functions to file
bfilename = "%s_refined_boundaries" % (fname)
fb = HDF5File(mesh_new.mpi_comm(),os.path.join(data_dir, bfilename+".hdf5"), 'w')
fb.write(bound_new, bfilename)

velfilename = "velocityREADfine2_%s" % (fname)
fv = HDF5File(mesh_new.mpi_comm(),os.path.join(data_dir, velfilename+".hdf5"), 'w')
fv.write(vel_new, velfilename)

ufilename = "concentration_initial2_%s" % (fname)
fu = HDF5File(mesh_new.mpi_comm(),os.path.join(data_dir, ufilename+".hdf5"), 'w')
fu.write(conc_init, ufilename)

mfilename = "refined%s" % (fname)
fm = HDF5File(mesh_new.mpi_comm(),os.path.join(data_dir, mfilename+".hdf5"), 'w')
fm.write(mesh_new, mfilename)
