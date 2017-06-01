"""
    This function refines a domain wide coarse mesh according to actual concentration patterns
    using a concentration threshold, cell area tolerance and refinement time step

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
from cfg_datadir import data_dir, decimals_mesh, times_refine, lc_param, startpoint, timestep_times, cmaxF, c_tol


def M_define(c0, meshR, velocity, velmax_old, dt_ADV):

    """
    Define coordinates of dofs where concentration exceeds threshold and update .geo file with this info

    :param c0: initial concentration
    :param meshR: refined mesh
    :return: velmax maximum velocity in refined area; test.geo file for generation of refined mesh using gmsh
    """

    # Define Function Spaces and their characteristics
    Q = FunctionSpace(meshR, "CG", 1)

    # Compute maximum movement and R_cells
    max_move = velmax_old * dt_ADV
    R_cells_write = max_move * times_refine

    # Get all dof coordinates
    gdim = meshR.geometry().dim()
    dofs_coors = Q.tabulate_dof_coordinates().reshape((-1, gdim))
    # Index array of dofs where concentration exceeds threshold c_tol and list of these dof's coordinates
    index_c0 = np.where(c0.vector()[:] > c_tol)
    dofs_c0 = dofs_coors[index_c0]

    # Round dof coordinates
    M_coordsN = np.ndarray.round(dofs_c0, decimals=decimals_mesh)
    # Delete non-unique coordinates - line 54 - 57 based on http://stackoverflow.com/a/31097280
    sorted_idx = np.lexsort(M_coordsN.T)
    sorted_Mcoords = M_coordsN[sorted_idx, :]
    row_mask = np.append([True], np.any(np.diff(sorted_Mcoords, axis=0),1))
    M_coords = sorted_Mcoords[row_mask]

    # Determine max velocity in refinement area and shift coordinates points in vel-direction
    velmax_vals = np.zeros(len(M_coords))
    M_new = np.zeros((len(M_coords), 3))
    for m in range(len(M_coords[:,0])):
        M_new[m,0] = M_coords[m,0] + velocity(M_coords[m,0], M_coords[m,1], M_coords[m,2])[0] * dt_ADV * 0.75
        M_new[m,1] = M_coords[m,1] + velocity(M_coords[m,0], M_coords[m,1], M_coords[m,2])[1] * dt_ADV * 0.75
        M_new[m,2] = M_coords[m,2] + velocity(M_coords[m,0], M_coords[m,1], M_coords[m,2])[2] * dt_ADV * 0.75
        velmax1 = velocity(M_coords[m,0], M_coords[m,1], M_coords[m,2])
        velmax2 = np.sqrt(np.sum(np.square(velmax1)))
        velmax_vals[m] = velmax2

    # Compute absolute velocity max in this area of refinement
    velmax = np.max(velmax_vals)

    # -----USER HAS TO SET PARAMETER--------If mesh refinement shall be performed without emphasis in streamline direction uncomment this line
    #M_new = M_coords.copy()

    # Define max distance for refinement
    R_cells_Max_write = R_cells_write * 10.0
    len_p = len(M_new[:,0])
    m1 = startpoint

    # Update coarse .geo file with info about refinement coordinates and R_cells
    with open("test.geo", "a") as myfile:
        for m in range(m1, m1 + len_p):
            myfile.write("Point(%d) = {%.2f, %.2f, %.2f, lc};\n" % (m, M_new[m - m1,0], M_new[m - m1,1], M_new[m - m1,2]))

        myfile.write("Field[1] = Attractor; \n")
        myfile.write("Field[1].NodesList = {%d" % (m1))
        for m in range(m1 + 1, m1 + len_p):
            myfile.write(",%d" % (m))
        myfile.write("};\n")
        myfile.write("Field[2] = Threshold;\nField[2].IField = 1;\nField[2].LcMin = %.2f;\n" % (lc_param))
        myfile.write("Field[2].LcMax = lc;\nField[2].DistMin = %.2f;\nField[2].DistMax = %.2f;\nField[3] = Min;\n" % (R_cells_write, R_cells_Max_write))
        myfile.write("Field[3].FieldsList ={2};\nBackground Field = 3;\n")

    return[velmax]


# +++++++++++++++++++++++++++++++ Main program ++++++++++++++++++++++++++++++++
# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters['allow_extrapolation'] = True
velmax_old = params['velmax']
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

# Read concentration and velocity from file
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

# Call M_define to save new mesh .geo file
[velmax] = M_define(c0, mesh_a, velocity, velmax_old, dt_ADV)


# Update parameters and save to file
params['velmax'] = velmax
with open(os.path.join(data_dir, "data.pkl"), "w") as f:
    pickle.dump(params, f)

