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
from cfg_datadir import data_dir, c_tol, cmaxF, c_VOL, times_refine, decimals_mesh, timestep_times



def M_define(c0, meshR, meshC, velocity, velmax_old, dt_ADV):

    """
    Mark all cells where concentration exceeds threshold c_tol and all cells within distance R_cells

    :param c0: initial concentration
    :param meshC: coarse mesh
    :param meshR: refined mesh
    :return: cell_markers2: CellFunction that indicates where the coarse mesh has to be refined,
            M_cells: Function representing the area to be refined, velmax: maximum velocity in refined area
    """

    # Define Function Spaces and their characteristics
    Q = FunctionSpace(meshR, "CG", 1)
    Q_coarse = FunctionSpace(meshC, "CG", 1)
    dofmap_coarse = Q_coarse.dofmap()
    cell_dofs = []

    # Compute maximum movement and R_cells
    max_move = velmax_old * dt_ADV
    R_cells = max_move * times_refine

    # Get all dof coordinates
    gdim = meshR.geometry().dim()
    dofs_coors = Q.tabulate_dof_coordinates().reshape((-1, gdim))
    # Index array of dofs where concentration exceeds threshold c_tol and list of these dof's coordinates
    index_c0 = np.where(c0.vector()[:] > c_tol)
    dofs_c0 = dofs_coors[index_c0]

    # Round dof coordinates
    M_coordsN = np.ndarray.round(dofs_c0, decimals=decimals_mesh)
    # Delete non-unique coordinates - line 62 - 65 based on http://stackoverflow.com/a/31097280
    sorted_idx = np.lexsort(M_coordsN.T)
    sorted_Mcoords = M_coordsN[sorted_idx, :]
    row_mask = np.append([True], np.any(np.diff(sorted_Mcoords, axis=0),1))
    M_coords = sorted_Mcoords[row_mask]

    # Determine max velocity in refinement area and shift coordinates points in vel-direction
    velmax_vals = np.zeros(len(M_coords))
    M_new = np.zeros((len(M_coords), 2))
    for m in range(len(M_coords[:,0])):
        M_new[m,0] = M_coords[m,0] + velocity(M_coords[m,0], M_coords[m,1])[0] * dt_ADV * 0.75
        M_new[m,1] = M_coords[m,1] + velocity(M_coords[m,0], M_coords[m,1])[1] * dt_ADV * 0.75
        velmax1 = velocity(M_coords[m,0], M_coords[m,1])
        velmax2 = np.sqrt(np.sum(np.square(velmax1)))
        velmax_vals[m] = velmax2

    # Compute absolute velocity max in this area of refinement
    velmax = np.max(velmax_vals)

    # -----USER HAS TO SET PARAMETER--------If mesh refinement shall be performed without emphasis in streamline direction uncomment this line
    # M_new = M_coords.copy()

    # Mark cells that are within R_cells of any of these coordinates and write their dofs to a list
    cell_markers2 = CellFunction("bool", meshC)
    for cell in cells(meshC):
        for i in range(len(M_coords)):
            if cell.distance(Point(M_new[i])) < R_cells:
                cell_dofs.extend(dofmap_coarse.cell_dofs(cell.index()))
                cell_markers2[cell] = True

    # Set M_cells to 100.0 for respective dofs
    M_cells = Function(Q_coarse)
    M_cells.vector()[cell_dofs] = 100.0

    return[M_cells, velmax, cell_markers2]


def refine_ui(M_cells, meshC, velocity, cell_markers2):

    """
    Refine mesh according to marked cells and function M_cells until cell-area based threshold is reached

    :param M_cells, cell_markers2: indicator for where to refine mesh
    :param meshC: coarse mesh
    :param velocity: velocity
    :return: mesh_new - refined mesh; conc_new, vel_new, bound_new - concentration, velocity and boundaries on
            refined mesh

    """

    # First refinement using cell_markers2 function
    adapt(meshC, cell_markers2)
    meshC = meshC.child()
    flag_ref = 1

    # Further refinement using M_cells until all marked cell area is smaller than c_VOL
    while flag_ref == 1:
        cell_markers = CellFunction("bool", meshC)
        cell_markers.set_all(False)
        flag_ref = 0
        for cell in cells(meshC):
            p = cell.midpoint()
            cv = cell.volume()
            if (abs(M_cells(p))) > 90.0 and cv > c_VOL:
                cell_markers[cell] = True
                flag_ref = 1

        # Refine mesh according to indicators and tolerances
        adapt(meshC, cell_markers)
        meshC = meshC.child()

    # Set refined mesh to mesh_new
    mesh_new = meshC

    # Adapt functions to final mesh
    adapt(c0, mesh_new)
    conc_init = c0.child()
    adapt(velocity, mesh_new)
    vel_new = velocity.child()

    # Mark exterior boundary
    bound_new = FacetFunction("size_t", mesh_new)
    for f in facets(mesh_new):
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
velmax_old = params['velmax']
dt_ADV = params['dt_ADV']
t1 = params['t1']

# Load mesh from file
mesh_vel = Mesh()
mname = "%s" % (fname)
hdf = HDF5File(mesh_vel.mpi_comm(), os.path.join(data_dir, mname+".hdf5"), "r")
hdf.read(mesh_vel, mname, False)
hdf.close()
mesh_coarse = Mesh()
m2name = "%s_coarse" % (fname)
hd2f = HDF5File(mesh_coarse.mpi_comm(), os.path.join(data_dir, m2name+".hdf5"), "r")
hd2f.read(mesh_coarse, m2name, False)
hd2f.close()
mesh_a = Mesh()
maname = "refined%s" % (fname)
hdf2 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, maname+".hdf5"), "r")
hdf2.read(mesh_a, maname, False)
hdf2.close()
mesh = mesh_vel

# Define Function Spaces
Q = FunctionSpace(mesh_a, "CG", 1)
V2 = VectorFunctionSpace(mesh_vel, "CG", 2)
V_coarse = VectorFunctionSpace(mesh_coarse, "CG", 2)

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
velocity_coarse = Function(V_coarse)
velin2name = "velocity_coarse_%s" % (fname)
hdf5 = HDF5File(mesh_coarse.mpi_comm(), os.path.join(data_dir, velin2name+".hdf5"), "r")
hdf5.read(velocity_coarse, velin2name)
hdf5.close()

# Perform mesh refinement
[M_cells, velmax, cell_markers2] = M_define(c0, mesh_a, mesh_coarse, velocity, velmax_old, dt_ADV)
[mesh_new, conc_init, vel_new, bound_new] = refine_ui(M_cells, mesh_coarse, velocity, cell_markers2)

# Check mass conservation for refinement and interpolation on newly refined mesh
m_previous = assemble(c0 * dx)
m_new = assemble(conc_init * dx)
print "Masses old: ", m_previous, "and new: ",  m_new
print "due to interpolation the new mass is ", (100. / m_previous) * m_new, "percent"
if m_new < (0.95 * m_previous) or m_new > (1.05 * m_previous):
    print "NEW MESH REFINEMENT has led to mass loss/gain > 5%"

# Define time step for subsequent advection computations
dt1 = ((cmaxF * mesh_new.hmin()) / (velmax))
if t1 == 0:
    Q_1 =FunctionSpace(mesh_new, "CG", 1)
    mass_now = assemble(conc_init*dx)
    params['mass_start'] = mass_now

# Update parameters and save to File
params['velmax'] = velmax
params['dt'] = dt1
with open(os.path.join(data_dir, "data.pkl"), "w") as f:
    pickle.dump(params, f)

# Save updated mesh and functions fo file
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

