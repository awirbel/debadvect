"""
    This function saves all the mesh files as .hdf5 to be handled by the other modules and enable parallel computation

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
import pickle
import os
from cfg_datadir import data_dir


# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']

# Load meshes from file
mesh_a = Mesh(os.path.join(data_dir, ("%s.xml" % (fname))))
mesh_c = Mesh(os.path.join(data_dir, ("%s_fine.xml" % (fname))))
mesh_b = Mesh(os.path.join(data_dir, ("%s_coarse.xml" % (fname))))

# Read subdomains and boundaries from mesh file
subdomains_a = MeshFunction("size_t", mesh_a, os.path.join(data_dir, ("%s_physical_region.xml" % (fname))))
boundaries_a = MeshFunction("size_t", mesh_a, os.path.join(data_dir, ("%s_facet_region.xml" % (fname))))
boundaries_c = MeshFunction("size_t", mesh_c, os.path.join(data_dir, ("%s_fine_facet_region.xml" % (fname))))

# Save subdomains and boundaries to files that can be read also in parallel
b1filename = "%s_fine_boundaries" % (fname)
fu = HDF5File(mesh_c.mpi_comm(),os.path.join(data_dir, b1filename+".hdf5"), 'w')
fu.write(boundaries_c, b1filename)

bfilename = "%s_boundaries" % (fname)
fu0 = HDF5File(mesh_a.mpi_comm(),os.path.join(data_dir, bfilename+".hdf5"), 'w')
fu0.write(boundaries_a, bfilename)

mfilename = "%s" % (fname)
fu1 = HDF5File(mesh_a.mpi_comm(),os.path.join(data_dir, mfilename+".hdf5"), 'w')
fu1.write(mesh_a, mfilename)
m2filename = "%s_coarse" % (fname)
fu2 = HDF5File(mesh_b.mpi_comm(),os.path.join(data_dir, m2filename+".hdf5"), 'w')
fu2.write(mesh_b, m2filename)
m3filename = "%s_fine" % (fname)
fu3 = HDF5File(mesh_c.mpi_comm(),os.path.join(data_dir, m3filename+".hdf5"), 'w')
fu3.write(mesh_c, m3filename)

# Save coarse mesh to VTK for visualisation
meshfile = File(os.path.join(data_dir, ("coarsemesh.pvd")))
meshfile << mesh_b

