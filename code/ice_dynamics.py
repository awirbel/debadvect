"""
    icetools code Newton iteration version

    original version can be found: https://github.com/alexjarosch/icetools

    Copyright 2015-2019 Alexander H. Jarosch, modified by Anna Wirbel

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
parameters["form_compiler"]["quadrature_degree"] = 2

# Load mesh from file
mesh = Mesh()
mname = "%s" % (fname)
hdfin3 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, mname+".hdf5"), "r")
hdfin3.read(mesh, mname, False)

# Load boundaries, free surface tag(2); bed tag(1)
boundaries = FacetFunction("size_t", mesh)
binname = "%s_boundaries" % (fname)
hdfin2 = HDF5File(mesh.mpi_comm(), os.path.join(data_dir, binname+".hdf5"), "r")
hdfin2.read(boundaries, binname)

# Define some general constants
g = -9.81               # gravitational constant
rho = 917.0             # fluid density
Aglen = 2.4e-24         # Glen flow parameter for temperate ice (Cuffey & Paterson,2010 p. 73)
nglen = 3.0             # Glen's n
glen_fact = 0.5 * Aglen**(-1.0/nglen)

# Define the body force f, i.e. gravity, driving the flow
f_x0 = 0.0
# In 2D use lines 58-59 and comment lines 61-63:
f_x1 = g*rho
f = Constant((f_x0, f_x1))
# In 3D use lines 61-63 and comment lines 58-59:
#f_x1 = 0.0
#f_x2 = g*rho
#f = Constant((f_x0, f_x1, f_x2))

# Define Function Spaces
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define solution function
w = Function(W)
(u, p) = split(w)
(v, q) = TestFunctions(W)

# Apply a no-slip boundary condition for velocity
# In 2D use line 78 and comment line 80:
noslip = Constant((0.0,0.0))
# In 3D use line 80 and comment line 78:
#noslip = Constant((0.0, 0.0, 0.0))
bc0 = DirichletBC(W.sub(0), noslip, boundaries, 1)
# Collect boundary conditions
bcs = [bc0]

##### New Code ######

nu = 8e13 # linear viscosity for for first linear guess

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 5
solver.parameters["newton_solver"]["linear_solver"] = 'mumps'
solver.solve()

# Define the nonlinear stokes equation directly:
def visc(u):

    eps_dot = sqrt(0.5*inner(sym(grad(u)),sym(grad(u)))) # second invariant of strain
    nu_out = glen_fact * eps_dot**((1.0 - nglen)/nglen)
    # return nu_out
    return Min(nu_out, 2e15)    # would introduce a viscosity limit

nu = visc(u)

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 80
solver.parameters["newton_solver"]["error_on_nonconvergence"] = False
solver.parameters["newton_solver"]["relaxation_parameter"] = 0.6
solver.parameters["newton_solver"]["relative_tolerance"] = 1E-4
solver.parameters["newton_solver"]["linear_solver"] = 'mumps'
solver.solve()

(u, p) = w.split(deepcopy=True)

# Save solution to file
ufilename = "velocity_%s" % (fname)
fu = HDF5File(mesh.mpi_comm(),os.path.join(data_dir, ufilename+".hdf5"), 'w')
fu.write(u, ufilename)
uendfile_pvd = File(os.path.join(data_dir, ("velocity_%s.pvd" % (fname))))
uendfile_pvd << u
