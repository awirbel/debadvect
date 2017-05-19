"""
    This model solves the time dependent advection-diffusion equation using FEM
    SUPG is used in terms of numerical stability issues

    Copyright 2017 Anna Wirbel and Alexander H. Jarosch
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

    This module uses portions of on the undocumented fenics demo:
    demo_advection-diffusion.py Copyright (C) 2007 by Kristian Oelgaard as part of DOLFIN
    under the GNU Lesser Public License <http://www.gnu.org/licenses/>.

"""

from dolfin import *
import numpy as np
import pickle
import os
from cfg_datadir import data_dir, cmaxF, timestep_times

# Load parameters
with open(os.path.join(data_dir, "data.pkl"), "r") as f:
    params = pickle.load(f)

# Define parameters from file
fname = params['fname']
D_diff = params['D_diff']
velmax = params['velmax']
t1 = params['t1']
dt_total_old = params['dt_total_old']
dt = params['dt']
dt_ADV = params['dt_ADV']
mass_initial = params['mass_initial']
mass_previous = params['mass_previous']
cmax = cmaxF
parameters['allow_extrapolation'] = True
# set the form compiler quadrature degree lower in case of 3D modelling
# parameters["form_compiler"]["quadrature_degree"] = 2

# time stepping
dt = dt * timestep_times
t = dt
dt_form = Constant(dt)

# Load mesh from file
mesh_a = Mesh()
mainname = "refined%s" % (fname)
hdfin1 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, mainname+".hdf5"), "r")
hdfin1.read(mesh_a, mainname, False)
minmeshsize = mesh_a.hmin()

# Required in SUPG stabilisation
h = CellSize(mesh_a)
# Define a Function Space for the concentration (Scalar Function)
Q = FunctionSpace(mesh_a, "CG", 1)
# Define VectorFunctionSpace for velocities
V2 = VectorFunctionSpace(mesh_a, "CG", 2)

# Zero source function
f = Constant(0.0)
# Read velocity and starting concentration from file
velocity = Function(V2)
c0 = Function(Q)
velfilename = "velocityREADfine2_%s" % (fname)
hdfin2 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir,velfilename+".hdf5"), 'r')
hdfin2.read(velocity, velfilename)
hdfin2.close()
uinname = "concentration_initial2_%s" % (fname)
hdfin3 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, uinname+".hdf5"), "r")
hdfin3.read(c0, uinname)
hdfin3.close()

# ********************Set up FEM problem********************
# Test and trial functions
c, v = TrialFunction(Q), TestFunction(Q)

# Mid-point solution
c_mid = 0.5 * (c0 + c)

# Residual
r = c - c0 + dt_form * (dot(velocity, grad(c_mid)) - D_diff * div(grad(c_mid)) - f)

# Galerkin variational problem
F1 = v * (c - c0) * dx + dt_form * (v * dot(velocity, grad(c_mid)) * dx + D_diff * dot(grad(v), grad(c_mid)) * dx - f * v * dx)

# Add SUPG stabilisation terms
vnorm = sqrt(dot(velocity, velocity))
F = F1 + (h / (2.0 * vnorm)) * dot(velocity, grad(v)) * r * dx

# Create bilinear and linear forms
a = lhs(F)
L = rhs(F)

# Set Dirichlet boundary condition and read domain boundaries from file
boundary_DC = FacetFunction("size_t", mesh_a)
binname = "%s_refined_boundaries" % (fname)
hdfin4 = HDF5File(mesh_a.mpi_comm(), os.path.join(data_dir, binname+".hdf5"), "r")
hdfin4.read(boundary_DC, binname)
hdfin4.close()
# Define values and set up Dirichlet boundary condition(s)
preC = Constant(0.0)
bc = DirichletBC(Q, preC, boundary_DC, 1)

# Assemble matrix, only once before time step and apply boundary condition(s)
A = assemble(a)
bc.apply(A)
# Create linear solver and factorize matrix
parameters['lu_solver']["reuse_factorization"] = True

# Set intial condition
c = c0
# Compute initial mass of this time step
mass_initdt = assemble(c0*dx)

# Define parameters for time stepping loop
t_ref = 1
T_ref = int(dt_ADV / dt) + 2
dt_final = dt_ADV * (t1+1)
flag_continue = 0
flag_stop = 0

# Time-stepping loop
while flag_stop < 0.5:

    print "advection time step: ", t1, " in computation: ", t_ref, "of", T_ref, "dt is: ", dt
    # Assemble vector and apply boundary condition(s)
    b = assemble(L)
    bc.apply(b)

    if flag_continue == 1:
        # Change dt to finish this refinement time step
        dt_form.assign(dt_new)
        dt_total = dt_ADV * (t1+1)
        # Assemble matrix again as dt has changed and assemble vector
        A = assemble(a)
        b = assemble(L)
        # Apply boundary condition(s)
        bc.apply(A)
        bc.apply(b)
        solve(A, c.vector(), b, "mumps", "default")
        # Set stopping flag
        flag_stop = 1
    else:
        # Perform computations
        solve(A, c.vector(), b, "mumps", "default")
        dt_total = t_ref * dt + dt_total_old

    # Copy solution from previous interval
    c0 = c
    t_ref = t_ref + 1

    print "total time step: ", dt_total, "mass in % of this: ", (100. /mass_initdt) *assemble(c*dx), "mass in % of initial: ",(100. / mass_initial) * assemble(c*dx)
    # Check if refinement time step exceeded within next computation, if TRUE change dt
    if (dt_total + dt) > dt_final:
        dt_new = dt_final - dt_total
        dt_form.assign(dt_new)
        flag_continue = 1


# Write solution to file, to be used in next refinement time step
ufilename = "concentration_refine_%s" % (fname)
fu = HDF5File(mesh_a.mpi_comm(),os.path.join(data_dir, ufilename+".hdf5"), 'w')
fu.write(c, ufilename)

# Write solution to VTK for visualisation
initialbfile_pvd = File(os.path.join(data_dir, ("concentration_t%d.pvd" % (t1))))
initialbfile_pvd << c

# Compute mass conservation parameters
mass_total = assemble(c * dx)
mass_percent = (100. / mass_initdt) * mass_total
mass_percent_initial = (100. / mass_initial) * mass_total
mass_percent_dt = (100. / mass_previous) * mass_initdt

print "mass conservation parameters: "
print "total mass: ", mass_total, " in percent of this time step starting mass: ", mass_percent
print "and compared to initial mass: ", mass_percent_initial

# Update parameters and save to file
params['t1'] = t1 + 1
params['dt_total_old'] = dt_total
params['mass_previous'] = mass_total
with open(os.path.join(data_dir, 'data.pkl'), "w") as f:
    pickle.dump(params, f)
