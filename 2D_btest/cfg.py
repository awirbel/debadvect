"""
    configuration script

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

import pickle
import os
from cfg_datadir import data_dir
import numpy as np

"""

    Save parameters that have to be passed between functions in a pickle as dictionary

    :param
    fname           name prefix of gmsh mesh files
    D_diff          Diffusion constant
    mass_initial    mass of initial concentration at t=0
    mass_previous   mass of previous time step
    velmax          max velocity in refined area
    dt_total_old    total time step in previous computation step
    dt_ADV          advection time step
    dt              actual time step used in advection
    t1              number of actual advection time step
"""

run_params = dict(
    fname = "box",
    D_diff = 1.e-9,
    mass_initial = 0.0,
    mass_previous = 0.0,
    velmax = 0.0,
    dt_total_old = 0.0,
    dt_ADV = 0.1 * np.pi,
    dt = 0.0,
    t1 = 0
    )


with open(os.path.join(data_dir, "data.pkl"), "wb" ) as f:
    pickle.dump(run_params, f)
