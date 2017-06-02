"""
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

    Save path of directory of data files and parameters required in different functions which are constant

    :param
    cmaxF               Courant-Friedrich-Lewy condition parameter
    c_tol               concentration threshold for refinement
    decimals_mesh       round coordinates for mesh refinement (0 - only cut decimals, -1 - full tens)
    c_VOL               threshold of maximum volume (area in 2D) of refined cells in mesh
    timestep_times      factor to multiply time step computed with CFL condition (1 for CFL number=0.5)
    times_refine        factor distance of coordinates for refinement
    lc_param            parameter to set mesh size in 3D refinement
    startpoint          number of point to start list of points for refinement

"""

data_dir = '/PATH/TO/YOUR/CURRENT/WORKING/DIRECTORY'
# The following parameters have to be set according to the problem at hand
cmaxF = 0.5
c_tol = 1.e-2
decimals_mesh = 0
c_VOL = 0.075
timestep_times = 1
times_refine = 1.0
lc_param = 0.15
startpoint = 9
