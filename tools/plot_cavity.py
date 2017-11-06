#
#  PCMSolver, an API for the Polarizable Continuum Model
#  Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
#
#  This file is part of PCMSolver.
#
#  PCMSolver is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PCMSolver is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
#
#  For information on the complete list of contributors to the
#  PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
#

# -*- python -*-
# -*- coding: utf-8 -*-
# vim:filetype=python:

# Written by Roberto Di Remigio <roberto.d.remigio@uit.no>
# University of Tromso, 2017
"""
Plot molecular cavity from a compressed NumPy format file.
Color map the finite elements according to a surface function, saved
to NumPy format file.
"""

import os
import sys
sys.path.append(os.path.dirname(__file__))

try:
    import docopt
except:
    sys.path.append('cmake/lib/docopt')
    import docopt

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

options = """
Usage:
    ./plot_cavity.py <cavity_npz> [--map-by <npy>]
    ./plot_cavity.py (-h | --help)

Options:
  <cavity_npz>   Compressed NumPy file with cavity specifications.
  --map-by <npy> NumPy format file with surface function to color-map finite elements.
  -h --help      Show this screen.
"""


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(
            0.0, midpoint, 128, endpoint=False), np.linspace(
                midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def plot(cavity_npz, surf_func_npy=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    cavity = np.load(cavity_npz)

    nElements = cavity['elements']
    centroids = cavity['centers']

    # Plot collocation points
    ax.scatter(centroids[0, :], centroids[1, :], centroids[2, :], c='black', alpha=0.5)

    # Generate color mapping
    colors = (.5, .1, .3, 0.3)
    if surf_func_npy:
        surf_func = np.load(surf_func_npy)
        shifted_cmap = shiftedColorMap(cm.coolwarm, midpoint=0.75, name='shifted')
        mappable = cm.ScalarMappable(cmap=shifted_cmap)
        mappable.set_array(surf_func.flatten())
        plt.colorbar(mappable)
        # Provide colors for Poly3DCollection
        colors = mappable.to_rgba(surf_func.flatten())
    # Generate list of vertices
    vertices = [
        zip(cavity['vertices_' + str(i)][0, :], cavity['vertices_' + str(i)][1, :], cavity['vertices_' + str(i)][2, :])
        for i in range(nElements)
    ]
    elements = Poly3DCollection(vertices, facecolors=colors)
    ax.add_collection3d(elements)
    ax.set_axis_off()
    plt.show()


def main():
    try:
        arguments = docopt.docopt(options, argv=None)
    except docopt.DocoptExit:
        sys.stderr.write('ERROR: bad input to %s\n' % sys.argv[0])
        sys.stderr.write(options)
        sys.exit(-1)
    cavity_npz = arguments['<cavity_npz>']
    surf_func_npy = arguments['--map-by']
    plot(cavity_npz, surf_func_npy)


if __name__ == '__main__':
    main()
