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


class CODATA:
    """A class holding conversion factors from atomic units

       Attributes:
           ToAngstrom conversion factor from AU of length to Angstrom
           ToFemtoseconds conversion factor from AU of time to Femtoseconds
    """

    def __init__(self, toang, tofs):
        self.ToAngstrom = toang
        self.ToFemtosecond = tofs

# Dictionary with the useful conversion factors for various CODATA sets of constants
CODATA2010 = CODATA(0.52917721092, 2.418884326502e-02)
CODATA2006 = CODATA(0.52917720859, 2.418884326505e-02)
CODATA2002 = CODATA(0.5291772108, 2.418884326505e-02)
CODATA1998 = CODATA(0.5291772083, 2.418884326500e-02)
CODATAdict = dict([(2010, CODATA2010), (2006, CODATA2006), (2002, CODATA2002), (1998, CODATA1998)])
