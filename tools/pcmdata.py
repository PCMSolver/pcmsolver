#
#  PCMSolver, an API for the Polarizable Continuum Model
#  Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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
CODATAdict = {
    2010: CODATA(0.52917721092, 2.418884326502e-02),
    2006: CODATA(0.52917720859, 2.418884326505e-02),
    2002: CODATA(0.52917721080, 2.418884326505e-02),
    1998: CODATA(0.52917720830, 2.418884326500e-02),
}

# Dictionary of known solvents, this is a mirror of the contents of Solvent.cpp
allowedSolvents = {
    'Water': ('WATER', 'H2O'),
    'Propylene Carbonate': ('PROPYLENE CARBONATE', 'C4H6O3'),
    'Dimethylsulfoxide': ('DIMETHYLSULFOXIDE', 'DMSO'),
    'Nitromethane': ('NITROMETHANE', 'CH3NO2'),
    'Acetonitrile': ('ACETONITRILE', 'CH3CN'),
    'Methanol': ('METHANOL', 'CH3OH'),
    'Ethanol': ('ETHANOL', 'CH3CH2OH'),
    'Acetone': ('ACETONE', 'C2H6CO'),
    '1,2-Dichloroethane': ('1,2-DICHLOROETHANE', 'C2H4CL2'),
    'Methylenechloride': ('METHYLENECHLORIDE', 'CH2CL2'),
    'Tetrahydrofurane': ('TETRAHYDROFURANE', 'THF'),
    'Aniline': ('ANILINE', 'C6H5NH2'),
    'Chlorobenzene': ('CHLOROBENZENE', 'C6H5CL'),
    'Chloroform': ('CHLOROFORM', 'CHCL3'),
    'Toluene': ('TOLUENE', 'C6H5CH3'),
    '1,4-Dioxane': ('1,4-DIOXANE', 'C4H8O2'),
    'Benzene': ('BENZENE', 'C6H6'),
    'Carbon Tetrachloride': ('CARBON TETRACHLORIDE', 'CCL4'),
    'Cyclohexane': ('CYCLOHEXANE', 'C6H12'),
    'N-heptane': ('N-HEPTANE', 'C7H16'),
    'Explicit': ('EXPLICIT', 'DUMMY')
}
