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

# Written by Jonas Juselius <jonas.juselius@chem.uit.no>
# University of Tromso, 2008
# Adapted to PCMSolver by Luca Frediani <luca.frediani@uit.no>
# University of Tromso, 2011
# Various modifications by Roberto Di Remigio <roberto.d.remigio@uit.no>
# University of Pisa, 2011-2012
# University of Tromso 2013-2017
# Virginia Tech 2017-
"""
Module collecting functions to parse the input to PCMSolver.
This is based on the Getkw library by J. Juselius.

Conventions:
routine names my_perfect_routine
keyword names MYPERFECTKEYWORD
"""

import sys
import tempfile
import os
from copy import deepcopy
import re

from .getkw import Section, GetkwParser
from .pcmdata import CODATAdict, allowedSolvents

isAngstrom = False
CODATAyear = 2010


def parse_pcm_input(inputFile, write_out=False):
    """
    Parses human-readable input to PCMSolver.

    The human-readable input file is first read and converted to uppercase.
    Parsing occurs after case conversion, so that input reading is case-insensitive.
    Optionally, save the result of parsing to file.

    Parameters
    ----------
    inputFile: str
        Full path to human-readable PCMSolver input.
    write_out: bool, optional
        Whether to save the result of parsing to file.
        If True, the name for the new file will be the name of the input file,
        prefixed by '@'.

    Returns
    -------
    parsed: str
        The parsed, machine-readable input.

    Example
    -------

    parsed = pcmsolver.parse_pcm_input(inp, write_out=True)

    ...
    """
    # Set up valid keywords.
    valid_keywords = setup_keywords()

    # Convert to uppercase and get the path to the temporary file
    uppercased = convert_to_upper_case(inputFile)

    # Set up a GetKw object and let it parse our input:
    # here is where the magic happens.
    inkw = GetkwParser().parseFile(uppercased)
    # Remove temporary file
    os.remove(uppercased)
    inkw.sanitize(valid_keywords)
    inkw.run_callbacks(valid_keywords)

    if write_out:
        # The parsed, machine-readable file is now saved.
        parsedFile = os.path.join(os.path.dirname(inputFile), '@' + os.path.basename(inputFile))
        with open(parsedFile, 'w') as tmp:
            tmp.write(str(inkw.top))

    return str(inkw.top)


def iequal(a, b):
    """
    Checks that strings a and b are equal, regardless of case.
    """
    try:
        return a.upper() == b.upper()
    except AttributeError:
        return a == b


def convert_to_upper_case(filename):
    """
    Reads contents of filename and converts them to uppercase.
    The case-converted contents are written back to a temporary file.
    """
    contents = ''
    final = ''
    with open(filename, 'r') as inputFile:
        contents = inputFile.readlines()
    # In case the "restart" field is present the case conversion
    # has to be done a bit more carefully. The filename argument to
    # "restart" should not be touched, whereas the rest of the file
    # has to be put in uppercase.
    npz_file = re.compile('NPZFILE', re.IGNORECASE)
    for line in contents:
        if npz_file.search(line):
            rst_split = re.split(r"=", line)
            rst_split = rst_split[0].upper() + ' = ' + rst_split[1]
            line = ''.join(rst_split)
        else:
            line = line.upper()
        final += line
    # Convert also occurences of TRUE/FALSE to True/False,
    # Python will otherwise interpret them as strings!
    final = re.sub('FALSE', 'False', final)
    final = re.sub('TRUE', 'True', final)

    # Write to temporary file in current working directory
    # The temporary file has the name of the input file, joined with a unique id
    temp, path = tempfile.mkstemp(prefix=filename + '.', dir=os.getcwd())
    with open(path, 'w') as outputFile:
        outputFile.write(final)
    os.close(temp)
    return path


def setup_keywords():
    """
    Sets up sections, keywords and respective callback functions.
    """
    # Top-level section
    top = Section('toplevel', callback=verify_top)
    top.set_status(True)
    # Define units of measure
    # Valid values: AU (atomic units) or ANGSTROM
    # Default: AU
    top.add_kw('UNITS', 'STR', 'AU')
    # Define set of CODATA constants
    # Valid values: 2010, 2006, 2002, 1998
    # Default: 2010
    top.add_kw('CODATA', 'INT', 2010)

    # Cavity section
    cavity = Section('CAVITY', callback=verify_cavity)
    # Type of the cavity
    # Valid values: GEPOL and RESTART
    cavity.add_kw('TYPE', 'STR')
    # Name of the file containing data for restarting a cavity
    # Valid for: Restart cavity
    # Default: empty string
    cavity.add_kw('NPZFILE', 'STR', '')
    # Average area (in au^2)
    # Valid for: GePol
    # Valid values: double strictly greater than 0.01 au^2
    # Default: 0.3 au^2
    cavity.add_kw('AREA', 'DBL', 0.3)
    # Scaling of the atomic radii
    # Valid for: GePol
    # Valid values: boolean
    # Default: True
    cavity.add_kw('SCALING', 'BOOL', True)
    # Built-in radii set
    # Valid for: GePol
    # Valid values: BONDI, UFF or ALLINGER
    # Default: BONDI
    cavity.add_kw('RADIISET', 'STR', 'BONDI')
    # Minimal radius of added spheres (in au)
    # Valid for: GePol
    # Valid values: double greater than 0.4 au
    # Default: 100.0 au (no added spheres)
    cavity.add_kw('MINRADIUS', 'DBL', 100.0)
    # Spheres geometry creation mode
    # Valid for: GePol
    # Valid values: EXPLICIT, ATOMS or IMPLICIT
    # Default: IMPLICIT
    cavity.add_kw('MODE', 'STR', 'IMPLICIT')
    # List of atoms with custom radius
    # Valid for: GePol in MODE=ATOMS
    # Valid values: array of integers
    cavity.add_kw('ATOMS', 'INT_ARRAY')
    # List of custom radii
    # Valid for: GePol in MODE=ATOMS
    # Valid values: array of doubles
    cavity.add_kw('RADII', 'DBL_ARRAY')
    # List of spheres
    # Valid for: GePol in MODE=EXPLICIT
    # Valid values: array of doubles in format [x, y, z, R]
    cavity.add_kw('SPHERES', 'DBL_ARRAY', callback=verify_spheres)
    top.add_sect(cavity)

    # Medium section
    medium = Section('MEDIUM', callback=verify_medium)
    # Type of solver
    # Valid values: IEFPCM, CPCM
    medium.add_kw('SOLVERTYPE', 'STR', 'IEFPCM')
    # Whether nonequilibrium response is to be used
    # Valid for: IEFPCM, CPCM
    # Valid values: boolean
    # Default: False
    medium.add_kw('NONEQUILIBRIUM', 'BOOL', False)
    # Solvent
    # Valid for: IEFPCM, CPCM
    # Valid values: string
    # Default: Explicit
    medium.add_kw('SOLVENT', 'STR', 'EXPLICIT')
    # Symmetrization of the PCM matrix
    # Valid for: IEFPCM
    # Valid values: boolean
    # Default: True
    medium.add_kw('MATRIXSYMM', 'BOOL', True)
    # CPCM correction factor
    # Valid for: CPCM
    # Valid values: positive double greater than 0.0
    medium.add_kw('CORRECTION', 'DBL', 0.0)
    # Type of integrator for the diagonal of the boundary integral operators
    # Valid for: IEFPCM, CPCM
    # Valid values: COLLOCATION
    # Default: COLLOCATION
    # Notes: In future releases we will add PURISIMA and NUMERICAL as additional options
    medium.add_kw('DIAGONALINTEGRATOR', 'STR', 'COLLOCATION')
    # Scaling factor for diagonal of collocation matrix
    # Valid for: IEFPCM, CPCM
    # Valid values: positive double greater than 0.0
    # Default: 1.07
    # Notes: values commonly used in the literature are 1.07 and 1.0694
    medium.add_kw('DIAGONALSCALING', 'DBL', 1.07)
    # Radius of the solvent probe (in au)
    # Valid for: IEFPCM, CPCM
    # Valid values: double in [0.1, 100.0] au
    # Default: 1.0
    medium.add_kw('PROBERADIUS', 'DBL', 1.0)
    top.add_sect(medium)

    # Green's function section
    green = Section('GREEN', callback=verify_green)
    # Green's function type
    # Valid values: VACUUM, UNIFORMDIELECTRIC, SPHERICALDIFFUSE, SPHERICALSHARP
    # Default: VACUUM
    green.add_kw('TYPE', 'STR', 'VACUUM')
    # Green's function derivative calculation strategy
    # Valid values: NUMERICAL, DERIVATIVE, GRADIENT, HESSIAN
    # Default: DERIVATIVE
    # Notes: all other values for debug purposes only
    green.add_kw('DER', 'STR', 'DERIVATIVE')
    # Static dielectric permittivity
    # Valid for: UNIFORMDIELECTRIC
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPS', 'DBL', 1.0)
    # Dynamic dielectric permittivity
    # Valid for: UNIFORMDIELECTRIC
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPSDYN', 'DBL', 1.0)
    # Dielectric profile type
    # Valid for: SPHERICALDIFFUSE
    # Valid values: TANH, ERF, LOG
    # Default: LOG
    green.add_kw('PROFILE', 'STR', 'LOG')
    # Static dielectric permittivity inside the interface
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPS1', 'DBL', 1.0)
    # Dynamic dielectric permittivity inside the interface
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPSDYN1', 'DBL', 1.0)
    # Static dielectric permittivity outside the interface
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPS2', 'DBL', 1.0)
    # Dynamic dielectric permittivity outside the interface
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double greater than 1.0
    # Default: 1.0
    green.add_kw('EPSDYN2', 'DBL', 1.0)
    # Center of the diffuse profile
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double
    # Default: 100.0
    # Notes: for SPHERICALDIFFUSE and SPHERICALSHARP corresponds to the sphere radius
    green.add_kw('CENTER', 'DBL', 100.0)
    # Width of the diffuse profile
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: positive double
    # Default: 5.0
    # Notes: this is used differently for different profiles
    green.add_kw('WIDTH', 'DBL', 5.0)
    # Center of the diffuse interface
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: array of doubles
    green.add_kw('INTERFACEORIGIN', 'DBL_ARRAY', [0.0, 0.0, 0.0])
    # Maximum angular momentum value
    # Valid for: SPHERICALDIFFUSE, SPHERICALSHARP
    # Valid values: integer
    green.add_kw('MAXL', 'INT', 50)
    medium.add_sect(green)

    green_part = deepcopy(green)
    green.add_sect(green_part)

    # Molecule section
    molecule = Section('MOLECULE')
    # List of geometry and classical point charges
    # Valid values: array of doubles in format [x, y, z, Q]
    # Notes: charges are always in atomic units
    molecule.add_kw('GEOMETRY', 'DBL_ARRAY', callback=verify_geometry)
    # Calculate the molecular electrostatic potential (MEP) at the cavity for the given molecule
    # Valid values: boolean
    # Default: True
    molecule.add_kw('MEP', 'BOOL', True)
    top.add_sect(molecule)

    # ChargeDistribution section
    # Set a classical charge distribution, inside or outside the cavity
    # No additional spheres will be generated.
    charge_distribution = Section('CHARGEDISTRIBUTION', callback=verify_charge_distribution)
    # Monopoles
    # Valid values: array of doubles in format [x, y, z, Q]
    # Notes: charges are always in atomic units
    charge_distribution.add_kw('MONOPOLES', 'DBL_ARRAY')
    # Dipoles
    # Valid values: array of doubles in format [x, y, z, mu_x, mu_y, mu_z]
    # Notes: dipole moment components are always in atomic units
    charge_distribution.add_kw('DIPOLES', 'DBL_ARRAY')
    top.add_sect(charge_distribution)

    return top


def verify_top(section):
    global isAngstrom, CODATAyear
    allowed_units = ('AU', 'ANGSTROM')
    key = section.get('UNITS')
    val = key.get()
    if (val not in allowed_units):
        print(('Units requested {} are not among the allowed units: {}'.format(val, allowed_units)))
        sys.exit(1)
    isAngstrom = True if (val == 'ANGSTROM') else False
    allowed_codata = (2010, 2006, 2002, 1998)
    CODATAyear = section.get('CODATA').get()
    if (CODATAyear not in allowed_codata):
        print(('CODATA set requested {} is not among the allowed sets: {}'.format(CODATAyear, allowed_codata)))
        sys.exit(1)


def verify_cavity(section):
    allowed = ('GEPOL', 'RESTART')
    type = section.get('TYPE')
    if (type.get() not in allowed):
        print(('Requested {} cavity is not among the allowed types: {}'.format(type, allowed)))
        sys.exit(1)

    # Convert units if input was given in Angstrom
    # The conversion functions check by themselves if the conversion is necessary or not!!
    if (section['AREA'].is_set()):
        convert_area_scalar(section['AREA'])
    if (section['MINRADIUS'].is_set()):
        convert_length_scalar(section['MINRADIUS'])

    if (type.get() == 'GEPOL'):
        area = section.get('AREA')
        a = area.get()
        if (a < 0.01):
            print(('Requested area value too small: {}. Minimal value is: 0.01 au^2'.format(a)))
            sys.exit(1)
        minRadius = section.get('MINRADIUS')
        mr = minRadius.get()
        if (mr < 0.4):
            print(('Requested minimal radius for added spheres too small: {}. Minimal value is: 0.4 au'.format(mr)))
            sys.exit(1)
    elif (type.get() == 'RESTART'):
        npzfile = section.get('NPZFILE')
        # Restart string is the filename, with extension, where the cavity specs were saved.
        # We only get the filename here, either an empty or a non-empty string,
        # further management of input is then done C++-side.
        if npzfile.get() == '':
            print('You need to specify a .npz filename for a restart...')
            sys.exit(1)

    radiiSet = section.get('RADIISET')
    allowed_sets = ('BONDI', 'UFF', 'ALLINGER')
    if (radiiSet.get() not in allowed_sets):
        print(('Radii set requested {} is not among the allowed sets: {}'.format(radiiSet.get(), allowed_sets)))
        sys.exit(1)
    allowed_modes = ('EXPLICIT', 'ATOMS', 'IMPLICIT')
    mode = section.get('MODE')
    if (mode.get() not in allowed_modes):
        print(('Cavity creation mode requested {} is not among the allowed modes: {}'.format(
            mode.get(), allowed_modes)))
        sys.exit(1)

    atoms = section.get('ATOMS')
    at = atoms.get()
    radii = section.get('RADII')
    convert_length_array(radii)
    r = radii.get()

    if (mode.get() == 'ATOMS'):
        if (len(r) != len(at) or len(at) == 0):
            print('Incoherent input for Atoms keyword. Check that Atoms and Radii are consistent.')
            sys.exit(1)
        else:
            for i, v in enumerate(at):
                if (at.count(v) > 1):
                    print('Incoherent input for Atoms keyword. Too many spheres on the same atom(s).')
                    sys.exit(1)


def verify_medium(section):
    solvent = section.get('SOLVENT')
    explicitSolvent = solvent.get() in allowedSolvents['Explicit']
    if (explicitSolvent):
        PRF = section.is_set('PROBERADIUS')
        GIF = section.is_set('GREEN<INSIDE>')
        GOF = section.is_set('GREEN<OUTSIDE>')
        if (not PRF):
            print('Error: Explicit solvent chosen but ProbeRadius not specified')
        if (not GIF):
            print('Error: Explicit solvent chosen but Green<inside> not specified')
        if (not GOF):
            print('Error: Explicit solvent chosen but Green<outside> not specified')
        if (not GIF or not GOF or not PRF):
            sys.exit(1)
    solventFound = False
    for i, v in allowedSolvents.items():
        if (solvent.get() in v):
            # Set name to the first value in the value pair
            # C++ will look for this name!
            solvent.set(v[0])
            solventFound = True
            break
    if (not solventFound):
        print(('Unknown solvent {}'.format(solvent.get())))
        print(('Choose a solvent from the following list:\n{}\n or specify the solvent data explicitly.'.format(
            list(allowedSolvents.keys()))))
        sys.exit(1)

    correction = section.get('CORRECTION')
    if (correction.get() < 0.0):
        print('Correction for CPCM solver must be greater than 0.0')
        sys.exit(1)

    integrator = section.get('DIAGONALINTEGRATOR')
    if (integrator.get() != 'COLLOCATION'):
        print('Only the collocation integrator is available')
        sys.exit(1)

    scaling = section.get('DIAGONALSCALING')
    if scaling.get() == 0.0:
        print('Scaling of diagonal for collocation matrices cannot be zero')
        sys.exit(1)

    convert_length_scalar(section.get('PROBERADIUS'))
    radius = section.get('PROBERADIUS')
    if (radius.get() < 0.1 or radius.get() > 100):
        print('Probe radius has to be within [0.1,100] Atomic Units')
        sys.exit(1)

    allowed_types = ('IEFPCM', 'CPCM')
    key = section.get('SOLVERTYPE')
    val = key.get()
    if (val not in allowed_types):
        print(('Allowed types are: {}'.format(allowed_types)))
        sys.exit(1)


def verify_green(section):
    allowed = ('VACUUM', 'UNIFORMDIELECTRIC', 'SPHERICALDIFFUSE', 'SPHERICALSHARP')
    allowed_der = ('NUMERICAL', 'DERIVATIVE', 'GRADIENT', 'HESSIAN')
    allowed_profiles = ('TANH', 'ERF', 'LOG')

    eps = section.get('EPS')
    epsdyn = section.get('EPSDYN')

    convert_length_scalar(section.get('CENTER'))
    convert_length_scalar(section.get('WIDTH'))
    convert_length_array(section.get('INTERFACEORIGIN'))

    type = section.get('TYPE')
    if (type.get() not in allowed):
        print(('Allowed Green\'s functions are: {}'.format(allowed)))
        sys.exit(1)

    der = section.get('DER')
    if (der.get() not in allowed_der):
        print(('Allowed derivatives strategies are: {}'.format(allowed)))
        sys.exit(1)

    if (type.get() == 'UNIFORMDIELECTRIC'):
        if not eps.is_set():
            print('Eps not defined for UniformDielectric')
            sys.exit(1)
        if not epsdyn.is_set():
            print('EpsDyn not defined for UniformDielectric')
            sys.exit(1)

    profile = section.get('PROFILE')
    if (profile.get() not in allowed_profiles):
        print(('Allowed profiles are: {}'.format(allowed_profiles)))
        sys.exit(1)


def check_array(name, array, offset):
    dim = len(array)
    if (dim % offset != 0):
        print(('Empty or incoherent {0} array'.format(name)))
        sys.exit(1)
    # Convert only geometry to AU. Leave the rest untouched
    if (isAngstrom):
        j = 0
        for i in range(dim // offset):
            array[j] /= CODATAdict[CODATAyear].ToAngstrom
            array[j + 1] /= CODATAdict[CODATAyear].ToAngstrom
            array[j + 2] /= CODATAdict[CODATAyear].ToAngstrom
            j += offset


def verify_geometry(keyword):
    data = keyword.get()
    check_array('GEOMETRY', data, 4)


def verify_charge_distribution(section):
    mono = section.get('MONOPOLES').get()
    check_array('MONOPOLES', mono, 4)
    dipole = section.get('DIPOLES').get()
    check_array('DIPOLES', dipole, 6)


def verify_spheres(keyword):
    length = len(keyword.get())
    if (length % 4 != 0):
        print('Empty or incoherent Spheres list.')
        sys.exit(1)
    convert_length_array(keyword)


def convert_length_array(keyword):
    length = len(keyword.get())
    if (isAngstrom):
        for i in range(length):
            keyword[i] /= CODATAdict[CODATAyear].ToAngstrom


def convert_length_scalar(keyword):
    if (isAngstrom):
        keyword[0] /= CODATAdict[CODATAyear].ToAngstrom


def convert_area_scalar(keyword):
    if (isAngstrom):
        keyword[0] /= (CODATAdict[CODATAyear].ToAngstrom * CODATAdict[CODATAyear].ToAngstrom)
