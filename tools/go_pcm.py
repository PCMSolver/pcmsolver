#
#  PCMSolver, an API for the Polarizable Continuum Model
#  Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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
#
# Python front-end to run_pcm standalone PCMSolver executable
#
# Written by Roberto Di Remigio <roberto.d.remigio@uit.no>
"""
Parses the PCMSolver input file and launches a standalone calculation with the
run_pcm executable.
"""

import docopt
import os
import platform
import subprocess
import sys

from pcmsolver import parse_pcm_input

options = """
Usage:
    ./pcmsolver.py --inp=<input_file> --exe=<path_to_exe>
    ./pcmsolver.py (-h | --help)

Options:
  --exe=<path_to_exe>  Path to run_pcm executable.
  --inp=<input_file>   PCMSolver input file.
  -h --help            Show this screen.
"""


def main():
    try:
        arguments = docopt.docopt(options, argv=None)
    except docopt.DocoptExit:
        sys.stderr.write('ERROR: bad input to {}\n'.format(sys.argv[0]))
        sys.stderr.write(options)
        sys.exit(-1)
    input_file = arguments['--inp']
    path_to_exe = arguments['--exe']
    # Do the parsing and return validated input
    parsed = parse_pcm_input(input_file)
    # The parsed, machine-readable file is now saved.
    parsedFile = os.path.join(os.path.dirname(input_file), '@' + os.path.basename(input_file))
    # Any leftover parsedFile will be overwritten
    with open(parsedFile, 'w') as tmp:
        tmp.write(parsed)
    # Execute standalone
    run_pcm(path_to_exe, parsedFile)


def run_pcm(path_to_exe, parsedFile):
    suffix = '.exe' if platform.system() == 'Windows' else ''
    executable = os.path.join(path_to_exe, 'run_pcm' + suffix)
    p = subprocess.Popen(
        [executable, parsedFile], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_coded, stderr_coded = p.communicate()
    stdout = stdout_coded.decode('UTF-8')
    stderr = stderr_coded.decode('UTF-8')

    print(stdout)
    # Print the contents of stderr, without exiting
    sys.stderr.write(stderr)


if __name__ == '__main__':
    main()
