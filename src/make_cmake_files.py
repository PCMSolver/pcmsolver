#!/usr/bin/python

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
# Create CMakeLists.txt template for leaf directories
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import os
import sys
import glob

sys.path.append('cmake/lib/docopt')
import docopt

options = """
Usage:
  ./make_cmake_files.py [options]
  ./make_cmake_files.py (-h | --help)

Options:
  --libname=<LIBNAME>                    Name of the library to be created [default: ''].
  --lang=<LANGUAGE>                      Source file language <CXX/C/F> [default: CXX].
  -h --help                              Show this screen.
"""


def glob_sources_cxx(dir_name):
    """Create a list of C++ headers and sources to be used in a CMakeLists.txt file."""
    headers = 'list(APPEND headers_list '
    headers += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sorted(glob.glob(dir_name + '/*.hpp')))
    headers += ')\n\n'

    sources = 'list(APPEND sources_list '
    sources += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sorted(glob.glob(dir_name + '/*.cpp')))
    sources += ')\n\n'

    message = '# List of headers\n' + headers \
            + '# List of sources\n' + sources
    return message


def glob_sources_c(dir_name):
    """Create a list of C headers and sources to be used in a CMakeLists.txt file."""
    headers = 'list(APPEND headers_list '
    headers += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sorted(glob.glob(dir_name + '/*.h')))
    headers += ')\n\n'

    sources = 'list(APPEND sources_list '
    sources += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sorted(glob.glob(dir_name + '/*.c')))
    sources += ')\n\n'

    message = '# List of headers\n' + headers \
            + '# List of sources\n' + sources
    return message


def glob_sources_fortran(dir_name):
    """Create a list of Fortran sources to be used in a CMakeLists.txt file."""
    types = ('*.f', '*.F', '*.f77', '*.F77', '*.f90', '*.F90')
    list_of_sources = []
    for ftype in types:
        list_of_sources.extend(sorted(glob.glob(dir_name + '/' + ftype)))
    sources = 'list(APPEND sources_list '
    sources += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in list_of_sources)
    sources += ')\n\n'

    message = '# List of sources\n' + sources
    return message


try:
    arguments = docopt.docopt(options, argv=None)
except docopt.DocoptExit:
    sys.stderr.write('ERROR: bad input to %s\n' % sys.argv[0])
    sys.stderr.write(options)
    sys.exit(-1)

# Grab command-line arguments
libname = arguments['--libname']
lang = arguments['--lang']

root_directory = os.getcwd()
dname = os.path.join(root_directory, libname)
fname = os.path.join(dname, 'CMakeLists.txt.try')
f = open(fname, 'w')
language = ''
if (lang == 'CXX'):
    f.write(glob_sources_cxx(dname))
    language = 'CXX'
elif (lang == 'C'):
    f.write(glob_sources_c(dname))
    language = 'C'
else:
    f.write(glob_sources_fortran(dname))
    language = 'Fortran'

if (lang == 'CXX'):
    f.write('add_library(' + libname + ' OBJECT ${sources_list} ${headers_list})\n')
    f.write('set_target_properties(' + libname + ' PROPERTIES POSITION_INDEPENDENT_CODE 1 )\n')
    f.write('add_dependencies(' + libname + ' generate-config-hpp)\n')
    f.write('target_compile_options(' + libname + ' PRIVATE "$<$<CONFIG:DEBUG>:${EXDIAG_CXX_FLAGS}>")\n')
    f.write('# Sets install directory for all the headers in the list\n')
    f.write('foreach(_header ${headers_list})\n')
    f.write('    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/' + libname + ')\n')
    f.write('endforeach()\n')
elif (lang == 'C'):
    f.write('add_library(' + libname + ' OBJECT ${sources_list} ${headers_list})\n')
    f.write('set_target_properties(' + libname + ' PROPERTIES POSITION_INDEPENDENT_CODE 1 )\n')
    f.write('# Sets install directory for all the headers in the list\n')
    f.write('foreach(_header ${headers_list})\n')
    f.write('    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/' + libname + ')\n')
    f.write('endforeach()\n')
else:
    f.write('add_library(' + libname + ' OBJECT ${sources_list})\n')
    f.write('set_target_properties(' + libname + ' PROPERTIES POSITION_INDEPENDENT_CODE 1 )\n')

print('Template for {} created'.format(libname))
print('Don\'t forget to fix excluded files and dependencies!!!')

# vim:et:ts=4:sw=4
