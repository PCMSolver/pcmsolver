#! /usr/bin/python

import string
import os
import re

#-------------------------------------------------------------------------------

pcmsolver_copyright_fortran_start = '!pcmsolver_copyright_start'
pcmsolver_copyright_fortran_end   = '!pcmsolver_copyright_end'

pcmsolver_copyright_fortran = '''
!      PCMSolver, an API for the Polarizable Continuum Model
!      Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
!      
!      This file is part of PCMSolver.
!
!      PCMSolver is free software: you can redistribute it and/or modify       
!      it under the terms of the GNU Lesser General Public License as published by
!      the Free Software Foundation, either version 3 of the License, or
!      (at your option) any later version.
!                                                                           
!      PCMSolver is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!      GNU Lesser General Public License for more details.
!                                                                           
!      You should have received a copy of the GNU Lesser General Public License
!      along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
!
!      For information on the complete list of contributors to the
!      PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
'''

pcmsolver_copyright_c_start = '/\* pcmsolver_copyright_start \*/'
pcmsolver_copyright_c_end   = '/\* pcmsolver_copyright_end \*/'

pcmsolver_copyright_c = '''
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *                                                                          
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *                                                                          
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
'''

#-------------------------------------------------------------------------------

def from_file(file):

    f = open(file, 'r')
    s = f.read()
    f.close()
    return s

#-------------------------------------------------------------------------------

def to_file(s, file):

    f = open(file, 'w')
    f.write(s)
    f.close()

#-------------------------------------------------------------------------------

fortran_file_l = []
c_file_l       = []
h_file_l       = []
cpp_file_l     = []
hpp_file_l     = []

pwd = os.getcwd()

for path, dir, files in os.walk(pwd):
    for file in files:
        if file[-4:] == '.F90':
           fortran_file_l.append(path + '/' + file)
        if file[-4:] == '.f90':
           fortran_file_l.append(path + '/' + file)
        if file[-2:] == '.F':
           fortran_file_l.append(path + '/' + file)
        if file[-2:] == '.c':
           c_file_l.append(path + '/' + file)
        if file[-2:] == '.h':
           h_file_l.append(path + '/' + file)
        if file[-4:] == '.cpp':
           cpp_file_l.append(path + '/' + file)
        if file[-4:] == '.hpp':
           hpp_file_l.append(path + '/' + file)

# fortran files pcmsolver copyright
p = re.compile('(?<=' + pcmsolver_copyright_fortran_start + ').*(?=' + pcmsolver_copyright_fortran_end + ')', re.DOTALL)
for file in fortran_file_l:
    s = from_file(file)
    if pcmsolver_copyright_fortran_start in s:
       s = re.sub(p, pcmsolver_copyright_fortran, s)
       to_file(s, file)

# c files pcmsolver copyright
p = re.compile('(?<=' + pcmsolver_copyright_c_start + ').*(?=' + pcmsolver_copyright_c_end + ')', re.DOTALL)
for file in c_file_l:
    s = from_file(file)
    if string.replace(pcmsolver_copyright_c_start, '\\', '') in s:
       s = re.sub(p, pcmsolver_copyright_c, s)
       to_file(s, file)

# h files pcmsolver copyright
p = re.compile('(?<=' + pcmsolver_copyright_c_start + ').*(?=' + pcmsolver_copyright_c_end + ')', re.DOTALL)
for file in h_file_l:
    s = from_file(file)
    if string.replace(pcmsolver_copyright_c_start, '\\', '') in s:
       s = re.sub(p, pcmsolver_copyright_c, s)
       to_file(s, file)

# cpp files pcmsolver copyright
p = re.compile('(?<=' + pcmsolver_copyright_c_start + ').*(?=' + pcmsolver_copyright_c_end + ')', re.DOTALL)
for file in cpp_file_l:
    s = from_file(file)
    if string.replace(pcmsolver_copyright_c_start, '\\', '') in s:
       s = re.sub(p, pcmsolver_copyright_c, s)
       to_file(s, file)

# hpp files pcmsolver copyright
p = re.compile('(?<=' + pcmsolver_copyright_c_start + ').*(?=' + pcmsolver_copyright_c_end + ')', re.DOTALL)
for file in hpp_file_l:
    s = from_file(file)
    if string.replace(pcmsolver_copyright_c_start, '\\', '') in s:
       s = re.sub(p, pcmsolver_copyright_c, s)
       to_file(s, file)
