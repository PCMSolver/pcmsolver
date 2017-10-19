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


def glob_sources(dname):
    import glob
    import os
    """Create a sorted list of test sources to be used in a CMakeLists.txt file."""
    return sorted(glob.glob(os.path.join(dname, '*.cpp')))


def extract_labels(dname, srcs):
    import os
    import re
    """Extract Catch labels from unit test source file."""
    return [re.findall('\[(.*?)\]', content) for content in (open(os.path.join(dname, src)).read() for src in sources)]


dname = os.getcwd()
fname = os.path.join(dname, 'CMakeLists.txt.try')
f = open(fname, 'w')
sources = glob_sources(dname)
labels = extract_labels(dname, sources)
src_lbl = dict(zip([os.path.basename(src) for src in sources], labels))

dirname = os.path.basename(os.path.normpath(dname))
message = 'add_library(' + dirname + '-tests OBJECT '
message += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sources)
message += ')\n'
message += 'if(BUILD_CUSTOM_BOOST)\n'
message += '  add_dependencies(' + dirname + '-tests custom_boost)\n'
message += 'endif()\n\n'
f.write(message)

for src, lbl in src_lbl.iteritems():
    f.write('# %s test\n' % src)
    f.write('add_Catch_test(%s "%s")\n\n' % (os.path.splitext(src)[0], ';'.join(lbl)))

print('Template created')
print('Don\'t forget to fix labels, excluded files and dependencies!!!')

# vim:et:ts=4:sw=4
