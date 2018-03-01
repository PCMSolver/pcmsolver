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

# Create CMakeLists.txt template for leaf directories
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import glob
import os
import re


def glob_sources(dname):
    """Create a sorted list of test sources to be used in a CMakeLists.txt file."""
    return [os.path.basename(x) for x in sorted(glob.glob(os.path.join(dname, '*.cpp')))]


def extract_labels(dname, srcs):
    """Extract Catch labels from unit test source file."""
    return [re.findall('\[(.*?)\]', content) for content in (open(os.path.join(dname, src)).read() for src in sources)]


dname = os.getcwd()
fname = os.path.join(dname, 'CMakeLists.txt.try')
sources = glob_sources(dname)
labels = extract_labels(dname, sources)
src_lbl = dict(zip(sources, labels))
dirname = os.path.basename(os.path.normpath(dname))
message = """add_library({0}-tests OBJECT
  {1}
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies({0}-tests custom_boost)
endif()
""".format(dirname, '\n  '.join(sources))
add_test = []
for src, lbl in src_lbl.items():
    add_test.append("""# {0} test
add_Catch_test(
  NAME
    {1}
  LABELS
    {2}
  )\n""".format(src, os.path.splitext(src)[0], '\n'.join(lbl)))
with open(fname, 'w') as f:
    f.write(message)
    f.write('\n'.join(add_test))

print('Template created')
print('Don\'t forget to fix labels, excluded files and dependencies!!!')
