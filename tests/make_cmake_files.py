#!/usr/bin/python
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
labels  = extract_labels(dname, sources)
src_lbl = dict(zip([os.path.basename(src) for src in sources], labels))
#message = 'list(APPEND tests_sources '
#message += ' '.join('%s' % ''.join(map(str, os.path.basename(x))) for x in sources)
#message += ')\n\n'
#f.write(message)

for src, lbl in src_lbl.iteritems():
    f.write('# %s test\n' % src)
    f.write('set_property(GLOBAL APPEND PROPERTY UnitTestsSources ${CMAKE_CURRENT_LIST_DIR}/%s)\n' % src)
    f.write('add_Catch_test(%s "%s")\n\n' % (os.path.splitext(src)[0], ';'.join(lbl)))

print('Template created')
print('Don\'t forget to fix labels, excluded files and dependencies!!!')

# vim:et:ts=4:sw=4
