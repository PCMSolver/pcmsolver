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

import sys
import os
import re
from recommonmark.parser import CommonMarkParser
import sphinx_rtd_theme
import subprocess
import shutil

sys.path.append(os.path.abspath(os.curdir))
import cloc_tools

sys.path.append(os.path.abspath('../tools'))
import metadata

# We assume that the sphinx-build command is issued from within the doc
# directory and that the build directory is called _build:
# >> sphinx-build . _build
# This simplifies handling local and RTD build considerably.

# project_root_dir -- the root of the project, i.e. <checkout_dir>/pcmsolver
project_root_dir = os.path.abspath(os.pardir)
# project_doc_dir  -- .rst location, i.e. <checkout_dir>/pcmsolver/doc
project_doc_dir = os.getcwd()
# project_src_dir  -- source code location, i.e. <checkout_dir>/pcmsolver/src
project_src_dir = os.path.join(project_root_dir, 'src')
# Finally, this is where the documentation will be built
doc_build_dir = os.path.join(project_doc_dir, '_build')
print('Project root directory {}'.format(project_root_dir))
print('Project doc directory {}'.format(project_doc_dir))
print('Project src directory {}'.format(project_src_dir))
print('Documentation build directory {}'.format(doc_build_dir))

extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.todo', 'sphinx.ext.coverage', 'sphinx.ext.mathjax', 'sphinx.ext.ifconfig',
    'sphinxcontrib.bibtex', 'breathe'
]

breathe_projects = {'PCMSolver': os.path.join(doc_build_dir, 'xml')}
breathe_default_project = 'PCMSolver'
breathe_default_members = ('members', 'protected-members', 'private-members')
source_parsers = {
    '.md': CommonMarkParser,
}
source_suffix = ['.rst', '.md']
master_doc = 'index'
project = 'PCMSolver'
copyright = '2018, Roberto Di Remigio, Luca Frediani and contributors'
author = 'Roberto Di Remigio, Luca Frediani and contributors'
pcmsolver_version = metadata.__version_long
version = pcmsolver_version
sane_describe = re.compile(
    """^v?(?P<tag>(?P<forwardseries>(?P<major>\d+)\.(?P<minor>\d+))[\.]?(?P<patch>\d+)?[-]?(?P<prere>((a)|(b)|(rc))\d+)?)?[+](?P<sha>\w+)?$"""
)
mobj = sane_describe.match(pcmsolver_version)
major = mobj.group('major')
minor = mobj.group('minor')
patch = mobj.group('patch')
if mobj.group('prere'):
    tweak = '+'.join([mobj.group('prere'), mobj.group('sha')])
else:
    tweak = sha = mobj.group('sha')
language = 'en'
exclude_patterns = ['_build']
pygments_style = 'sphinx'
todo_include_todos = True
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_logo = 'gfx/logo.jpg'
html_use_smartypants = True
html_show_sphinx = True
html_show_copyright = True
html_search_language = 'en'
htmlhelp_basename = 'PCMSolverdoc'


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""
    try:
        retcode = subprocess.call('doxygen {0}'.format(os.path.join(folder, 'Doxyfile')), shell=True)
        if retcode < 0:
            sys.stderr.write('doxygen terminated by signal {0}'.format(-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: {0}".format(e))


def remove(obj):
    if os.path.isdir(obj):
        try:
            shutil.rmtree(obj)
        except OSError:
            pass
    else:
        try:
            os.remove(obj)
        except OSError:
            pass


def configure_file(rep, fname, **kwargs):
    r''' Configure a file.

    :param rep:
       a (placeholder : replacement) dictionary
    :param fname:
       name of the file to be configured, without suffix
    :param \**kwargs:
       See below

    :Keyword arguments:
       * *in_path*  -- directory for the unconfigured file
       * *suffix*   -- suffix of the unconfigured file, with separators
       * *prefix*   -- prefix for the configured file
       * *out_path* -- directory for the configured file
    '''
    in_path = kwargs.get('in_path', os.getcwd())
    suffix = kwargs.get('suffix', '.in')
    out_path = kwargs.get('out_path', in_path)
    prefix = kwargs.get('prefix', '')
    fname_in = fname + suffix
    filedata = ''
    with open(os.path.join(in_path, fname_in), 'r') as fin:
        filedata = fin.read()
    rep = dict((re.escape(k), v) for k, v in rep.items())
    pattern = re.compile("|".join(list(rep.keys())))
    filedata = pattern.sub(lambda m: rep[re.escape(m.group(0))], filedata)
    fname_out = prefix + fname
    with open(os.path.join(out_path, fname_out), 'w+') as fout:
        fout.write(filedata)


def generate_bar_charts(perl_exe, count_dirs, scratch_dir, output_dir):
    """Generate lines-of-code bar charts.

    :param perl_exe:
       the perl executable
    :param count_dirs:
       a list of directories where to count lines of code
    :param scratch_dir:
       where intermediate files (JSON, matplotlib scripts) are to be saved
    :param output_dir:
       where the bar charts are to be saved
    """
    # Create directory to save bar charts SVG
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Generate per-directory and total charts
    for count_dir in count_dirs:
        script = cloc_tools.bar_chart(perl_exe, count_dir, scratch_dir, output_dir)
        exec(compile(open(script).read(), script, 'exec'))
        os.remove(script)
    # Generate total charts
    script = cloc_tools.bar_chart(perl_exe, project_root_dir, scratch_dir, output_dir, is_total=True)
    exec(compile(open(script).read(), script, 'exec'))
    os.remove(script)


def setup(app):
    # Clean up leftovers, these are identified based on the .gitignore in this directory
    print('Clean up leftovers from previous build')
    [remove(os.path.join(project_doc_dir, x.strip())) for x in open(os.path.join(project_doc_dir, '.gitignore'))]
    # Configure Doxyfile.in
    dot_path = os.path.split(which('dot'))[0] if which('dot') else ''
    perl_exe = which('perl')
    rep = {
        '@PROJECT_VERSION_MAJOR@': major,
        '@PROJECT_VERSION_MINOR@': minor,
        '@PROJECT_VERSION_PATCH@': patch,
        '@PROJECT_VERSION_TWEAK@': tweak,
        '@project_root_dir@': project_root_dir,
        '@doc_build_dir@': doc_build_dir,
        '@DOXYGEN_DOT_PATH@': dot_path,
        '@PERL_EXECUTABLE@': perl_exe
    }
    configure_file(rep, 'Doxyfile', in_path=project_doc_dir, suffix='.in')
    # Make a copy of api/pcmsolver.h and strip it of all
    # PCMSolver_API markers in front of function signatures
    rep = {'PCMSolver_API ': '', 'pcmsolver.h': 'mock_pcmsolver.h'}
    configure_file(
        rep,
        'pcmsolver.h',
        in_path=os.path.join(project_root_dir, 'api'),
        suffix='',
        prefix='mock_',
        out_path=doc_build_dir)
    # Generate directories list (full paths), remove bin using filter
    dirs = [os.path.join(r, x) for r, d, _ in os.walk(project_src_dir) for x in d]
    # Filter bin and utils/getkw from from the directories list
    exclude = [os.path.join(project_src_dir, 'bin'), os.path.join(project_src_dir, 'utils', 'getkw')]
    dirs = list(filter(lambda x: x not in exclude, dirs))
    generate_bar_charts(perl_exe, dirs, os.path.join(doc_build_dir), os.path.join(project_doc_dir, 'gfx',
                                                                                  'bar_charts'))

    if (os.environ.get('READTHEDOCS', None) == 'True'):

        def generate_doxygen_xml(app):
            """Run the doxygen make commands on the ReadTheDocs server"""
            run_doxygen(os.getcwd())

        # Add hook for building doxygen xml when needed
        app.connect("builder-inited", generate_doxygen_xml)
    else:
        run_doxygen(project_doc_dir)
