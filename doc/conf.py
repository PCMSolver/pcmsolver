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

import sys
import os
import re
from recommonmark.parser import CommonMarkParser
import sphinx_rtd_theme
import subprocess
import shutil
from pygments.lexers import get_lexer_for_filename

sys.path.insert(0, os.path.abspath('../tools'))
import metadata

extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.todo', 'sphinx.ext.coverage', 'sphinx.ext.mathjax', 'sphinx.ext.ifconfig',
    'sphinxcontrib.bibtex', 'breathe'
]

breathe_projects = {'PCMSolver': 'xml'}
breathe_default_project = 'PCMSolver'
breathe_default_members = ('members', 'protected-members', 'private-members')
templates_path = ['_templates']
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
prere = mobj.group('prere')
sha = mobj.group('sha')
tweak = '+'.join([prere, sha])
language = 'en'
exclude_patterns = ['_build']
pygments_style = 'sphinx'
todo_include_todos = True
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_logo = 'gfx/logo.jpg'
html_static_path = ['_static']
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


def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

    if read_the_docs_build:
        run_doxygen(os.getcwd())


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
    import os
    import re
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
    f = open(os.path.join(in_path, fname_in), 'r')
    filedata = f.read()
    f.close()
    rep = dict((re.escape(k), v) for k, v in rep.items())
    pattern = re.compile("|".join(list(rep.keys())))
    filedata = pattern.sub(lambda m: rep[re.escape(m.group(0))], filedata)
    fname_out = prefix + fname
    f = open(os.path.join(out_path, fname_out), 'w+')
    f.write(filedata)
    f.close()


def generate_bar_charts(mod_dir, dir_lang, savedir):
    r'''Generate lines-of-code bar charts.

    :param mod_dir:
       location of the cloc_tools module
    :param dir_lang:
       a (directory : language) dictionary
    :param savedir:
       location of the YAML files
    '''
    import sys
    sys.path.append(mod_dir)
    from cloc_tools import bar_chart
    # Generate scripts and list of scripts (absolute paths)
    list_of_scripts = [bar_chart(root_dir, language, savedir) for root_dir, language in dir_lang.items()]
    # Generate charts
    for fname in list_of_scripts:
        exec(compile(open(fname).read(), fname, 'exec'))


def setup(app):
    # We first need to define some directories:
    # project_root_dir -- the root of the project
    # project_src_dir  -- source code location: os.path.join(project_root_dir, 'src')
    # project_doc_dir  -- .rst location: os.path.join(project_root_dir, 'doc')
    if (os.environ.get('READTHEDOCS', None) == 'True'):
        project_root_dir = os.path.abspath(os.pardir)
        project_doc_dir = os.getcwd()
        project_src_dir = os.path.join(project_root_dir, 'src')
    else:
        project_root_dir = os.getcwd()
        project_doc_dir = os.path.join(project_root_dir, 'doc')
        project_src_dir = os.path.join(project_root_dir, 'src')
    print(('Project root directory {}'.format(project_root_dir)))
    print(('Project doc directory {}'.format(project_doc_dir)))
    print(('Project src directory {}'.format(project_src_dir)))

    # Clean up leftovers
    print('Clean up leftovers from previous build')
    [remove(os.path.join(project_doc_dir, x.strip())) for x in open(os.path.join(project_doc_dir, '.gitignore'))]
    # Configure Doxyfile.in
    dot_path = os.path.split(which('dot'))[0] if which('dot') else ''
    rep = {
        '@PROJECT_VERSION_MAJOR@': major,
        '@PROJECT_VERSION_MINOR@': minor,
        '@PROJECT_VERSION_PATCH@': patch,
        '@PROJECT_VERSION_TWEAK@': tweak,
        '@PROJECT_SOURCE_DIR@': project_root_dir,
        '@DOXYGEN_DOT_PATH@': dot_path,
        '@PERL_EXECUTABLE@': which('perl')
    }
    configure_file(rep, 'Doxyfile', in_path=project_doc_dir, suffix='.in')
    # Make a copy of api/pcmsolver.h and strip it of all
    # PCMSolver_API markers in front of function signatures
    rep = {'PCMSolver_API ': ''}
    configure_file(
        rep,
        'pcmsolver.h',
        in_path=os.path.join(project_root_dir, 'api'),
        suffix='',
        prefix='mock_',
        out_path=project_doc_dir)
    # Configure cloc_tools.py.in
    rep = {
        '@PYTHON_EXECUTABLE@': sys.executable,
        '@PROJECT_SOURCE_DIR@': project_root_dir,
        '@PROJECT_BINARY_DIR@': project_root_dir,
        '@PERL_EXECUTABLE@': which('perl')
    }
    configure_file(
        rep, 'cloc_tools.py', in_path=os.path.join(project_doc_dir, 'gfx'), suffix='.in', out_path=project_doc_dir)
    # Generate directories list (full paths), remove bin using filter
    d = [
        y for y in [os.path.join(root, x) for root, dirs, _ in os.walk(project_src_dir) for x in dirs]
        if y != os.path.join(project_src_dir, 'bin')
    ]
    # Remove 'CMakeLists.txt' from sublists using filter
    f = [[
        z for z in
        [y for y in [x for x in os.listdir(l) if os.path.isfile(os.path.join(l, x))] if y != 'CMakeLists.txt']
        if not z.endswith('.mod')
    ] for l in d]
    # Take first element in each sublist
    f = [x[0] for x in f]
    # Apply map to get language name
    l = [get_lexer_for_filename(x).name for x in f]
    # Finally zip d and f into the dir_lang dictionary
    dir_lang = dict(list(zip(d, l)))
    generate_bar_charts(project_doc_dir, dir_lang, project_doc_dir)

    if (os.environ.get('READTHEDOCS', None) == 'True'):
        # Add hook for building doxygen xml when needed
        app.connect("builder-inited", generate_doxygen_xml)
    else:
        run_doxygen(project_doc_dir)
