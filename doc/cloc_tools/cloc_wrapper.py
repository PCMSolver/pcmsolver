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

# Execute cloc.pl Perl script and wrap results in a nicer format.
# cloc script by Al Danial, available at http://cloc.sourceforge.net/
# licensed under the GNU General Public License v2
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import subprocess
import os
import time
import sys
import json
import shlex


def cloc_command(perl, cloc_script, args):
    """Wrapper to the cloc.pl Perl script.

    Keyword arguments:
    perl -- perl executable
    cloc_script -- the cloc.pl script
    files -- list of files to parse
    args -- additional list of command line arguments to cloc.pl
    """
    command = [perl, cloc_script] + args
    try:
        retcode = subprocess.call(command, shell=False)
        if retcode < 0:
            sys.stderr.write('{0} terminated by signal {1}'.format(command, -retcode))
    except OSError as e:
        sys.stderr.write("{0} execution failed: {1}".format(command, e))


def bar_chart(perl, count_dir, scratch_dir, output_dir, is_total=False):
    """Generates matplotlib script for lines of code bar chart.

    Keyword arguments:
    perl -- perl executable
    count_dir -- directory where to count lines of code
    scratch_dir -- where intermediate files (JSON, matplotlib scripts) are to be saved
    output_dir -- where the bar charts are to be saved
    is_total -- whether this is the global counting or not
    """
    if is_total:
        template = os.path.join(os.path.dirname(__file__), 'total_bar_chart.txt')
        with open(template, 'r') as tmp:
            script = '\n' + tmp.read()
        annotation = '\'PCMSolver\\nFiles: {0:d}\'.format(nr_files)'
    else:
        template = os.path.join(os.path.dirname(__file__), 'bar_chart.txt')
        with open(template, 'r') as tmp:
            script = '\n' + tmp.read()
        annotation = '\'Folder: {0}\\nLanguage: {1}\\nFiles: {2:d}\'.format(tag, language, nr_files)'
    tag = os.path.basename(count_dir)
    cloc_data = count_lines_of_code(perl, count_dir, scratch_dir, tag, output_dir)
    plot_script = os.path.join(scratch_dir, tag + '.py')
    with open(plot_script, 'w+') as bar_plot:
        bar_plot.write(header())
        bar_plot.write(script.format(svg_dir=output_dir, tag=tag, data=cloc_data, annotation=annotation))
    return plot_script


def count_lines_of_code(perl, count_dir, scratch_dir, plot_name, svg_dir):
    """Runs cloc.pl to get total count and returns data for matplotlib script.

    Keyword arguments:
    perl -- perl executable
    count_dir -- directory where to run cloc.pl
    scratch_dir -- where intermediate files (JSON and matplotlib script) are to be saved
    plot_name -- name of the plot (e.g. cavity or solver)
    svg_dir -- where the SVG will be saved
    """
    # yapf: disable
    json_report = os.path.join(scratch_dir, 'cloc_temp.json')
    cloc_options = [
        count_dir
      , '--json'
      , '--report-file={}'.format(json_report)
      , '--include-lang="C++","C","Fortran 77","Fortran 90","Fortran 95","Python"'
      , '--exclude-dir="external,examples,doc,cmake"'
      , '--fullpath', '--not-match-d="src/utils/getkw"'
      , '--fullpath', '--not-match-d=\'/build\S[a-zA-Z0-9]+/\''
      , '--force-lang="C++",hpp'
      , '--force-lang="C",h'
      , '--force-lang="Fortran 90",f'
      , '--force-lang="Fortran 90",F'
      , '--force-lang="Fortran 90",f95'
      , '--force-lang="Fortran 90",F95'
      , '--quiet'
    ]
    # yapf: enable
    cloc_command(perl, os.path.join(os.path.dirname(__file__), 'cloc.pl'), shlex.split(' '.join(cloc_options)))
    with open(json_report, 'r') as cloc_out:
        data = json.load(cloc_out)
    os.remove(json_report)
    # Prepare data for bar chart
    cloc_data = {
        'nr_blanks': data['SUM']['blank'],
        'nr_comments': data['SUM']['comment'],
        'nr_files': data['SUM']['nFiles'],
        'nr_cpp_code': data['C++']['code'] if 'C++' in data else 0,
        'nr_c_code': data['C']['code'] if 'C' in data else 0,
        'nr_fortran_code': data['Fortran 90']['code'] if 'Fortran 90' in data else 0,
        'nr_python_code': data['Python']['code'] if 'Python' in data else 0
    }
    return cloc_data


def header():
    """Print header of matplotlib script generating bar chart.

    Keyword arguments:
    plot_script -- name of the plotting script
    """

    header = """#!/usr/bin/env python
# Automatically generated on {now}
# Data obtained from the cloc.pl Perl script.
# cloc script by Al Danial, available at http://cloc.sourceforge.net/
# licensed under the GNU General Public License v2
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License
    """.format(now=time.strftime('%c'))
    return header
