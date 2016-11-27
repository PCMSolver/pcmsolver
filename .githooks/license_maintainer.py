#!/usr/bin/python
"""
Add or update the license information in the header of source files.

The script reads in the .gitattributes file, located in the project root
directory, to figure out which files need to be inspected and which license
header they need to have.
For example:

src/pedra/pedra_dlapack.F90 !licensefile
src/solver/*.hpp licensefile=.githooks/LICENSE-C++

The first line specifies that the file in src/pedra/pedra_dlapack.F90 should
not be touched, while the second line states that all .hpp files in src/solver
should get an header from the template in .githooks/LICENSE-C++
Location of files in .gitattributes are always specified with respect
to the project root directory.

The script reads in the appropriate license header template and prepares an
header with the correct year and authors information. These are encoded in the
variables YEAR and AUTHORS. The latter has to be modified by hand.
"""

from datetime import date
import glob
import mmap
import os
import re
import shutil
import tempfile

def add_header(filepath, header, YEAR, AUTHORS):
    """
    Add or update header in source file
    """
    tmpdir = tempfile.gettempdir()
    tmpfil = os.path.join(tmpdir, os.path.basename(filepath) + '.bak')
    shutil.copy2(filepath, tmpfil)
    with open(tmpfil, 'r+b') as tmp:
        inpt = tmp.readlines()
        output = []

        # Check if header is already present
        present = re.compile('PCMSolver, an API for the Polarizable Continuum Model')
        if filter(present.search, inpt):
            # Check if year and authors in current file are up to date
            toupdate = re.compile(r'{0} (?!{1} {2}).*\n'.format('Copyright \(C\)', YEAR, AUTHORS))
            if filter(toupdate.search, inpt):
                print('Updating header in {}'.format(filepath))
                # Check to preserve '#!' at the top of the file
                if len(inpt) > 0 and inpt[0].startswith('#!'):
                    output.append(inpt[0] + '\n')
                    inpt = inpt[1:]
                regex = re.compile(r'Copyright \(C\).*\n')
                repl  = r'Copyright (C) ' + YEAR + ' ' + AUTHORS + '\n'
                output.extend(map(lambda x: re.sub(regex, repl, x), inpt))
        else:
            print('Adding header in {}'.format(filepath))
            # Check to preserve '#!' at the top of the file
            if len(inpt) > 0 and inpt[0].startswith('#!'):
                output.append(inpt[0] + '\n')
                inpt = inpt[1:]
            output.append(header)
            for line in inpt:
                output.append(line)

        if output:
            try:
                f = open(filepath, 'w')
                f.writelines(output)
            except IOError, err:
                print('Something went wrong trying to add header to {}: {}'.format(filepath, err))
            finally:
                f.close()
        os.remove(tmpfil)


def prepare_header(stub, YEAR, AUTHORS):
    """
    Update year and author information in license header template
    """
    with open(stub, 'r+b') as l:
        header = l.read()
        # Insert correct YEAR and AUTHORS in stub
        rep = {'YEAR' : YEAR, 'AUTHORS' : AUTHORS}
        rep = dict((re.escape(k), v) for k, v in rep.iteritems())
        pattern = re.compile("|".join(rep.keys()))
        header = pattern.sub(lambda m: rep[re.escape(m.group(0))], header)
    return header


def file_license(attributes):
    """
    Obtain dictionary { file : license } from .gitattributes
    """
    file_license = {}
    with open(attributes, 'r+b') as f:
        # Read in .gitattributes
        tmp_mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        # Removing all comment lines
        gitattributes = re.sub(r'(?m)^\#.*\n?', '', tmp_mm).split()
        # Obtain list of files
        fil = [x for x in gitattributes if not 'licensefile' in x]
        # Remove licensefile= from strings
        lic = [re.sub(r'licensefile\=', '', x) for x in gitattributes if 'licensefile' in x]
        # Create list of blacklisted files
        blacklist = [fname for key, value in dict(zip(fil, lic)).items()
                       if value == '!licensefile'
                       for fname in glob.glob(key)]
        # Now create a dictionary with the files to be considered for
        # license header manipulation
        file_license = { key : value
                  for k, value in dict(zip(fil, lic)).items()
                  for key in glob.glob(k)
                  if key not in blacklist}
    return file_license


def license_maintainer():
    """
    Maintain license header in source files
    """
    YEAR    = str(date.today().year)
    AUTHORS = 'Roberto Di Remigio, Luca Frediani and collaborators.'

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    headerize = file_license(os.path.join(project_root_dir, '.gitattributes'))

    for fname, license in headerize.items():
        # Prepare header
        header = prepare_header(os.path.join(project_root_dir, license), YEAR, AUTHORS)
        add_header(fname, header, YEAR, AUTHORS)


if __name__ == '__main__':
    license_maintainer()
