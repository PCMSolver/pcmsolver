#!/usr/bin/env python

import sys, subprocess, re

#! path to archiver
CMAKE_AR = sys.argv[1]
#! library to be created
OUTFILE = sys.argv[2]
#! where is file with list of object files
#! to be included in library
OBJLISTFILERPATH = sys.argv[3]

#! create list of files from file names stored at OBJLISTFILERPATH
#! and on some systems remove __.SYMDEF SORTED
files = []
for line in open(OBJLISTFILERPATH, 'r').readlines():
        line = re.sub(r'__.SYMDEF SORTED', '', line)
        files.append(line[0:len(line)-1])

#! print message to console
print('Running: ' + CMAKE_AR + ' ru ' + OUTFILE + ' ' + ' '.join(files))

#! every file in the list will be added to library
for file in files:
    subprocess.call([CMAKE_AR, 'ru', OUTFILE, file])
