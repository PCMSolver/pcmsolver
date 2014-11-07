#!/usr/bin/env python

import os, sys, subprocess, re

OBJDIR = sys.argv[1]
CMAKE_AR = sys.argv[2]
OUTFILE = sys.argv[3]
OBJLISTFILERPATH = sys.argv[4]

os.chdir(OBJDIR)

files = []
for line in open(OBJLISTFILERPATH, 'r').readlines():
    line = re.sub(r'__.SYMDEF SORTED', '', line)
    files.append(line[0:len(line)-1])


files = ' '.join(files)

command = CMAKE_AR + ' ru ' + OUTFILE + ' ' + files

print('Running: ' + command)
subprocess.call(command, shell=True)

exit(0)
