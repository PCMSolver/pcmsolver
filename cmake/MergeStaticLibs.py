#!/usr/bin/env python

import os, sys, subprocess, re

OBJDIR = sys.argv[1]
CMAKE_COMMAND = sys.argv[2]
CMAKE_AR = sys.argv[3]
OUTFILE = sys.argv[4]
OBJLISTFILERPATH = sys.argv[5]

os.chdir(OBJDIR)

files = []
for line in open(OBJLISTFILERPATH, 'r').readlines():
    line = re.sub(r'__.SYMDEF SORTED', '', line)
    files.append(line[0:len(line)-1])


files = ' '.join(files)
command = CMAKE_AR + ' ru ' + OUTFILE + ' ' + files

subprocess.call(CMAKE_COMMAND + ' -E echo ' + command)
subprocess.call(command)

exit(0)
