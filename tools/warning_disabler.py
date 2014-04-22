#! /usr/bin/python

import string
import os
import re

#-------------------------------------------------------------------------------

disabler_start = '''/* warning-disabler-start */

'''
disabler_end   = '''
/* warning-disabler-end */

'''

disabler_top = '''#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#endif
'''

top = disabler_start + disabler_top + disabler_end

disabler_bottom = '''#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif
'''

bottom = disabler_start + disabler_bottom + disabler_end

#-------------------------------------------------------------------------------

def from_file(file):

    f = open(file, 'r')
    s = f.read()
    f.close()
    return s

#-------------------------------------------------------------------------------

def to_file(s, file):

    f = open(file, 'w')
    f.write(s)
    f.close()

#-------------------------------------------------------------------------------

h_file_l = []
c_file_l = []

pwd = os.getcwd()

for path, dir, files in os.walk(pwd):
    for file in files:
        if file[-2:] == '.c':
           c_file_l.append(path + '/' + file)
        if file[-2:] == '.h':
           h_file_l.append(path + '/' + file)

for file in h_file_l:
    s = from_file(file)
    s = top + s + bottom
    to_file(s, file)

for file in c_file_l:
    s = from_file(file)
    s = top + s + bottom
    to_file(s, file)
