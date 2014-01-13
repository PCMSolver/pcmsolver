#!/usr/bin/env python
# vim:ft=python
# Read a .txt file containing arrays of doubles into
# a NumPy object. The script then writes the contents
# of the .txt file to a binary .npy file.
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import sys
import getopt
import os
import string
import numpy


#numpy.loadtxt("filename")

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'txt_to_npy.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'txt_to_npy.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print 'Input file is "', inputfile
   print 'Output file is "', outputfile

   array = numpy.loadtxt(inputfile, dtype=float)
   print(array)

   numpy.save(outputfile, array)


if __name__ == "__main__":
   main(sys.argv[1:])


# vim:et:ts=4:sw=4
