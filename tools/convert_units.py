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
from tempfile import TemporaryFile

toAngstrom =  0.52917721092; # CODATA 2010 recommended data
toAtomicUnits = 1.0/toAngstrom;

def main(argv):
    """ Convert a .txt file containing an array to a .txt file.
    Contents of the array are laid out in memory in Fortran order.
    Unit conversions and transposition of data are performed according
    to the following flags:
    -l convert a length from angstrom to au
    -a convert an area from angstrom**2 to au**2
    -v convert a volume from angstrom**3 to au**3
    -t transpose the array
    -i <input_file.txt>
    -o <output_file.npy>
    """
    convert_length_to_au = False
    convert_area_to_au = False
    convert_volume_to_au = False
    transpose = False

    inputfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv,"hlavti:o:",["l", "a", "v", "t", "ifile=", "ofile="])
    except getopt.GetoptError:
        print('convert_units.py -a -l -v -t -i <input_file> -o <output_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('convert_units.py -a -l -v -t -i <input_file> -o <output_file>')
            sys.exit()
        elif opt == '-l':
            convert_length_to_au = True
        elif opt == '-a':
            convert_area_to_au = True
        elif opt == '-v':
            convert_volume_to_au = True
        elif opt == '-t':
            transpose = True
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
    print('Input file is {}'.format(input_file))
    print('Output file is {}'.format(output_file))

    # Load the .txt file
    contents = numpy.loadtxt(input_file, dtype=float, unpack=transpose)
    if convert_length_to_au:
        contents *= toAtomicUnits
    elif convert_area_to_au:
        contents *= (toAtomicUnits * toAtomicUnits)
    elif convert_volume_to_au:
        contents *= (toAtomicUnits * toAtomicUnits * toAtomicUnits)
    # Reshape to get array into Fortran order
    numpy.reshape(contents, contents.shape, 'F')

    # Now save everything into a .npy file
    numpy.savetxt(output_file, contents)

if __name__ == "__main__":
    main(sys.argv[1:])

# vim:et:ts=4:sw=4
