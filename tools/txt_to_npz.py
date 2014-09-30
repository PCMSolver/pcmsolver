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

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hw:e:c:n:o:",["wfile=", "efile=", "cfile=", "nfile", "ofile="])
    except getopt.GetoptError:
        print('txt_to_npy.py -w <weights_file> -e <elRadius_file> -c <centers_file> -n <normals_file> -o <output_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('txt_to_npy.py -w <weights_file> -e <elRadius_file> -c <centers_file> -n <normals_file> -o <output_file>')
            sys.exit()
        elif opt in ("-w", "--wfile"):
            weights_file = arg
        elif opt in ("-e", "--efile"):
            elRadius_file = arg
        elif opt in ("-c", "--cfile"):
            centers_file = arg
        elif opt in ("-n", "--nfile"):
            normals_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
    print('Weights file is {}'.format(weights_file))
    print('Element radius file is {}'.format(elRadius_file))
    print('Centers file is {}'.format(centers_file))
    print('Normals file is {}'.format(normals_file))
    print('Output file is {}'.format(output_file))

    # First of all, get the number of elements
    raw_w = numpy.loadtxt(weights_file, dtype=numpy.intc)
    elements = numpy.array(raw_w[0], ndmin=1)

    # Load the weights and perform sanity check.
    weights = numpy.loadtxt(weights_file, dtype=float, skiprows=1)
    if (weights.size != elements):
        print('Data saved in weights array is not consistent!')
        sys.exit(1)
    else:
        # Reshape to get array into Fortran order
        numpy.reshape(weights, (elements, 1), 'F')

    # Load the element radius and perform sanity check.
    elRadius = numpy.loadtxt(elRadius_file, dtype=float)
    if (elRadius.size != elements):
        print('Data saved in element radius array is not consistent!')
        sys.exit(1)
    else:
        # Reshape to get array into Fortran order
        numpy.reshape(elRadius, (elements, 1), 'F')

    # Load the centers and perform sanity check.
    centers = numpy.loadtxt(centers_file, dtype=float)
    if (centers.size != 3*elements):
        print('Data saved in centers array is not consistent!')
        sys.exit(1)
    else:
        # Reshape to get array into Fortran order
        numpy.reshape(centers, (3, elements), 'F')

    # Load the normals and perform sanity check.
    normals = numpy.loadtxt(normals_file, dtype=float)
    if (normals.size != 3*elements):
        print('Data saved in normals array is not consistent!')
        sys.exit(1)
    else:
        # Reshape to get array into Fortran order
        numpy.reshape(normals, (3, elements), 'F')

    # Now save everything into a single .npz file
    numpy.savez(output_file, elements=elements, weights=weights, elRadius=elRadius, centers=centers, normals=normals)

if __name__ == "__main__":
    main(sys.argv[1:])

# vim:et:ts=4:sw=4
