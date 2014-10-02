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
    """ Get normal vectors to tesserae representative points.
    Inputs: .npy file containing the coordinates of the representative points
    .npy file containing the coordinates of the centers of the spheres to which the tessera belongs.
    Outputs: .txt file containing the normal vectors
    """

    try:
        opts, args = getopt.getopt(argv,"hc:s:o:",["cfile=", "sfile=", "ofile="])
    except getopt.GetoptError:
        print('get_normals.py -c <centers_file.npy> -s <element_sphere_center_file.npy> -o <normals_file.txt>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('get_normals.py -c <centers_file.npy> -s <element_sphere_center_file.npy> -o <normals_file.txt>')
            sys.exit()
        elif opt in ("-c", "--cfile"):
            centers_file = arg
        elif opt in ("-s", "--sfile"):
            element_sphere_center_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
    print('Centers file is {}'.format(centers_file))
    print('Element sphere center file is {}'.format(element_sphere_center_file))
    print('Output file is {}'.format(output_file))

    # Load the .npy file
    centers = numpy.load(centers_file)
    element_sphere_center = numpy.load(element_sphere_center_file)
    # tmp contains unnormalized normal vectors
    tmp = centers - element_sphere_center
    # normalize
    for i in range(centers.shape[1]):
        tmp[:, i] /= numpy.linalg.norm(tmp[:, i])
    # transfer to the proper array
    normals = tmp
    # Reshape to get array into Fortran order
    numpy.reshape(normals, normals.shape, 'F')

    # Now save everything into a .npy file
    numpy.savetxt(output_file, normals)

if __name__ == "__main__":
    main(sys.argv[1:])

# vim:et:ts=4:sw=4
