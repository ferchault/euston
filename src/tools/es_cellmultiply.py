#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Multiplies a cell along the lattice vectors.

Comes with support for scaled input and output, h matrices and arbitrary cell shapes.

Supported File Formats
----------------------

====== ====== ====== =========
Format Input  Output Extension
====== ====== ====== =========
XYZ    Native Native .xyz
====== ====== ====== =========

Command Line Interface
----------------------
.. program:: es_cellmultiply.py 

.. option:: input

   XYZ input filename.

.. option:: output

   XYZ output filename.

.. option:: --X

   X repeat count. Default: 1. Only positive integers supported.

.. option:: --Y

   Y repeat count. Default: 1. Only positive integers supported.   

.. option:: --Z

   Z repeat count. Default: 1. Only positive integers supported.   

.. option:: --sc_in

   Whether the input file contains scaled coordinates.

.. option:: --sc_out

   Whether the out file contains scaled coordinates.

.. option:: --hmat

   H matrix (cell vectors in columns). Comma-separated in row-first notation without spaces.

.. option:: --abc

   a, b, c, alpha, beta, gamma. Comma-separated without spaces.

.. option:: --radians

   Whether angles in --abc are given in radians.

Implementation
--------------
"""

import argparse
import numpy as np
import euston.io as io
import euston.geometry as geo

parser = argparse.ArgumentParser(description='Multiplies a cell along the lattice vectors.')
parser.add_argument('input', type=str, help='XYZ input filename')
parser.add_argument('output', type=str, help='XYZ output filename')
parser.add_argument('--X', type=int, help='X repeat', default=1)
parser.add_argument('--Y', type=int, help='Y repeat', default=1)
parser.add_argument('--Z', type=int, help='Z repeat', default=1)
parser.add_argument('--sc_in', action='store_true', help='Whether the input file contains scaled coordinates.')
parser.add_argument('--sc_out', action='store_true', help='Whether the output file should contain scaled coordinates.')
parser.add_argument('--hmat', type=str, help='H matrix (cell vectors in columns). Comma-separated in row-first notation without spaces.')
parser.add_argument('--abc', type=str, help='a, b, c, alpha, beta, gamma. Comma-separated without spaces.')
parser.add_argument('--radians', action='store_true', help='Whether angles in --abc are given in radians.')

def main(parser):
    """
    Main routine wrapper.

    :param argparse.ArgumentParser parser: Argument parser
    """
    args = parser.parse_args()

    xyz = io.XYZ(args.input)

    hmat = None
    if not (args.sc_in and args.sc_out):
        if (args.hmat is None and args.abc is None) or (args.hmat is not None and args.abc is not None):
            print 'Please specify either the H matrix or lattice constants unless input and output are scaled.'
            exit(1)

        if args.hmat is not None:
            try:
                hmat = map(float, args.hmat.split(','))
            except:
                print 'Invalid H matrix entries.'
                exit(2)
            hmat = np.array(hmat).reshape((3,3)).T

        if args.abc is not None:
            try:
                abc = map(float, args.abc.split(','))
            except:
                print 'Invalid abc entries.'
                exit(2)
            hmat = geo.abc_to_hmatrix(*abc, degrees=(not args.radians))

    for images in (args.X, args.Y, args.Z):
        if images < 1:
            print 'Invalid multiplier setting - has to be at least one.'
            exit(3)

    multiplied = geo.cell_multiply(xyz.get_coordinates(), args.X, args.Y, args.Z, h_matrix=hmat, scaling_in=args.sc_in, scaling_out=args.sc_out)
    factor = args.X*args.Y*args.Z

    output = io.XYZ()
    output.set_data(xyz.get_labels()*factor, multiplied)
    io.write_lines(args.output, output.to_string())

if __name__ == '__main__':
    main(parser)
