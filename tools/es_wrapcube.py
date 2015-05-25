#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wraps all coordinates in a cube file back into the system cell.

.. autofunction:: main

Command Line Interface
----------------------

.. program:: es_wrapcube.py

.. option:: filename

   The cubefile to read from. Both Bohr and Angstrom units are supported. May be gzipped.

.. option:: output

   The cubefile to save to.


Implementation
--------------
"""

# system modules
import argparse
import numpy as np

# custom modules
import euston.io as io
import euston.geometry as geom

parser = argparse.ArgumentParser(description='Wraps all coordinates in a cube file back into the system cell.')
parser.add_argument('filename', type=str, help='The cube file.')
parser.add_argument('output', type=str, help='The output file name.')

def main(parser):
    """
    Main routine wrapper.

    :param argparse.ArgumentParser parser: Argument parser
    """
    args = parser.parse_args()

    print 'Reading cubefile...                 ',
    cube = io.CubeFile(args.filename)
    print 'Completed, %d atoms %d voxels.      ' % (cube.count_atoms(), cube.count_voxels())

    print 'Wrapping atom coordinates..         ',
    h_mat = cube.get_h_matrix()
    scaled = geom.cartesian_to_scaled_coordinates(cube.get_coordinates(), h_mat)
    scaled = np.mod(scaled, 1)
    coord = geom.scaled_to_cartesian_coordinates(scaled, h_mat)
    cube.set_coordinates(coord)
    print 'Completed.'

    print 'Writing data to new cube file...    ',
    io.write_lines(args.output, cube.to_string())
    print 'Completed.'

if __name__ == '__main__':
    main(parser)
