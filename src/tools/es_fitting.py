#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the Wigner-Seitz projection of cube file data of periodic data.

.. autofunction:: main

Command Line Interface
----------------------

.. program:: es_fitting.py 

.. option:: filename

   The cubefile to read from. Both Bohr and Angstrom units are supported.

.. option:: --periodic

   Whether to treat input file in periodic boundary conditions. In the current implementation, this requires copying all atom positions for the 26 direct neighbouring cubes and therefore carries substantial cost in terms of memory requirements. The runtime scales O(log n) with the number of atoms.

.. option:: --leafsize

   Number of points at which brute-force nearest neighbour search is employed. Increasing this number introduces higher memory requirements, but may speed up processing. Default: 5000.

Implementation
--------------
"""

# system modules
import argparse

# third-party modules
from scipy.spatial import KDTree
import numpy as np

# custom modules
import euston.io as io
import euston.geometry as geom

parser = argparse.ArgumentParser(description='Calculates the Wigner-Seitz projection of cube file data of periodic data.')
parser.add_argument('filename', type=str, help='The cube file.')
parser.add_argument('--periodic', action='store_true', help='Treats cube data periodically.')
parser.add_argument('--leafsize', type=int, help='Number of points at which brute-force nearest neighbour search is employed.', default=5000)

def main(parser):
    """
    Main routine wrapper.

    :param argparse.ArgumentParser parser: Argument parser
    """
    args = parser.parse_args()

    print 'Reading cubefile...                 ',
    cube = io.CubeFile(args.filename)
    print 'Completed, %d atoms %d voxels.' % (cube.count_atoms(), cube.count_voxels())

    print 'Calculating cell dimensions...      ',
    h_mat = cube.get_h_matrix()
    print 'Completed.'

    print 'Adding image atoms...               ',
    coord = cube.get_coordinates()
    images = np.zeros(((len(coord)*27), 3))
    counter = 0
    for x in range(-1, 2):
        for y in range(-1, 2):
            for z in range(-1, 2):
                if x*y*z == 0:
                    continue
                shift = (h_mat.transpose()*np.array([x, y, z])).sum(axis=0)
                images[counter*cube.count_atoms():(counter+1)*cube.count_atoms(), :] = coord+shift
                counter += 1
    print 'Completed, %d atoms added.' % (cube.count_atoms()*26)

    images[-cube.count_atoms():] = coord

    if args.periodic:
        kdt = KDTree(images, leafsize=args.leafsize)
    else:
        kdt = KDTree(coord, leafsize=args.leafsize)
    vals = np.zeros(cube.count_atoms())
    for x in range(cube.get_xlen()):
        print x
        for y in range(cube.get_ylen()):
            for z in range(cube.get_zlen()):
                d, i = kdt.query(cube.get_voxel_pos(x, y, z, centered=True))
                i = i % cube.count_atoms()
                vals[i] += cube.get_val(x, y, z)

    print 'Results (atom - value)'
    vals *= cube.get_voxel_volume()
    for idx, val in enumerate(vals):
        print idx, val

if __name__ == '__main__':
    main(parser)
