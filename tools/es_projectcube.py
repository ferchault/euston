#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the projection of the cube file contents on a lattice vector.

.. autofunction:: main

Command Line Interface
----------------------

.. program:: es_projectcube.py

.. option:: filename

   The cubefile to read from. Both Bohr and Angstrom units are supported. May be gzipped.

.. option:: index

   Axis vector index to sum over. 0 = first, 2 = last.

.. option:: --absolute

   Whether to sum absolute values or raw values.

.. option:: --pervolume

   Whether to normalise the sum of the data points per volume.

.. option:: --perpoint

   Whether to normalise the sum of the data points per number of points.

Implementation
--------------
"""

# system modules
import argparse

# custom modules
import euston.io as io
import euston.geometry as geom

parser = argparse.ArgumentParser(description='Calculates the projection of the cube file contents on a lattice vector.')
parser.add_argument('filename', type=str, help='The cube file.')
parser.add_argument('output', type=str, help='The output file name.')
parser.add_argument('index', type=int, help='The axis index.')
parser.add_argument('--absolute', action='store_true', help='Whether to sum absolute values or raw values.')
parser.add_argument('--pervolume', action='store_true', help='Give projected value per slice volume.')
parser.add_argument('--perpoint', action='store_true', help='Give projected value per data point in slice.')

def main(parser):
    """
    Main routine wrapper.

    :param argparse.ArgumentParser parser: Argument parser
    """
    args = parser.parse_args()

    if args.index not in range(3):
        raise ValueError('Axes index invalid.')

    print 'Reading cubefile...                 ',
    cube = io.CubeFile(args.filename)
    print 'Completed, %d atoms %d voxels.      ' % (cube.count_atoms(), cube.count_voxels())

    print 'Calculating slice volume...         ',
    h_mat = cube.get_h_matrix()
    abc = geom.hmatrix_to_abc(h_mat)
    if args.index == 0:
        slicecount = cube.get_xlen()
    elif args.index == 1:
        slicecount = cube.get_ylen()
    else:
        slicecount = cube.get_zlen()
    slicevolume = geom.cell_volume(h_mat) / slicecount
    print 'Completed, %d slices of %f Angstrom^3 each.' % (slicecount, slicevolume)

    fh = open(args.output, 'w')
    print 'Summing voxel along axis...         ',
    proj = cube.get_projection(args.index, args.absolute)
    if args.pervolume:
        fh.write('# Normalised by slice volume.\n')
        proj /= slicevolume
    if args.perpoint:
        fh.write('# Normalised by data point count per slice.\n')
        proj /= (cube.count_voxels()/len(proj))
    print 'Completed.'

    fh.write('# index, centered bin position along cell vector, cube data unit per aforementioned units\n')
    for idx, e in enumerate(proj):
        fh.write('%d %f %f\n' % (idx, abc[args.index] / slicecount * (idx+0.5), e))
    fh.close()

if __name__ == '__main__':
    main(parser)
