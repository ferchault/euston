#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extracts coordinates from CP2K input format into XYZ files.

Supports scaled CP2K coordinates and scaled XYZ output.

Command Line Interface
----------------------
.. program:: es_cp2k2xyz.py 

.. option:: input

   CP2K input filename.

.. option:: output

   XYZ output filename.

.. option:: --scaled

   Whether the output should be in fractions of the box.

Implementation
--------------
"""

import argparse
import numpy as np
import euston.io as io
import euston.geometry as geo

parser = argparse.ArgumentParser(description='Converts from CP2K input format into XYZ files.')
parser.add_argument('input', type=str, help='Input file name.')
parser.add_argument('output', type=str, help='Output file name.')
parser.add_argument('--scaled', action='store_true', help='Whether the output should be in fractions of the box.')


def main(parser):
    """
    Main routine wrapper.

    :param argparse.ArgumentParser parser: Argument parser
    """

    args = parser.parse_args()
    cp2k = io.Cp2kInput(args.input)

    retval = cp2k.get_cell_vectors()
    if retval is None:
        print 'Unable to find cell information. Aborting.'
        return
    a, b, c = retval

    scaled = cp2k.boolean(cp2k.get_path('FORCE_EVAL / SUBSYS / COORD / SCALED'), False)

    coord_lines = [_ for _ in cp2k.get_path('FORCE_EVAL / SUBSYS / COORD / *') if
                   not _.startswith('SCALED ') and not _.startswith('UNIT ')]
    coordinates = []
    names = []
    for line in coord_lines:
        parts = line.split()
        if len(parts) < 4:
            print 'Ignoring input line: %s' % line
            continue
        names.append(parts[0])
        try:
            coordinates.append(map(float, parts[1:4]))
        except:
            print 'Ignoring input line: %s' % line
            names = names[:-1]
    coordinates = np.array(coordinates)

    if scaled != args.scaled:
        h = np.zeros((3, 3))
        h[0, :] = a
        h[1, :] = b
        h[2, :] = c
        if args.scaled == True:
            coordinates = geo.cartesian_to_scaled_coordinates(coordinates, h)
        else:
            coordinates = geo.scaled_to_cartesian_coordinates(coordinates, h)

    # XYZ output
    lines = []
    lines.append('%d' % len(coordinates))
    lines.append('Unit cell vectors: %f %f %f, %f %f %f, %f %f %f' % tuple(list(a) + list(b) + list(c)))
    for i in range(len(coordinates)):
        lines.append('%s %f %f %f' % (names[i], coordinates[i][0], coordinates[i][1], coordinates[i][2]))
    io.write_lines(args.output, lines)

    print 'Success.'

if __name__ == '__main__':
    main(parser)