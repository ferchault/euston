#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the RDF / g(r) for a periodic monoatomic lattice.

.. autofunction:: main

Command Line Interface
----------------------

.. program:: es_latticerdf.py

.. option:: maxr

   The maximum distance for the RDF.

.. option:: lattice

   The lattice type. Supported: fcc

.. option:: abc

   a, b, c, alpha, beta, gamma. Comma-separated without spaces.

.. option:: --radians

   Whether angles in --abc are given in radians.
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

parser = argparse.ArgumentParser(
	description='Calculates the RDF / g(r) for a periodic monoatomic lattice.')
parser.add_argument('maxr', type=float, help='The maximum distance for the RDF.')
parser.add_argument('lattice', type=str, help='Lattice type.', choices='fcc bcc cubic'.split())
parser.add_argument('abc', type=str, help='a, b, c, alpha, beta, gamma. Comma-separated without spaces.')
parser.add_argument('--radians', action='store_true', help='Whether angles in --abc are given in radians.')


def main(parser):
	"""
	Main routine wrapper.

	:param argparse.ArgumentParser parser: Argument parser
	"""
	args = parser.parse_args()

	# get h mat
	try:
		abc = map(float, args.abc.split(','))
	except:
		print 'Invalid abc entries.'
		exit(3)
	if len(abc) != 6:
		print 'Not enough entries for cell lengths.'
		exit(5)
	h_matrix = geom.abc_to_hmatrix(*abc, degrees=(not args.radians))

	# count repetitions
	maxr = (args.maxr + geom.cell_longest_diameter(h_matrix)) * 2
	max_a = geom.vector_repetitions(maxr, h_matrix, 0)
	max_b = geom.vector_repetitions(maxr, h_matrix, 1)
	max_c = geom.vector_repetitions(maxr, h_matrix, 2)

	# build primitive unit cell
	if args.lattice == 'fcc':
		pos = np.array(((0.,0.,0.), (0.,0.,1.), (0.,1.,1.), (0.,1.,0.), (0.,0.5,0.5),
			(1.,0.,0.), (1.,0.,1.), (1.,1.,1.), (1.,1.,0.), (1.,0.5,0.5),
			(0.5,0.5,0.), (0.5,0.,0.5), (0.5, 1.,0.5), (0.5, 0.5, 1.)))
	elif args.lattice == 'bcc':
		pos = np.array(((0.,0.,0.), (0.,0.,1.), (0.,1.,1.), (0.,1.,0.),
			(1.,0.,0.), (1.,0.,1.), (1.,1.,1.), (1.,1.,0.), (0.5,0.5,0.5)))
	elif args.lattice == 'cubic':
		pos = np.array(((0.,0.,0.), (0.,0.,1.), (0.,1.,1.), (0.,1.,0.),
			(1.,0.,0.), (1.,0.,1.), (1.,1.,1.), (1.,1.,0.)))

	# repeat
	multiplied = geom.cell_multiply(pos, max_a, max_b, max_c, h_matrix=h_matrix, scaling_in=True)
	reference = geom.cell_multiply(pos, 1, 1, 1, h_matrix=h_matrix, scaling_in=True)
	reference += geom.repeat_vector(h_matrix, max_a / 2, max_b / 2, max_c / 2)

	# calculating distances
	print 'Calculating all %d distances' % (len(reference) * len(multiplied))
	distances = []
	for i in range(len(reference)):
		for j in range(len(multiplied)):
			# workaround: potential bug in cell_multiply giving identical coordinates
			dist = geom.distance_pbc(reference[i], multiplied[j], h_matrix)
			if dist < 0.01:
				continue
			distances.append(dist)

	# binning
	hist, bins = np.histogram(distances, bins=1000, range=(0, args.maxr))
	bins = (bins[1:]+bins[:-1])/2
	for b, h in zip(bins, hist):
		print b, h

if __name__ == '__main__':
	main(parser)
