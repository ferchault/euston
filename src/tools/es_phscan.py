#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates possible pH values for a given water volume.

Implementation
--------------
"""

import argparse
import numpy as np
import euston.io as io
import euston.geometry as geo
import math

parser = argparse.ArgumentParser(description='Calculates possible pH values at 25°C.')
parser.add_argument('--hmat', type=str,
					help='H matrix (cell vectors in columns). Comma-separated in row-first notation without spaces.')
parser.add_argument('--abc', type=str, help='a, b, c, alpha, beta, gamma. Comma-separated without spaces.')
parser.add_argument('--radians', action='store_true', help='Whether angles in --abc are given in radians.')


def main(parser):
	"""
	Main routine wrapper.

	:param argparse.ArgumentParser parser: Argument parser
	"""
	args = parser.parse_args()

	hmat = None
	if (args.hmat is None and args.abc is None) or (args.hmat is not None and args.abc is not None):
		print 'Please specify either the H matrix or lattice constants unless input and output are scaled.'
		exit(1)

	if args.hmat is not None:
		try:
			hmat = map(float, args.hmat.split(','))
		except:
			print 'Invalid H matrix entries.'
			exit(2)
		hmat = np.array(hmat).reshape((3, 3)).T

	if args.abc is not None:
		try:
			abc = map(float, args.abc.split(','))
		except:
			print 'Invalid abc entries.'
			exit(2)
		hmat = geo.abc_to_hmatrix(*abc, degrees=(not args.radians))

	volume = geo.cell_volume(hmat)
	base_count = 6.022 * 10 ** 23 * volume * 10 ** (-30) * 1000
	print 'Volume:', volume
	print 'Expecting this many water molecules:', base_count * 55.4
	print 'Expecting this many H+ ions:', base_count * 10e-7
	# ph = -lg(c(H+))
	print -math.log(3 / base_count)


if __name__ == '__main__':
	main(parser)
