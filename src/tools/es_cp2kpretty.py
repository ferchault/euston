#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reads a CP2K input file and reformats it such that

* Every section is properly indented.
* Every section is terminated with a &END line that includes the name of the section.
* Sections are enclosed in blank lines.

Command Line Interface
----------------------
.. program:: es_fitting.py 

.. option:: input

   CP2K input filename.

.. option:: output

   CP2K output filename.

Implementation
--------------
"""

import argparse
import euston.io as io

parser = argparse.ArgumentParser(description='Pretty-prints CP2K input files.')
parser.add_argument('input', type=str, help='CP2K input filename')
parser.add_argument('output', type=str, help='CP2K output filename')


def main(parser):
	"""
	Main routine wrapper.

	:param argparse.ArgumentParser parser: Argument parser
	"""
	args = parser.parse_args()

	cp2k = io.Cp2kInput(args.input)
	lines = cp2k.to_string(args.output)
	io.write_lines(args.output, lines)


if __name__ == '__main__':
	main(parser)
