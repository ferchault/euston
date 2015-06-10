#!/usr/bin/env python

from setuptools import setup, find_packages
import sys
import os

sys.path.append(os.path.abspath('..'))

setup(name='euston',
	  version='0.1a1',
	  description='FileIO and geometry functions for analysis of quantum chemical simulations mostly to work with CP2K',
	  author='Guido Falk von Rudorff',
	  author_email='guido@vonrudorff.de',
	  url='https://github.com/ferchault/euston',
	  packages=['euston', ],
	  test_suite="tests",
	  license='LGPL',
	  classifiers=['Development Status :: 3 - Alpha', ],
	  scripts=['tools/es_cellmultiply.py', 'tools/es_cp2k2xyz.py', 'tools/es_cp2kperf.py', 'tools/es_cp2kpretty.py',
			   'tools/es_fitting.py', 'tools/es_phscan.py', 'tools/es_projectcube.py', 'tools/es_wrapcube.py'],
)
