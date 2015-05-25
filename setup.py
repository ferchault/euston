#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='euston',
      version='0.1a1',
      description='FileIO and geometry functions for analysis of quantum chemical simulations mostly to work with CP2K',
      author='Guido Falk von Rudorff',
      author_email='guido@vonrudorff.de',
      url='https://github.com/ferchault/euston',
      packages=['euston', ],
      test_suite="tests",
      license='LGPL',
      classifiers=['Development Status :: 3 - Alpha',]
     )
