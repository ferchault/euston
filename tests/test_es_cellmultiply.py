import unittest
import shlex
import tempfile
import os

from tools import es_cellmultiply


class Test_Cellmultiply(unittest.TestCase):
	def _build_args(self, command):
		return es_cellmultiply.parser.parse_args(shlex.split(command))

	def _get_data_file(self, filename):
		return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', filename)

	def _get_tempfilename(self):
		fh = tempfile.NamedTemporaryFile(delete=False)
		fn = fh.name
		fh.close()
		return fn

	def assertSystemExit(self, args, code):
		try:
			es_cellmultiply.main(args)
			self.fail('No SystemExit risen.')
		except SystemExit as e:
			self.assertEqual(code, e.code)

	def test_simple(self):
		"""Tests for required filenames."""
		self.assertRaises(SystemExit, self._build_args, '')

	def test_loadunscaled(self):
		"""Tests whether loading and writing the same file gives identical coordinates."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('--abc 1,1,1,90,90,90 %s %s' % (infile, outfile))
		es_cellmultiply.main(args)
		fh = open(outfile)
		lines = fh.readlines()
		self.assertEqual('2', lines[0].strip())
		self.assertTrue(lines[2].startswith('C '))
		self.assertTrue(lines[3].startswith('H '))
		self.assertEqual([1., 2., 3.], map(float, lines[2].split()[1:]))
		self.assertEqual([4., 5., 6.], map(float, lines[3].split()[1:]))

	def test_loadunscaledtoscaled(self):
		"""Tests the transformation of unscaled to scaled XYZ."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('--abc 1,2,3,90,90,90 %s %s --sc_out' % (infile, outfile))
		es_cellmultiply.main(args)
		fh = open(outfile)
		lines = fh.readlines()
		self.assertEqual('2', lines[0].strip())
		self.assertTrue(lines[2].startswith('C '))
		self.assertTrue(lines[3].startswith('H '))
		self.assertEqual([1., 1., 1.], map(float, lines[2].split()[1:]))
		self.assertEqual([4., 2.5, 2.], map(float, lines[3].split()[1:]))

	def test_hmaterror(self):
		"""Tests for requiring the cell information only for not-scaled setups."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('%s %s' % (infile, outfile))
		self.assertSystemExit(args, 1)
		args = self._build_args('%s %s --sc_in' % (infile, outfile))
		self.assertSystemExit(args, 1)
		args = self._build_args('%s %s --sc_out' % (infile, outfile))
		self.assertSystemExit(args, 1)
		args = self._build_args('%s %s' % (infile, outfile))
		self.assertSystemExit(args, 1)
		args = self._build_args('%s %s --sc_in --sc_out' % (infile, outfile))
		es_cellmultiply.main(args)

	def test_invalidhmat(self):
		"""Tests whether invalid hmatrix entries cancel the program."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('%s %s --sc_in --hmat a,b,c,d,e,f,g,h,i' % (infile, outfile))
		self.assertSystemExit(args, 2)

	def test_invalidabc(self):
		"""Tests whether invalid cell lengths and angles cancel the program."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('%s %s --sc_in --abc a,b,c,d,e,f' % (infile, outfile))
		self.assertSystemExit(args, 3)

		args = self._build_args('%s %s --sc_in --abc 1,2,3,4,5' % (infile, outfile))
		self.assertSystemExit(args, 5)

	def test_invalidrepeat(self):
		"""Tests whether invalid repeat settings cancel the program."""
		infile = self._get_data_file('unscaled.xyz')
		outfile = self._get_tempfilename()

		args = self._build_args('%s %s --sc_in --sc_out --X -1' % (infile, outfile))
		self.assertSystemExit(args, 4)
		args = self._build_args('%s %s --sc_in --sc_out --Y -1' % (infile, outfile))
		self.assertSystemExit(args, 4)
		args = self._build_args('%s %s --sc_in --sc_out --Z -1' % (infile, outfile))
		self.assertSystemExit(args, 4)
