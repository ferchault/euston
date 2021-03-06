import unittest
import numpy as np
import os

from euston.io import Cp2kInput
import StringIO

simple1 = '''&FORCE_EVAL
&END DFT'''

simple2 = '''&SECTIONA
&SECTIONB
FOOBAR 1
&END
SNAFU 2
&END SECTIONA'''

simple3 = '''&SECTIONA
# comment1
! comment2

&END SECTIONA'''

simple4 = '''&DFT


&END DFT'''


class TestCp2kInput(unittest.TestCase):
	def _get_data_file(self, filename):
		return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', filename)

	def test_load(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14.inp'))

	def test_boolean(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14.inp'))
		self.assertEqual(True, cp2k.boolean(None, True))
		self.assertEqual(False, cp2k.boolean(None, False))
		self.assertEqual(True, cp2k.boolean('.T.', False))
		self.assertEqual(True, cp2k.boolean('TRUE', False))
		self.assertEqual(True, cp2k.boolean('T', False))
		self.assertEqual(False, cp2k.boolean('F', True))

	def test_getpath(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14.inp'))
		self.assertEqual('1.0E-10', cp2k.get_path('FORCE_EVAL / DFT / QS / EPS_DEFAULT'))
		fh = StringIO.StringIO(simple1)
		cp2k = Cp2kInput(filehandle=fh)
		self.assertRaises(ValueError, cp2k.get_path, '')

	def test_getpath_multiple(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14-2.inp'))
		vals = cp2k.get_path('FORCE_EVAL / SUBSYS / COORD / *')
		self.assertEqual(['Ar 1 0 0', 'Ar 0 0 0'], vals)

	def test_getkeywordchecked(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14-2.inp'))
		self.assertEqual(None, cp2k.get_keyword_checked('FORCE_EVAL / SUBSYS / DUMMY', float))
		self.assertEqual(['a'], cp2k.get_keyword_checked('FORCE_EVAL / SUBSYS / DUMMY', str))
		self.assertEqual(['a'], cp2k.get_keyword_checked('FORCE_EVAL / SUBSYS / DUMMY'))

		fh = StringIO.StringIO(simple1)
		cp2k = Cp2kInput(filehandle=fh)
		self.assertEqual(None, cp2k.get_keyword_checked(''))

	def test_cellvectors(self):
		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14.inp'))
		a, b, c = cp2k.get_cell_vectors()
		self.assertTrue(np.allclose(a, np.array([4, 0, 0])))
		self.assertTrue(np.allclose(b, np.array([0, 5, 0])))
		self.assertTrue(np.allclose(c, np.array([0, 0, 6])))

		cp2k = Cp2kInput(self._get_data_file('regtest-ot-1-Ar-14.inp'))
		a, b, c = cp2k.get_cell_vectors()
		self.assertTrue(np.allclose(a, np.array([4, 0, 0])))
		self.assertTrue(np.allclose(b, np.array([0, 5, 0])))
		self.assertTrue(np.allclose(c, np.array([0, 0, 6])))

	def test_tostring(self):
		cp2k = Cp2kInput()
		cp2k.set_lines(simple2.split('\n'))
		expected = ['&SECTIONA', '  &SECTIONB', '    FOOBAR 1', '  &END SECTIONB', '', '  SNAFU 2', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string())
		expected = ['&SECTIONA', '  &SECTIONB', '    FOOBAR 1', '  &END', '', '  SNAFU 2', '&END']
		self.assertEqual(expected, cp2k.to_string(close_sections=False))
		cp2k.set_lines(simple1.split('\n'))
		self.assertRaises(ValueError, cp2k.to_string)
		cp2k.set_lines(simple3.split('\n'))
		expected = ['&SECTIONA', '  # comment1', '  ! comment2', '', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string(keep_empty=True))
		cp2k.set_lines(simple4.split('\n'))
		expected = ['&DFT', '&END DFT']
		self.assertEqual(expected, cp2k.to_string())

	def test_tostringcomments(self):
		cp2k = Cp2kInput()
		cp2k.set_lines(simple3.split('\n'))
		expected = ['&SECTIONA', '  # comment1', '  ! comment2', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string())
		expected = ['&SECTIONA', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string(keep_comments=False))
		expected = ['&SECTIONA', '# comment1', '! comment2', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string(indent_comments=False))
		expected = ['&SECTIONA', '&END SECTIONA']
		self.assertEqual(expected, cp2k.to_string(keep_comments=False, indent_comments=True))

	def test_setlines(self):
		cp2k = Cp2kInput()
		cp2k.set_lines(simple1.split('\n'))
		self.assertEqual(simple1, '\n'.join(cp2k._lines))