import unittest

import StringIO
import numpy as np
from euston.io import XYZ
simple1 = '''2
COMMENT
C 1 2 3
H 4 5 6'''

class TestXYZ(unittest.TestCase):
    def test_create(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)
        ref = np.linspace(1, 6, 6).reshape((2, 3))
        self.assertTrue(np.all(xyz.get_coordinates() == ref))

        copy = simple1.split('\n')
        copy[-1] += 'a'
        fh = StringIO.StringIO('\n'.join(copy))
        self.assertRaises(ValueError, XYZ, filehandle=fh)

        copy = simple1.split('\n')
        copy[-1] += ' a'
        fh = StringIO.StringIO('\n'.join(copy))
        self.assertRaises(ValueError, XYZ, filehandle=fh)

        copy = simple1.split('\n')
        copy[0] = 'a'
        fh = StringIO.StringIO('\n'.join(copy))
        self.assertRaises(ValueError, XYZ, filehandle=fh)

    def test_tostring_preserve_comment(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)
        self.assertEqual('COMMENT', xyz.to_string()[1])

    def test_tostring_preserve_coordinateslabels(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)
        ref = np.linspace(1, 6, 6).reshape((2, 3))
        ls = xyz.to_string()[-2:]
        labels = [_[0] for _ in simple1.split('\n')[-2:]]
        for i in range(2):
            p = ls[i].split()
            self.assertEqual(p[0], labels[i])
            vals = map(float, p[1:])
            self.assertTrue(np.all(vals == ref[i]))

    def test_invalidatomcount(self):
        copy = simple1.split('\n')
        copy[0] = '10'
        fh = StringIO.StringIO('\n'.join(copy))
        self.assertRaises(ValueError, XYZ, filehandle=fh)

    def test_nomultiframe(self):
        fh = StringIO.StringIO(simple1 + '\n' + simple1)
        self.assertRaises(NotImplementedError, XYZ, filehandle=fh)

    def test_countatoms(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)
        self.assertEqual(2, xyz.count_atoms())

    def test_getlabels(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)
        labels = [_[0] for _ in simple1.split('\n')[-2:]]
        self.assertEqual(labels, xyz.get_labels())

    def test_setdata(self):
        fh = StringIO.StringIO(simple1)
        xyz = XYZ(filehandle=fh)

        ref = np.linspace(1, 6, 6).reshape((2, 3))
        xyz.set_data(['C','C'], ref)
        self.assertRaises(ValueError, xyz.set_data, ['C',], ref)