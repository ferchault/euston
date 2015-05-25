from unittest import TestCase

import StringIO
import numpy as np
from euston.io import CubeFile

simple = '''HEADER 1
HEADER 2
1 0 0 0
-1 1 0 0
-1 0 2 0
-1 0 0 3
12 0 1 2 3
1
'''
simple2 = '''HEADER 1
HEADER 2
2 0 0 0
-1 1 0 0
-1 0 2 0
-1 0 0 3
12 0 1 2 3
12 0 2 3 4
1
'''
simple3 = '''HEADER 1
HEADER 2
2 0 0 0
-1 1 0 0
-1 0 2 0
-1 0 0 3
12 0 1 2 3
12 0 2 3 4
1 1 1 1 1 1 1
'''
simple4 = '''HEADER 1
HEADER 2
1 0 0 0
1 1 0 0
-1 0 1 0
1 0 0 1
12 0 1 2 3
1
'''

adv1 = '''-Quickstep-
 HARTREE POTENTIAL
  1    0.000000    0.000000    0.000000
  10    0.154626    0.000154    0.004382
  10   -0.077234    0.134457   -0.001403
  10    0.004438    0.000880    0.155833
  1 0 0 0 0
  '''
adv1 += '1 ' * 1000

class TestCubeFile(TestCase):
        def test_parse_error(self):
            fh = StringIO.StringIO(simple3)
            with self.assertRaises(IndexError):
                CubeFile(filehandle=fh)

        def test_atom_count(self):
            fh = StringIO.StringIO(simple)
            cube = CubeFile(filehandle=fh)
            self.assertEqual(cube.count_atoms(), 1)

            fh = StringIO.StringIO(simple2)
            cube = CubeFile(filehandle=fh)
            self.assertEqual(cube.count_atoms(), 2)

        def test_atom_coordinates(self):
            fh = StringIO.StringIO(simple)
            cube = CubeFile(filehandle=fh)
            coord = cube.get_coordinates()
            ref = np.array([1., 2., 3.])
            self.assertTrue((coord == ref).all())

        def test_voxel_volume(self):
            fh = StringIO.StringIO(simple)
            cube = CubeFile(filehandle=fh)
            self.assertEqual(cube.get_voxel_volume(), 6)

        def test_get_val(self):
            fh = StringIO.StringIO(simple)
            cube = CubeFile(filehandle=fh)
            self.assertEqual(cube.get_val(0, 0, 0), 1)
            with self.assertRaises(IndexError):
                cube.get_val(1, 0, 0)

        def test_get_voxel_pos(self):
            fh = StringIO.StringIO(simple)
            cube = CubeFile(filehandle=fh)
            ref = np.array([0.5, 1., 1.5])
            pos = cube.get_voxel_pos(0, 0, 0, centered=True)
            self.assertTrue((pos==ref).all())
            ref = np.array([0.0, 0.0, 0.0])
            pos = cube.get_voxel_pos(0, 0, 0, centered=False)
            self.assertTrue((pos==ref).all())

        def test_units(self):
            fh = StringIO.StringIO(simple4)
            cube = CubeFile(filehandle=fh)
            self.assertEqual(cube.get_xlen(), 1)
            self.assertEqual(cube.get_ylen(), 1)
            self.assertEqual(cube.get_zlen(), 1)
            h_mat = cube.get_h_matrix()

            self.assertAlmostEqual(h_mat[0, 0], 0.529, 3)
            self.assertEqual(h_mat[1, 1], 1)
            self.assertAlmostEqual(h_mat[2, 2], 0.529, 3)
            h_clean = h_mat-np.diag(np.diag(h_mat))
            off = np.sum(np.abs(h_clean))
            self.assertEqual(off, 0)

            self.assertAlmostEqual(cube.get_voxel_volume(), 0.280, 3)

            fh = StringIO.StringIO(adv1)
            cube = CubeFile(filehandle=fh)
            h_mat = cube.get_h_matrix()
            ref = np.array([[  8.18245553e-01,  -4.08704726e-01,   2.34848846e-02],
                [  8.14932903e-04,   7.11515801e-01,   4.65675945e-03],
                [  2.31885453e-02,  -7.42435626e-03,   8.24632722e-01]])
            for val in np.ravel(np.abs(h_mat-ref)):
                self.assertAlmostEqual(val, 0)

        def test_stringrepresentation(self):
            fh = StringIO.StringIO(simple4)
            cube = CubeFile(filehandle=fh)
            lines = '\n'.join(cube.to_string())
            fh2 = StringIO.StringIO(lines)
            print lines
            cube2 = CubeFile(filehandle=fh2)
            self.assertTrue(np.allclose(cube.get_coordinates(), cube2.get_coordinates()))
            self.assertTrue(np.allclose(cube._data, cube2._data))
            self.assertAlmostEqual(cube.get_voxel_volume(), cube2.get_voxel_volume(),delta=10e-6)
            self.assertEqual(cube.get_xlen(), cube2.get_xlen())
            self.assertEqual(cube.get_ylen(), cube2.get_ylen())
            self.assertEqual(cube.get_zlen(), cube2.get_zlen())
