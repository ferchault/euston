import unittest
import numpy as np
import euston.geometry as geo

class TestGeometry(unittest.TestCase):
    def test_anglebetween(self):
        zero = np.array([0, 0, 0])
        a = np.array([1, 0, 0])
        b = np.array([0, 1, 0])
        b2 = np.array([0, 2, 0])
        self.assertRaises(ValueError, geo._angle_between, zero, a)
        self.assertEqual(90, np.rad2deg(geo._angle_between(a, b)))
        self.assertEqual(90, np.rad2deg(geo._angle_between(a, b2)))
        self.assertEqual(0, np.rad2deg(geo._angle_between(a, a)))
        self.assertEqual(180, np.rad2deg(geo._angle_between(a, -a)))

    def test_hmatabc(self):
        """Tests whether alternating conversions are stable."""
        def _docheck(self, test):
            abc = np.array(test)
            result = geo.hmatrix_to_abc(geo.abc_to_hmatrix(*abc, degrees=True), degrees=True)
            self.assertTrue(np.allclose(np.array(test), result))

        # cubic: Sodium Chloride
        _docheck(self, [5.6402, 5.6402, 5.6402, 90, 90, 90])
        # tetragonal: Rutile
        _docheck(self, [4.5937, 4.5937, 2.9587, 90, 90, 90])
        # orthorhombic: alpha sulphur
        _docheck(self, [10.4646, 12.8660, 24.4860, 90, 90, 90])
        # hexagonal: Graphite
        _docheck(self, [2.461, 2.461, 6.708, 90, 90, 120])
        # trigonal
        _docheck(self, [4.5, 4.5, 4.5, 88, 88, 88])
        # monoclinic: Gypsum
        _docheck(self, [5.679, 15.202, 6.522, 90, 90, 118.43])
        # triclinic: Chalcanthite
        _docheck(self, [6.11, 10.673, 5.95, 97.58, 107.17, 77.55])

