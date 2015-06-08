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

    #def test_hmatabc(self):
    #    """Tests whether alternating conversions are stable."""
    #    def _docheck(self, test):
    #        abc = np.array(test)
    #        result = geo.hmatrix_to_abc(geo.abc_to_hmatrix(*abc, degrees=True), degrees=True)
    #        print np.array(test)-result
    #        self.assertTrue(np.allclose(np.array(test), result))

    #    _docheck(self, [1.1, 2.2, 3.3, 90,90,120])
    #    _docheck(self, [1.1, 2.2, 3.3, 90,90,90])
        # triclinic: Chalcanthite
    #    _docheck(self, [6.12, 10.72, 5.96, 82.4,107.3,102.6])
