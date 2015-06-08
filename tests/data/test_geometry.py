import unittest
import numpy as np
import euston.geometry as geo

# cubic: Sodium Chloride
cubic = [5.6402, 5.6402, 5.6402, 90, 90, 90]
# tetragonal: Rutile
tetragonal = [4.5937, 4.5937, 2.9587, 90, 90, 90]
# orthorhombic: alpha sulphur
orthorhombic = [10.4646, 12.8660, 24.4860, 90, 90, 90]
# hexagonal: Graphite
hexagonal = [2.461, 2.461, 6.708, 90, 90, 120]
# trigonal
trigonal = [4.5, 4.5, 4.5, 88, 88, 88]
# monoclinic: Gypsum
monoclinic = [5.679, 15.202, 6.522, 90, 90, 118.43]
# triclinic: Chalcanthite
triclinic = [6.11, 10.673, 5.95, 97.58, 107.17, 77.55]
# collection
bravais_lattices = (cubic, tetragonal, orthorhombic, hexagonal, trigonal, monoclinic, triclinic)

class TestGeometry(unittest.TestCase):
    def test_anglebetween(self):
        zero = np.array([0, 0, 0])
        a = np.array([1, 0, 0])
        b = np.array([0, 1, 0])
        b2 = np.array([0, 2, 0])
        nan = np.array([np.nan, 0, 0])
        self.assertRaises(ValueError, geo._angle_between, zero, a)
        self.assertEqual(90, np.rad2deg(geo._angle_between(a, b)))
        self.assertEqual(90, np.rad2deg(geo._angle_between(a, b2)))
        self.assertEqual(0, np.rad2deg(geo._angle_between(a, a)))
        self.assertRaises(ValueError, geo._angle_between, nan, nan)
        self.assertEqual(180, np.rad2deg(geo._angle_between(a, -a)))

    def test_hmatabc(self):
        """Tests whether alternating conversions are stable."""
        def _docheck(self, test):
            abc = np.array(test)
            result = geo.hmatrix_to_abc(geo.abc_to_hmatrix(*abc, degrees=True), degrees=True)
            self.assertTrue(np.allclose(np.array(test), result))

        for system in bravais_lattices:
            _docheck(self, system)

    def test_boxvertices(self):
        abc = [1, 1, 1, 90, 90, 90]
        hmat = geo.abc_to_hmatrix(*abc, degrees=True)
        refdata = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0], [1, 0, 1], [1, 0, 0]])
        self.assertTrue(np.allclose(geo.box_vertices(hmat, 0, 0, 0), refdata))
        refdata += np.array([0, 0, 1])
        self.assertTrue(np.allclose(geo.box_vertices(hmat, 0, 0, 1), refdata))
        refdata += np.array([0, 2, 0])
        self.assertTrue(np.allclose(geo.box_vertices(hmat, 0, 2, 1), refdata))
        refdata += np.array([3, 0, 0])
        self.assertTrue(np.allclose(geo.box_vertices(hmat, 3, 2, 1), refdata))

    def test_cellvolume(self):
        abc = [1, 1, 1, 90, 90, 90]
        hmat = geo.abc_to_hmatrix(*abc, degrees=True)
        self.assertEqual(1, geo.cell_volume(hmat))

    def test_scaledcartesian(self):
        """Tests whether alternating conversions are stable."""
        for fractions in ([0, 0, 0], [1, 1, 1], [0.1, 0.2, 0.3], [2.1, 2.2, 2.3]):
            base = np.array([fractions])
            for system in bravais_lattices:
                hmat = geo.abc_to_hmatrix(*system)
                self.assertTrue(np.allclose(geo.cartesian_to_scaled_coordinates(geo.scaled_to_cartesian_coordinates(base, hmat), hmat), base))

    def test_cellmultiply(self):
        base = np.array([[0.5, 0.5, 0.5]])
        self.assertRaises(ValueError, geo.cell_multiply, base, -1, -1, -1)
        self.assertRaises(ValueError, geo.cell_multiply, base, 1.1, 1.1, 1.1)

        for sc_in in (True, False):
            for sc_out in (True, False):
                for system in bravais_lattices:
                    hmat = geo.abc_to_hmatrix(*system)
                    result = geo.cell_multiply(base, 1, 1, 2, h_matrix=hmat, scaling_in=sc_in, scaling_out=sc_out)
                    reference = np.copy(base)
                    if sc_in:
                        reference = geo.scaled_to_cartesian_coordinates(reference, hmat)
                    if sc_out:
                        result = geo.scaled_to_cartesian_coordinates(result, hmat)
                    # single result part of multiplied result
                    result2 = geo.cell_multiply(base, 1, 1, 1, h_matrix=hmat, scaling_in=sc_in, scaling_out=False)
                    result2 *= np.array([1, 1, 0.5])
                    residuals = np.linalg.norm(result-result2[0], axis=1)
                    #self.assertTrue(min(residuals) < 10e-5)
