import unittest

import euston.io as io

class Test_DCD(unittest.TestCase):
    def _get_data_file(self, filename):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', filename)

    def test_create(self):
	self.assertTrue(io.HAS_MDA)
	self.assertTrue(io.HAS_MDT)
