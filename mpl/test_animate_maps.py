#!/usr/bin/env python3

import unittest
from animate_maps import rotate_longitude

class TestAnimateMaps(unittest.TestCase):

    def test_rotate_longitude(self):
        samples_per_day = 24
        dlon = -360. / samples_per_day
        self.assertEqual(rotate_longitude(0, samples_per_day), 0)
        self.assertEqual(rotate_longitude(samples_per_day, samples_per_day), 0) 
        self.assertEqual(rotate_longitude(1, samples_per_day), dlon)
        self.assertEqual(rotate_longitude(samples_per_day + 1, samples_per_day), dlon)

if __name__ == '__main__':
    unittest.main()
