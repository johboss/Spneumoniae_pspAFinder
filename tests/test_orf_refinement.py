import unittest

from orf_refinement import _map_region_coordinates


class MapRegionCoordinatesTests(unittest.TestCase):
    def test_maps_reverse_strand_orf_without_mirroring(self):
        self.assertEqual(
            _map_region_coordinates(window_start=1_631_525, local_start=1_035, local_end=3_338),
            (1_632_559, 1_634_862),
        )

    def test_maps_forward_strand_orf(self):
        self.assertEqual(
            _map_region_coordinates(window_start=39_000, local_start=740, local_end=2_647),
            (39_739, 41_646),
        )


if __name__ == "__main__":
    unittest.main()
