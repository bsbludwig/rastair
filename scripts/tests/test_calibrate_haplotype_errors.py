import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import calibrate_haplotype_errors as calibration


class ErrorCalibrationTests(unittest.TestCase):
    def test_long_tract_contractions_change_the_deletion_model(self):
        counts = {(1, 12): [12000, 10, 300, 400, 1000]}
        profile = calibration.fit_profile(counts, 'illumina', 'taps')
        row = next(row for row in profile['repeats'] if row['period'] == 1 and row['minimum_copies'] == 12)
        self.assertLess(row['deletion_open_phred'], row['insertion_open_phred'])
        self.assertLess(row['deletion_open_phred'], 20)
        self.assertEqual((profile['platform'], profile['library']), ('illumina', 'taps'))

    def test_stutter_steps_are_per_molecule_and_shrunk_toward_the_conservative_prior(self):
        counts = {(1, 14): [42000, 269, 596, 1087, 3000, 60, 536, 240, 29]}
        profile = calibration.fit_profile(counts, 'illumina', 'taps')
        row = next(row for row in profile['stutter'] if row['period'] == 1 and row['minimum_copies'] == 14)
        self.assertAlmostEqual(row['steps'][1], (536 + 100 * 14 * 10 ** -3.7 * 0.5) / 3100, places=9)
        self.assertLess(row['steps'][3], row['steps'][1])
        sparse = next(row for row in profile['stutter'] if row['period'] == 6 and row['minimum_copies'] == 20)
        self.assertLess(sum(sparse['steps']), 0.5)
        self.assertGreater(sum(sparse['steps']), 0.0)
        old = calibration.fit_profile({(1, 12): [12000, 10, 300, 400, 1000]}, 'illumina', 'taps')
        self.assertEqual(len(old['stutter']), 6 * 19)

    def test_conversion_is_not_an_indel_or_arbitrary_substitution(self):
        self.assertTrue(calibration.converted_match('ATGT', 'ACGC', True))
        self.assertFalse(calibration.converted_match('ATGT', 'ACGC', False))
        self.assertFalse(calibration.converted_match('AGGT', 'ACGC', True))

    def test_variant_exclusion_handles_half_open_intervals(self):
        intervals = calibration.merge([(10, 20), (18, 25), (30, 40)])
        starts = [start for start, _ in intervals]
        self.assertFalse(calibration.overlaps(intervals, starts, 25, 30))
        self.assertTrue(calibration.overlaps(intervals, starts, 24, 30))


if __name__ == '__main__':
    unittest.main()
