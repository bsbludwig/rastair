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

    def test_profiles_carry_no_stutter_table(self):
        profile = calibration.fit_profile({(1, 12): [12000, 10, 300, 400, 1000]}, 'illumina', 'taps')
        self.assertEqual(set(profile), {'version', 'platform', 'library', 'maximum_evidence_phred', 'repeats'})

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
