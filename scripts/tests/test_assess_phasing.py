import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from assess_phasing import Intervals, assess, compare_shared, phase_errors


class PhasingAssessmentTests(unittest.TestCase):
    def test_block_flips_do_not_count_and_single_site_flips_count_twice(self):
        report = phase_errors({1: [(10, 1), (20, 1)], 2: [(30, 0), (40, 1), (50, 0)], 3: [(60, 1)]})
        self.assertEqual(report['adjacent_pairs'], 3)
        self.assertEqual(report['switch_errors'], 2)
        self.assertEqual(report['minimum_orientation_disagreements'], 1)
        self.assertEqual(report['assessable_sites'], 5)

    def test_masks_genotype_mismatches_and_truth_blocks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            confidence = root / 'confidence.bed'
            confidence.write_text('chr1\t0\t100\n')
            mask = root / 'mask.txt'
            mask.write_text('Chr\tStart\tEnd\tGapType\nchr1\t19\t20\tcentromere\n')
            truth = root / 'truth.vcf'
            query = root / 'query.vcf'
            rows = ['chr1\t{}\t.\tA\tG\t40\tPASS\t.\tGT:PS\t{}:{}\n']
            truth.write_text(''.join(rows[0].format(pos, '0|1', 1 if pos < 40 else 2) for pos in [10, 20, 30, 40, 50]))
            query.write_text(''.join(rows[0].format(pos, gt, 10) for pos, gt in [(10, '0|1'), (20, '1|0'), (30, '1|1'), (40, '1|0'), (50, '1|0')]))
            report = assess(truth, query, 'chr1', 100, Intervals([confidence]), Intervals([mask]))
            self.assertEqual(report['confident_unmasked_truth_hets'], 4)
            self.assertEqual(report['genotype_disagreements_at_truth_hets'], 1)
            self.assertEqual(report['adjacent_pairs'], 1)
            self.assertEqual(report['switch_errors'], 0)

    def test_matched_comparison_uses_identical_links_in_joint_blocks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            bed = root / 'confidence.bed'
            bed.write_text('chr1\t0\t100\n')
            truth, a, b = [root / name for name in ['truth.vcf', 'a.vcf', 'b.vcf']]
            row = 'chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT:PS\t{}:{}\n'
            truth.write_text(''.join(row.format(pos, '0|1', 1) for pos in [10, 20, 30]))
            a.write_text(''.join(row.format(pos, gt, 10) for pos, gt in [(10, '0|1'), (20, '1|0'), (30, '0|1')]))
            b.write_text(''.join(row.format(pos, '0|1', ps) for pos, ps in [(10, 10), (20, 10), (30, 30)]))
            report = compare_shared(truth, [a, b], 'chr1', 100, Intervals([bed]), Intervals([]))
            self.assertEqual(report['shared_genotype_matching_phased_sites'], 3)
            self.assertEqual(report['arms'][str(a)]['adjacent_pairs'], 1)
            self.assertEqual(report['arms'][str(b)]['adjacent_pairs'], 1)
            self.assertEqual(report['arms'][str(a)]['switch_errors'], 1)
            self.assertEqual(report['arms'][str(b)]['switch_errors'], 0)


if __name__ == '__main__':
    unittest.main()
