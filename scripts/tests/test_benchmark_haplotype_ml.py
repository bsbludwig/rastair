import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import benchmark_haplotype_ml as benchmark


class PilotIntervalsTests(unittest.TestCase):
    def test_evaluation_excludes_uncalled_and_unconfident_bases(self):
        confident = [(10, 80), (90, 140), (200, 300)]
        targets = [(50, 110), (120, 130), (150, 180)]
        self.assertEqual(benchmark.intersect_intervals(confident, targets),
                         [(50, 80), (90, 110), (120, 130)])
        self.assertEqual(benchmark.intersect_intervals(confident, []), [])
        self.assertEqual(benchmark.intersect_intervals([(10, 20)], [(20, 30)]), [])

    def test_target_overlaps_merge_before_evaluation(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = pathlib.Path(directory) / 'targets.bed'
            bed.write_text('chr12\t30\t50\nchr19\t1\t9\nchr12\t10\t35\n')
            self.assertEqual(benchmark.intervals(bed, 'chr12'), [(10, 50)])
            with self.assertRaises(ValueError):
                benchmark.intervals(bed, 'chr12', strict=True)
            bed.write_text('chr12\t30\t20\n')
            with self.assertRaises(ValueError):
                benchmark.intervals(bed, 'chr12')


if __name__ == '__main__':
    unittest.main()
