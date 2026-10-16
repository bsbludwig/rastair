import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import benchmark_phasing as benchmark


class PhasingBenchmarkTests(unittest.TestCase):
    def inspect(self, contents):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / 'calls.vcf'
            path.write_text(contents)
            return benchmark.inspect_vcf(path)

    def test_hash_ignores_phase_orientation_but_detects_genotype_changes(self):
        plain = 'chr1\t10\t.\tA\tG\t40\tPASS\tDP=12\tGT:DP\t0/1:12\nchr1\t20\t.\tT\tC\t40\tPASS\tDP=12\tGT:DP\t0/1:12\n'
        phased = plain.replace('GT:DP\t0/1:12', 'GT:PS:DP\t0|1:10:12', 1).replace('GT:DP\t0/1:12', 'GT:PS:DP\t1|0:10:12')
        a = self.inspect(plain)
        b = self.inspect(phased)
        self.assertEqual(a['nonphase_sha256'], b['nonphase_sha256'])
        self.assertEqual(b['phased_sites'], 2)
        self.assertEqual(b['phase_blocks'], 1)
        changed = plain.replace('0/1:12', '1/1:12', 1)
        self.assertNotEqual(a['nonphase_sha256'], self.inspect(changed)['nonphase_sha256'])

    def test_orphan_phase_sets_and_phased_indels_are_rejected(self):
        orphan = 'chr1\t10\t.\tA\tG\t40\tPASS\tDP=12\tGT:PS\t0|1:10\n'
        with self.assertRaises(ValueError):
            self.inspect(orphan)
        insertion = orphan.replace('\tG\t40', '\tGT\t40')
        with self.assertRaises(ValueError):
            self.inspect(insertion)


if __name__ == '__main__':
    unittest.main()
