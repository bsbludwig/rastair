import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import benchmark_haplotype_likelihoods as benchmark


class BenchmarkTests(unittest.TestCase):
    def test_configured_arms_keep_distinct_search_limits_on_common_library_and_engine(self):
        manifest = {"arms": [
            {"name": "small", "library": "taps", "engine": "aligned", "arguments": ["--max-sites", "3"]},
            {"name": "shrinkage", "library": "taps", "engine": "aligned", "arguments": ["--max-haplotypes", "64"]},
        ]}
        self.assertEqual(benchmark.arms(manifest), [("small", "taps", "aligned"), ("shrinkage", "taps", "aligned")])
        self.assertEqual(benchmark.arm_arguments(manifest,"shrinkage"), ["--max-haplotypes", "64"])
        manifest["arms"][1]["name"] = "small"
        with self.assertRaises(ValueError):
            benchmark.arms(manifest)

    def test_panel_weights_reconstruct_the_candidate_population(self):
        candidates = {("chr1", p) for p in range(1, 101)}
        truth = {("chr1", 1), ("chr1", 2)}
        mandatory = {("chr1", 3), ("chr1", 4)}
        panel, sampling = benchmark.select_panel(candidates, truth, mandatory, 12, 7)
        self.assertEqual(len(panel), 16)
        self.assertAlmostEqual(sum(m["weight"] for m in panel.values()), 100)
        self.assertEqual(sampling["negative_pool"], 96)
        reordered, _ = benchmark.select_panel(set(reversed(sorted(candidates))), truth, mandatory, 12, 7)
        self.assertEqual(panel, reordered)
        self.assertEqual(panel[("chr1", 3)]["weight"], 1)

    def test_weighted_calibration_and_panel_counts_have_different_denominators(self):
        data = [{"posterior": ".9", "label": True, "weight": 1},
                {"posterior": ".8", "label": False, "weight": 10},
                {"posterior": ".01", "label": False, "weight": 10}]
        score = benchmark.allele_metrics(data)
        self.assertEqual((score["tp"], score["fp"], score["fn"]), (1, 1, 0))
        self.assertEqual(score["estimated_candidate_population_fp"], 10)
        self.assertEqual(score["weighted_alleles"], 21)
        self.assertAlmostEqual(score["weighted_brier"], (.1 ** 2 + 10 * .8 ** 2 + 10 * .01 ** 2) / 21)

    def test_genotypes_are_scored_as_unphased_base_pairs(self):
        lines = ["chr1\t10\t.\tA\tG\t.\tPASS\t.\tDP:GT\t8:1|0",
                 "chr1\t11\t.\tC\tA,G\t.\tPASS\t.\tGT\t2/1",
                 "chr1\t12\t.\tC\tCA\t.\tPASS\t.\tGT\t0/1",
                 "chr1\t13\t.\tCA\tC\t.\tPASS\t.\tGT\t1/1"]
        genotypes = benchmark.truth_genotypes(lines)
        self.assertEqual(genotypes, {("chr1", 10): ("A", "G"), ("chr1", 11): ("A", "G")})
        probabilities = [0.0] * 10
        probabilities[benchmark.GENOTYPES.index(("A", "G"))] = 1.0
        row = {"chrom": "chr1", "pos": "10", "ref": "A", "alt": "G",
               "truth_genotype": ("A", "G"), "genotype_posteriors": ",".join(map(str, probabilities))}
        score = benchmark.genotype_metrics([row, dict(row)])
        self.assertEqual(score["truth_snp_sites"], 1)
        self.assertEqual(score["genotype_errors"], 0)
        self.assertEqual(score["mean_log_loss"], 0)

    def test_probability_validation_rejects_missing_targets_and_bad_genotypes(self):
        rows = [{"chrom": "chr1", "pos": "10", "ref": "A", "alt": alt,
                 "posterior": ".4", "p_het": ".3", "p_hom_alt": ".1", "expected_fraction": ".25",
                 "site_posterior": ".4", "offered_linking_fragments": "2", "retained_linking_fragments": "1",
                 "genotype_posteriors": ",".join([".1"] * 10),
                 "site_genotype_posteriors": ",".join([".1"] * 10)} for alt in "CGT"]
        self.assertEqual(benchmark.validate(rows, {("chr1", 10)}), 3)
        for row in rows:
            row.update(observations="1", read_observations="2", overlapping_fragments="1", discordant_fragments="1", fragments="1")
        self.assertEqual(benchmark.validate(rows, {("chr1", 10)}), 3)
        rows[0]["observations"] = "2"
        with self.assertRaises(ValueError):
            benchmark.validate(rows, {("chr1", 10)})
        rows[0]["observations"] = "1"
        with self.assertRaises(ValueError):
            benchmark.validate(rows, {("chr1", 10), ("chr1", 11)})
        rows[0]["genotype_posteriors"] = ",".join([".2"] * 10)
        with self.assertRaises(ValueError):
            benchmark.validate(rows, {("chr1", 10)})


if __name__ == "__main__":
    unittest.main()
