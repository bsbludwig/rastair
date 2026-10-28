pub mod cpg;
pub mod denovo_cpg;
pub mod others;

pub use cpg::CpgFeatures;
pub use denovo_cpg::DenovoCpgFeatures;
pub use others::OthersFeatures;

#[cfg(test)]
mod tests {
    use super::{CpgFeatures, DenovoCpgFeatures};
    use crate::{
        call::tests::{
            pileup::{SimpleRead, SimpleReads},
            utils::{ReadLine, create_pileups},
        },
        metrics::PileupMetrics,
    };
    use hegel::{TestCase, generators as gs};
    use seqair_types::{Base, Base::*, Strand};

    const BASES: [Base; 4] = [A, C, G, T];

    /// Reads over two columns: each read's two bases and its strand.
    type Reads = Vec<([Base; 2], Strand)>;

    fn complement(base: Base) -> Base {
        match base {
            A => T,
            C => G,
            G => C,
            T => A,
            other => other,
        }
    }

    /// The same molecule read from the other strand: the reference and every
    /// read reverse-complemented, OT and OB swapped.
    fn mirror(reference: [Base; 2], reads: &Reads) -> ([Base; 2], Reads) {
        let flip = |[a, b]: [Base; 2]| [complement(b), complement(a)];
        let other = |strand| match strand {
            Strand::OT => Strand::OB,
            Strand::OB => Strand::OT,
        };
        (
            flip(reference),
            reads.iter().map(|&(bases, strand)| (flip(bases), other(strand))).collect(),
        )
    }

    /// The two columns, with this column's reads at `quality` where given.
    fn columns(
        reference: [Base; 2],
        reads: &Reads,
        quality: Option<(usize, u8)>,
    ) -> [PileupMetrics; 2] {
        let lines = reads
            .iter()
            .map(|&(bases, strand)| ReadLine { bases: bases.to_vec(), strand })
            .collect();
        let (segment, mut pileups) = create_pileups(reference.to_vec(), lines);
        if let Some((column, qual)) = quality
            && let Some(pileup) = pileups.get_mut(column)
        {
            pileup.reads = SimpleReads(
                pileup.reads.iter().map(|read| SimpleRead { qual, ..read.clone() }).collect(),
            );
        }
        [0, 1].map(|i| pileups[i].metrics(&segment).expect("metrics"))
    }

    /// `n` reads over two columns, the first carrying `alt` in the first.
    #[hegel::composite]
    fn reads(tc: &TestCase, alt: Base) -> Reads {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(12));
        (0..n)
            .map(|i| {
                let first = if i == 0 { alt } else { tc.draw(gs::sampled_from(BASES.as_slice())) };
                let second = tc.draw(gs::sampled_from(BASES.as_slice()));
                ([first, second], tc.draw(gs::sampled_from([Strand::OT, Strand::OB].as_slice())))
            })
            .collect()
    }

    #[track_caller]
    fn assert_close(name: &str, one: f32, other: f32, reads: &Reads) {
        assert!((one - other).abs() < 1e-5, "{name}: {one} against {other}, {reads:?}");
    }

    /// TAPS reads a methylated CpG as a `T` at its `C` on OT and as an `A` at
    /// its `G` on OB, so the CpG model's features at a `C>T` with the `G` after
    /// are those at the `G>A`, `C` before, of the same reads from the other
    /// strand.
    #[hegel::test]
    fn cpg_features_are_strand_symmetric(tc: TestCase) {
        let reads = tc.draw(reads(T));
        let [c, g] = columns([C, G], &reads, None);
        let ([mc, mg], mirrored) = mirror([C, G], &reads);
        let [before, at] = columns([mc, mg], &mirrored, None);

        let here = CpgFeatures::extract(&c.alt_metrics(T).expect("the alt"), None, Some(&g))
            .expect("features");
        let there = CpgFeatures::extract(&at.alt_metrics(A).expect("the alt"), Some(&before), None)
            .expect("features");
        assert_close("alt_ad_adj", here.alt_ad_adj, there.alt_ad_adj, &reads);
        assert_close("alt_score_adj", here.alt_score_adj, there.alt_score_adj, &reads);
        assert_close("beta_ratio", here.beta_ratio, there.beta_ratio, &reads);
    }

    /// The de-novo model's features at an `X>C` before a reference `G` are
    /// those at the complementary `X>G` after a reference `C`, of the same
    /// reads from the other strand.
    #[hegel::test]
    fn denovo_features_are_strand_symmetric(tc: TestCase) {
        let reference = [tc.draw(gs::sampled_from([A, G, T].as_slice())), G];
        let reads = tc.draw(reads(C));
        let [x, g] = columns(reference, &reads, None);
        let (mirrored_reference, mirrored) = mirror(reference, &reads);
        let [before, at] = columns(mirrored_reference, &mirrored, None);

        let here = DenovoCpgFeatures::extract(&x.alt_metrics(C).expect("the alt"), None, Some(&g))
            .expect("features");
        let there =
            DenovoCpgFeatures::extract(&at.alt_metrics(G).expect("the alt"), Some(&before), None)
                .expect("features");
        assert_close("alt_ad_adj", here.alt_ad_adj, there.alt_ad_adj, &reads);
        assert_close("alt_score_adj", here.alt_score_adj, there.alt_score_adj, &reads);
        assert_close("sb_adj", here.sb_adj, there.sb_adj, &reads);
        assert_close("beta_ratio", here.beta_ratio, there.beta_ratio, &reads);
    }

    /// A de-novo CpG's `beta_center` is the share of `T` among the `C` and `T`
    /// reads on OT at an `X>C` (of `A` among `G` and `A` on OB at an `X>G`),
    /// also where that base is the reference: at a `T>C` the methylated reads
    /// of the new `C` read as the reference `T`. A strand mirror cannot show
    /// this: it turns a `T>C` into an `A>G`, which loses the same reads.
    #[hegel::test]
    fn the_denovo_beta_counts_converted_reads_of_the_reference_base(tc: TestCase) {
        let reference = tc.draw(gs::sampled_from([A, G, T].as_slice()));
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(12));
        let reads: Vec<Base> = (0..n)
            .map(|i| if i == 0 { C } else { tc.draw(gs::sampled_from([C, T].as_slice())) })
            .collect();
        denovo_beta_counts_converted_reads(reference, &reads);
    }

    /// The case [`the_denovo_beta_counts_converted_reads_of_the_reference_base`]
    /// draws only a third of the time: a `T>C`, half its reads converted.
    #[test]
    fn a_t_to_c_denovo_beta_counts_the_reads_that_read_as_t() {
        denovo_beta_counts_converted_reads(T, &[C, T, C, T]);
    }

    /// OT reads showing `bases` at a de-novo `reference>C` CpG; at least one
    /// must be the `C`, so that there is an alt.
    fn denovo_beta_counts_converted_reads(reference: Base, bases: &[Base]) {
        let reference = [reference, G];
        // The partner `G` reads as itself, so it has no alt and `beta_ratio`
        // is `beta_center`'s alone.
        let reads: Reads = bases.iter().map(|&at| ([at, G], Strand::OT)).collect();
        let converted = bases.iter().filter(|&&at| at == T).count() as f32;
        let expected = (converted / bases.len() as f32 + 1.0).log2();

        let [x, g] = columns(reference, &reads, None);
        let here = DenovoCpgFeatures::extract(&x.alt_metrics(C).expect("the alt"), None, Some(&g))
            .expect("features");
        assert_close("beta_ratio at X>C", here.beta_ratio, expected, &reads);

        let (mirrored_reference, mirrored) = mirror(reference, &reads);
        let [before, at] = columns(mirrored_reference, &mirrored, None);
        let there =
            DenovoCpgFeatures::extract(&at.alt_metrics(G).expect("the alt"), Some(&before), None)
                .expect("features");
        assert_close("beta_ratio at X>G", there.beta_ratio, expected, &reads);
    }

    /// `alt_score_adj` is the CpG partner's alt score, its reads at its own
    /// base qualities, as the de-novo model has it: this column's base
    /// qualities do not enter it. A strand mirror cannot show this: the
    /// mirror takes this column's qualities along.
    #[hegel::test]
    fn the_partner_alt_score_ignores_this_columns_qualities(tc: TestCase) {
        let reads = tc.draw(reads(T));
        let quality = tc.draw(gs::integers::<u8>().min_value(2).max_value(39));
        let at_c = |quality| {
            let [c, g] = columns([C, G], &reads, quality);
            CpgFeatures::extract(&c.alt_metrics(T).expect("the alt"), None, Some(&g))
                .expect("features")
                .alt_score_adj
        };
        assert_close("alt_score_adj", at_c(None), at_c(Some((0, quality))), &reads);

        let ([mc, mg], mirrored) = mirror([C, G], &reads);
        let at_g = |quality| {
            let [before, at] = columns([mc, mg], &mirrored, quality);
            CpgFeatures::extract(&at.alt_metrics(A).expect("the alt"), Some(&before), None)
                .expect("features")
                .alt_score_adj
        };
        assert_close("alt_score_adj", at_g(None), at_g(Some((1, quality))), &reads);
    }
}
