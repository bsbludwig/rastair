use super::{
    ExpansionParams, Fragment, Haplotype, Indel, LikelihoodParams, PairEvidence, Reference, Scorer,
    Site,
};
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
pub(super) struct Candidate {
    pub site: Site,
    pub indels: Vec<Indel>,
}

pub(super) struct Inferred {
    pub sites: Vec<Site>,
    pub haplotypes: Vec<Haplotype>,
    pub pairs: Vec<PairEvidence>,
    /// Each retained chromosome's log likelihood of every scored fragment.
    pub likelihoods: Vec<Vec<f64>>,
    pub used: usize,
    pub pruning_mass: f64,
    pub fragments: usize,
    pub realigned: usize,
}

fn compatible(first: &Indel, second: &Indel) -> bool {
    let end = |edit: &Indel| u64::from(edit.start) + edit.deleted.len().max(1) as u64;
    u64::from(first.start) > end(second) || u64::from(second.start) > end(first)
}

fn extend(haplotype: &Haplotype, candidate: &Candidate) -> Vec<Haplotype> {
    let mut proposals = Vec::new();
    for &base in &candidate.site.alleles {
        if haplotype.indels.iter().any(|edit| edit.deletes(candidate.site.pos)) {
            continue;
        }
        let mut next = haplotype.clone();
        next.alleles.push((candidate.site.pos, base));
        proposals.push(next.clone());
        for edit in &candidate.indels {
            if next.indels.iter().all(|other| compatible(edit, other))
                && next.alleles.iter().all(|(pos, _)| !edit.deletes(*pos))
            {
                let mut alternate = next.clone();
                alternate.indels.push(edit.clone());
                alternate.indels.sort();
                proposals.push(alternate);
            }
        }
    }
    proposals
}

/// A prefix is committed only when every site's competing genotypes fit the beam.
/// The caller starts a fresh region at the first uncommitted candidate.
pub(super) fn infer(
    reference: &Reference,
    candidates: &[Candidate],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
    scorer: &mut Scorer,
) -> Result<Inferred> {
    ensure!(!candidates.is_empty(), "Empty regional candidate set");
    let mut haplotypes = vec![Haplotype { alleles: Vec::new(), indels: Vec::new() }];
    let mut holdout = Vec::new();
    let mut sites = Vec::new();
    let mut pruning_mass = 0.0_f64;
    for candidate in candidates.iter().take(usize::from(search.max_sites)) {
        // A spanning deletion needs a separate allele state at an interior locus.
        // Split here rather than deleting that chromosome from the earlier posterior.
        if haplotypes
            .iter()
            .chain(&holdout)
            .any(|h| h.indels.iter().any(|e| e.deletes(candidate.site.pos)))
        {
            break;
        }
        let proposals = haplotypes
            .iter()
            .chain(&holdout)
            .flat_map(|h| extend(h, candidate))
            .collect::<Vec<_>>();
        ensure!(!proposals.is_empty(), "No compatible regional haplotypes");
        let mut offered_sites = sites.clone();
        offered_sites.push(candidate.site.clone());
        if proposals.len() <= usize::from(search.max_haplotypes) {
            haplotypes = proposals;
            holdout.clear();
            sites = offered_sites;
            continue;
        }
        let indel = proposals.iter().flat_map(|h| &h.indels).next();
        let matrix = super::likelihood_matrix_geometry(
            reference,
            &offered_sites,
            &proposals,
            indel,
            fragments,
            params,
            true,
            scorer,
        )?;
        let pairs = super::posterior_pairs(
            &proposals,
            &offered_sites,
            &matrix.values,
            params,
            indel.is_some(),
        )?;
        let mut support = vec![0.0; proposals.len()];
        let mut snps = BTreeMap::new();
        let mut indels = BTreeMap::new();
        let starts: BTreeSet<_> =
            proposals.iter().flat_map(|h| &h.indels).map(|e| e.start).collect();
        for (index, pair) in pairs.iter().enumerate() {
            let first = proposals.get(pair.first).ok_or_eyre("Missing regional chromosome")?;
            let second = proposals.get(pair.second).ok_or_eyre("Missing regional chromosome")?;
            for chromosome in [pair.first, pair.second] {
                *support.get_mut(chromosome).ok_or_eyre("Invalid regional support")? +=
                    *pair.posterior;
            }
            let better = |previous: &mut usize| {
                if pairs.get(*previous).is_some_and(|p| {
                    p.log_likelihood + p.log_prior < pair.log_likelihood + pair.log_prior
                }) {
                    *previous = index;
                }
            };
            for site in &offered_sites {
                let a = first.allele_at(site.pos).ok_or_eyre("Missing regional SNP")? as u8;
                let b = second.allele_at(site.pos).ok_or_eyre("Missing regional SNP")? as u8;
                better(snps.entry((site.pos, a.min(b), a.max(b))).or_insert(index));
            }
            for &start in &starts {
                let a = first.indels.iter().find(|edit| edit.start == start);
                let b = second.indels.iter().find(|edit| edit.start == start);
                better(indels.entry((start, a.min(b), a.max(b))).or_insert(index));
            }
        }
        let mut keep = BTreeSet::new();
        for index in snps.into_values().chain(indels.into_values()) {
            let pair = pairs.get(index).ok_or_eyre("Missing protected regional genotype")?;
            keep.extend([pair.first, pair.second]);
        }
        if let Some(index) = proposals.iter().position(|h| {
            h.indels.is_empty()
                && offered_sites.iter().all(|s| h.allele_at(s.pos) == Some(s.reference))
        }) {
            keep.insert(index);
        }
        if keep.len() > usize::from(search.max_haplotypes) {
            break;
        }
        let mut ranked = (0..proposals.len()).collect::<Vec<_>>();
        ranked.sort_by(|a, b| {
            support
                .get(*b)
                .copied()
                .unwrap_or_default()
                .total_cmp(&support.get(*a).copied().unwrap_or_default())
                .then(a.cmp(b))
        });
        for &index in &ranked {
            if keep.len() >= usize::from(search.max_haplotypes) {
                break;
            }
            keep.insert(index);
        }
        let retained = pairs
            .iter()
            .filter(|p| keep.contains(&p.first) && keep.contains(&p.second))
            .map(|p| *p.posterior)
            .sum::<f64>();
        pruning_mass = pruning_mass.max((1.0 - retained).clamp(0.0, 1.0));
        holdout = ranked
            .into_iter()
            .filter(|i| !keep.contains(i))
            .take(usize::from(search.haplotype_holdout))
            .map(|i| proposals.get(i).cloned().ok_or_eyre("Missing regional holdout"))
            .collect::<Result<_>>()?;
        haplotypes = keep
            .into_iter()
            .map(|i| proposals.get(i).cloned().ok_or_eyre("Missing retained chromosome"))
            .collect::<Result<_>>()?;
        sites = offered_sites;
    }
    ensure!(!sites.is_empty(), "First regional candidate exceeds haplotype budget");
    let indel = haplotypes.iter().flat_map(|h| &h.indels).next();
    let matrix = super::likelihood_matrix_geometry(
        reference,
        &sites,
        &haplotypes,
        indel,
        fragments,
        params,
        true,
        scorer,
    )?;
    let pairs =
        super::posterior_pairs(&haplotypes, &sites, &matrix.values, params, indel.is_some())?;
    Ok(Inferred {
        used: sites.len(),
        sites,
        haplotypes,
        pairs,
        likelihoods: matrix.values,
        pruning_mass,
        fragments: matrix.fragments,
        realigned: matrix.realigned_fragments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqair::bam::cigar::CigarOpType::{Deletion, Insertion, Match};
    use seqair_types::Base;

    #[test]
    fn independent_indel_loci_keep_their_marginal_priors() -> Result<()> {
        let reference = Reference { start: 1000, bases: vec![Base::A; 100] };
        let candidates: Vec<_> = [1010, 1050]
            .into_iter()
            .map(|start| Candidate {
                site: Site { pos: start - 1, reference: Base::A, alleles: vec![Base::A] },
                indels: vec![Indel { start, deleted: vec![], inserted: vec![Base::C] }],
            })
            .collect();
        let params = super::super::features::parameters()?;
        let result = infer(
            &reference,
            &candidates,
            &[],
            &params,
            &ExpansionParams::default(),
            &mut Scorer::new(&params)?,
        )?;
        for start in [1010, 1050] {
            let mut probabilities = [0.0_f64; 3];
            for pair in &result.pairs {
                let copies = [pair.first, pair.second]
                    .into_iter()
                    .filter(|&i| {
                        result
                            .haplotypes
                            .get(i)
                            .is_some_and(|h| h.indels.iter().any(|e| e.start == start))
                    })
                    .count();
                *probabilities.get_mut(copies).ok_or_eyre("indel copies")? += *pair.posterior;
            }
            let prior = *params.indel_prior;
            for (observed, expected) in
                probabilities.into_iter().zip([1.0 - prior, prior / 2.0, prior / 2.0])
            {
                assert!((observed - expected).abs() < 1e-10);
            }
        }
        Ok(())
    }

    #[test]
    fn spanning_deletion_splits_without_erasing_the_first_locus() -> Result<()> {
        let reference = Reference { start: 1000, bases: vec![Base::A; 100] };
        let deletion = Indel { start: 1010, deleted: vec![Base::A; 10], inserted: vec![] };
        let candidates = [
            Candidate {
                site: Site { pos: 1009, reference: Base::A, alleles: vec![Base::A] },
                indels: vec![deletion.clone()],
            },
            Candidate {
                site: Site { pos: 1015, reference: Base::A, alleles: vec![Base::A] },
                indels: vec![Indel { start: 1016, deleted: vec![], inserted: vec![Base::T] }],
            },
        ];
        let params = super::super::features::parameters()?;
        let result = infer(
            &reference,
            &candidates,
            &[],
            &params,
            &ExpansionParams::default(),
            &mut Scorer::new(&params)?,
        )?;
        assert_eq!(result.used, 1);
        assert!(result.haplotypes.iter().any(|h| h.indels.contains(&deletion)));
        Ok(())
    }

    #[test]
    fn a_pruned_region_keeps_every_local_genotype() -> Result<()> {
        let reference = Reference { start: 1000, bases: vec![Base::A; 200] };
        let candidates = (0..8)
            .map(|i| {
                Ok(Candidate {
                    site: Site::new(1020 + 10 * i, Base::A, Base::KNOWN)?,
                    indels: vec![],
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let params = super::super::features::parameters()?;
        let search = ExpansionParams {
            max_haplotypes: 16,
            haplotype_holdout: 4,
            ..ExpansionParams::default()
        };
        let mut offset = 0;
        let mut regions = 0;
        while offset < candidates.len() {
            let result = infer(
                &reference,
                candidates.get(offset..).ok_or_eyre("candidates")?,
                &[],
                &params,
                &search,
                &mut Scorer::new(&params)?,
            )?;
            assert!(result.haplotypes.len() <= 16);
            for site in &result.sites {
                let genotypes: BTreeSet<_> = result
                    .pairs
                    .iter()
                    .map(|p| {
                        let (a, b) = super::super::pair_bases(p, &result.haplotypes, site.pos)?;
                        Ok(((a as u8).min(b as u8), (a as u8).max(b as u8)))
                    })
                    .collect::<Result<_>>()?;
                assert_eq!(genotypes.len(), 10);
            }
            offset += result.used;
            regions += 1;
        }
        assert!(regions > 1);
        Ok(())
    }

    #[test]
    fn two_indels_and_a_snp_share_one_diploid_posterior() -> Result<()> {
        let mut state = 91217_u32;
        let dna: Vec<u8> = (0..300)
            .map(|_| {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                *b"ACGT".get((state >> 30) as usize).expect("base")
            })
            .collect();
        let reference =
            Reference { start: 1000, bases: dna.iter().copied().map(Base::from).collect() };
        let insertion =
            Indel { start: 1070, deleted: Vec::new(), inserted: vec![Base::A, Base::C] };
        let deletion = Indel {
            start: 1120,
            deleted: reference.bases.get(120..123).ok_or_eyre("deletion")?.to_vec(),
            inserted: Vec::new(),
        };
        let alternate = if reference.base_at(1090) == Some(Base::A) { Base::G } else { Base::A };
        let mut changed = dna.clone();
        *changed.get_mut(90).ok_or_eyre("SNP")? = alternate as u8;
        let spelled = [
            changed.get(..70).ok_or_eyre("left")?,
            b"AC",
            changed.get(70..120).ok_or_eyre("middle")?,
            changed.get(123..).ok_or_eyre("right")?,
        ]
        .concat();
        let alt = super::super::reads::fixture_read(
            &reference,
            &spelled,
            &[(Match, 70), (Insertion, 2), (Match, 50), (Deletion, 3), (Match, 177)],
        )?;
        let wt = super::super::reads::fixture_read(&reference, &dna, &[(Match, 300)])?;
        let fragments: Vec<_> = [alt, wt]
            .into_iter()
            .flat_map(|read| std::iter::repeat_n(Fragment { reads: vec![read] }, 20))
            .collect();
        let candidates = [
            Candidate {
                site: Site {
                    pos: 1069,
                    reference: reference.base_at(1069).ok_or_eyre("anchor")?,
                    alleles: vec![reference.base_at(1069).ok_or_eyre("anchor")?],
                },
                indels: vec![insertion.clone()],
            },
            Candidate {
                site: Site::new(1090, reference.base_at(1090).ok_or_eyre("SNP")?, [alternate])?,
                indels: vec![],
            },
            Candidate {
                site: Site {
                    pos: 1119,
                    reference: reference.base_at(1119).ok_or_eyre("anchor")?,
                    alleles: vec![reference.base_at(1119).ok_or_eyre("anchor")?],
                },
                indels: vec![deletion.clone()],
            },
        ];
        let mut params = super::super::features::parameters()?;
        params.non_taps = true;
        let search = ExpansionParams::default();
        let inferred = infer(
            &reference,
            &candidates,
            &fragments,
            &params,
            &search,
            &mut Scorer::new(&params)?,
        )?;
        assert_eq!(inferred.used, 3);
        let correct: f64 = inferred
            .pairs
            .iter()
            .filter(|pair| {
                let Some(a) = inferred.haplotypes.get(pair.first) else { return false };
                let Some(b) = inferred.haplotypes.get(pair.second) else { return false };
                [a, b].iter().any(|h| {
                    h.indels == vec![insertion.clone(), deletion.clone()]
                        && h.allele_at(1090) == Some(alternate)
                }) && [a, b]
                    .iter()
                    .any(|h| h.indels.is_empty() && h.allele_at(1090) == reference.base_at(1090))
            })
            .map(|pair| *pair.posterior)
            .sum();
        assert!(correct > 0.99, "{correct}");
        let duplicated: Vec<_> = fragments
            .iter()
            .map(|fragment| Fragment {
                reads: vec![fragment.reads.first().expect("read").clone(); 2],
            })
            .collect();
        let paired = infer(
            &reference,
            &candidates,
            &duplicated,
            &params,
            &search,
            &mut Scorer::new(&params)?,
        )?;
        for (a, b) in inferred.pairs.iter().zip(&paired.pairs) {
            assert!((*a.posterior - *b.posterior).abs() < 1e-6);
        }
        let converted = |dna: &[u8]| {
            let mut bases = dna.to_vec();
            for (index, context) in dna.windows(2).enumerate() {
                if context == b"CG" {
                    *bases.get_mut(index).expect("CpG coordinate") = b'T';
                }
            }
            bases
        };
        let modified = super::super::reads::fixture_read(
            &reference,
            &converted(&spelled),
            &[(Match, 70), (Insertion, 2), (Match, 50), (Deletion, 3), (Match, 177)],
        )?;
        let unmodified =
            super::super::reads::fixture_read(&reference, &converted(&dna), &[(Match, 300)])?;
        let fragments: Vec<_> = [modified, unmodified]
            .into_iter()
            .flat_map(|read| std::iter::repeat_n(Fragment { reads: vec![read] }, 20))
            .collect();
        params.non_taps = false;
        let taps = infer(
            &reference,
            &candidates,
            &fragments,
            &params,
            &search,
            &mut Scorer::new(&params)?,
        )?;
        let correct = taps
            .pairs
            .iter()
            .filter(|pair| {
                let Some(a) = taps.haplotypes.get(pair.first) else { return false };
                let Some(b) = taps.haplotypes.get(pair.second) else { return false };
                [a, b].iter().any(|h| {
                    h.indels == vec![insertion.clone(), deletion.clone()]
                        && h.allele_at(1090) == Some(alternate)
                }) && [a, b]
                    .iter()
                    .any(|h| h.indels.is_empty() && h.allele_at(1090) == reference.base_at(1090))
            })
            .map(|p| *p.posterior)
            .sum::<f64>();
        assert!(correct > 0.99, "TAPS multi-indel posterior: {correct}");
        Ok(())
    }

    #[test]
    fn a_stutter_mixture_reads_a_lone_reference_molecule_as_a_stuttered_alt() -> Result<()> {
        let mut dna = b"GTTCTGGAACGCGCTTCT".to_vec();
        dna.extend(std::iter::repeat_n(b'A', 12));
        dna.extend(b"GTTAGGTACCGT");
        let local = Reference { start: 1100, bases: dna.iter().copied().map(Base::from).collect() };
        let mut padded = vec![Base::C; 100];
        padded.extend(&local.bases);
        padded.extend(vec![Base::G; 100]);
        let reference = Reference { start: 1000, bases: padded };
        let alternate = [&dna[..18], &dna[19..]].concat();
        let alt = super::super::reads::fixture_read(
            &local,
            &alternate,
            &[(Match, 18), (Deletion, 1), (Match, 23)],
        )?;
        let wt = super::super::reads::fixture_read(&local, &dna, &[(Match, 42)])?;
        let edit = Indel { start: 1118, deleted: vec![Base::A], inserted: vec![] };
        let anchor = reference.base_at(1117).ok_or_eyre("anchor")?;
        let candidates = [Candidate {
            site: Site { pos: 1117, reference: anchor, alleles: vec![anchor] },
            indels: vec![edit],
        }];
        let mut params = super::super::features::parameters()?;
        params.non_taps = true;
        params.error_profile.stutter = (2..=20)
            .map(|copies| super::super::error_model::StutterRate {
                period: 1,
                minimum_copies: copies,
                steps: [0.01, 0.1, 0.05, 0.005],
            })
            .collect();
        params.validate()?;
        let posterior = |alt_molecules: usize, wt_molecules: usize| -> Result<[f64; 3]> {
            let mut fragments = vec![Fragment { reads: vec![alt.clone()] }; alt_molecules];
            fragments.extend(vec![Fragment { reads: vec![wt.clone()] }; wt_molecules]);
            let inferred = infer(
                &reference,
                &candidates,
                &fragments,
                &params,
                &ExpansionParams::default(),
                &mut Scorer::new(&params)?,
            )?;
            let mut copies = [0.0; 3];
            for pair in &inferred.pairs {
                let count = [pair.first, pair.second]
                    .into_iter()
                    .filter(|&i| inferred.haplotypes.get(i).is_some_and(|h| !h.indels.is_empty()))
                    .count();
                *copies.get_mut(count).ok_or_eyre("copies")? += *pair.posterior;
            }
            Ok(copies)
        };
        // Dosage, conditional on the allele: the per-locus prior still weighs
        // on presence, which the forest decides in production.
        let [_, het, hom] = posterior(20, 1)?;
        assert!(hom / (het + hom) > 0.95, "lone reference molecule: het {het} hom {hom}");
        let [_, het, hom] = posterior(10, 10)?;
        assert!(het / (het + hom) > 0.99, "balanced: het {het} hom {hom}");
        Ok(())
    }

    #[test]
    fn padded_insertion_and_deletion_models_tolerate_one_opposing_molecule() -> Result<()> {
        let short = b"GTTCTGGAACGCGCTTCTATTAGGTA";
        let long = b"GTTCTGGAACGCGCATCTCTTCTATTAGGTA";
        for insertion in [true, false] {
            let (dna, alternate): (&[u8], &[u8]) =
                if insertion { (short, long) } else { (long, short) };
            let local =
                Reference { start: 1100, bases: dna.iter().copied().map(Base::from).collect() };
            let mut padded = vec![Base::A; 100];
            padded.extend(&local.bases);
            padded.extend(vec![Base::G; 100]);
            let reference = Reference { start: 1000, bases: padded };
            let cigar = if insertion {
                vec![(Match, 13), (Insertion, 5), (Match, 13)]
            } else {
                vec![(Match, 13), (Deletion, 5), (Match, 13)]
            };
            let alt = super::super::reads::fixture_read(&local, alternate, &cigar)?;
            let wt = super::super::reads::fixture_read(
                &local,
                dna,
                &[(Match, u32::try_from(dna.len())?)],
            )?;
            let edit = Indel {
                start: 1113,
                deleted: if insertion {
                    vec![]
                } else {
                    b"CATCT".iter().copied().map(Base::from).collect()
                },
                inserted: if insertion {
                    b"CATCT".iter().copied().map(Base::from).collect()
                } else {
                    vec![]
                },
            };
            let anchor = reference.base_at(1112).ok_or_eyre("anchor")?;
            let candidates = [Candidate {
                site: Site { pos: 1112, reference: anchor, alleles: vec![anchor] },
                indels: vec![edit],
            }];
            let mut fragments = vec![Fragment { reads: vec![alt] }; 20];
            fragments.push(Fragment { reads: vec![wt] });
            let mut params = super::super::features::parameters()?;
            params.non_taps = true;
            let inferred = infer(
                &reference,
                &candidates,
                &fragments,
                &params,
                &ExpansionParams::default(),
                &mut Scorer::new(&params)?,
            )?;
            let hom: f64 = inferred
                .pairs
                .iter()
                .filter(|p| {
                    p.first == p.second
                        && inferred.haplotypes.get(p.first).is_some_and(|h| !h.indels.is_empty())
                })
                .map(|p| *p.posterior)
                .sum();
            assert!(hom > 0.95, "insertion={insertion}: {hom}");
        }
        Ok(())
    }
}
