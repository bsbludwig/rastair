use super::{Fragment, Haplotype, Indel, LikelihoodParams, Reference, Site};
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, clap::Args, serde::Serialize)]
pub struct ExpansionParams {
    #[arg(skip)]
    #[serde(skip)]
    pub(crate) read_local: bool,
    /// Initial SNP discovery radius around the focal position
    #[arg(long, default_value_t = 150, value_parser = clap::value_parser!(u32).range(1..=10_000))]
    pub window_radius: u32,
    /// Maximum distance from the focal SNP reached through fragment links
    #[arg(long, default_value_t = 1000, value_parser = clap::value_parser!(u32).range(1..=10_000))]
    pub expansion_radius: u32,
    /// Safety limit on SNP sites, including the focal site
    #[arg(long, default_value_t = 32, value_parser = clap::value_parser!(u16).range(1..=64))]
    pub max_sites: u16,
    /// Maximum haplotypes retained for diploid inference
    #[arg(long, default_value_t = 64, value_parser = clap::value_parser!(u16).range(16..=256))]
    pub max_haplotypes: u16,
    /// Pruned haplotypes buffered for reconsideration after the next extension
    #[arg(long, default_value_t = 32, value_parser = clap::value_parser!(u16).range(0..=256))]
    pub haplotype_holdout: u16,
    /// Minimum independent fragments linking two candidate positions
    #[arg(long, default_value_t = 2, value_parser = clap::value_parser!(u16).range(1..))]
    pub min_linking_fragments: u16,
}

impl Default for ExpansionParams {
    fn default() -> Self {
        Self {
            read_local: false,
            window_radius: 150,
            expansion_radius: 1000,
            max_sites: 32,
            max_haplotypes: 64,
            haplotype_holdout: 32,
            min_linking_fragments: 2,
        }
    }
}

impl ExpansionParams {
    pub fn validate(&self) -> Result<()> {
        ensure!((1..=10_000).contains(&self.window_radius), "Initial radius must be 1..=10000");
        ensure!(
            (self.window_radius..=10_000).contains(&self.expansion_radius),
            "Expansion radius must be at least the initial radius and at most 10000"
        );
        ensure!((1..=64).contains(&self.max_sites), "Site safety limit must be 1..=64");
        ensure!((16..=256).contains(&self.max_haplotypes), "Haplotype budget must be 16..=256");
        ensure!(self.haplotype_holdout <= 256, "Haplotype holdout budget must be at most 256");
        ensure!(self.min_linking_fragments > 0, "Fragment link threshold must be positive");
        Ok(())
    }
}

pub(super) struct Expanded {
    pub sites: Vec<Site>,
    pub haplotypes: Vec<Haplotype>,
    pub disconnected_sites: usize,
    pub resource_limited_sites: usize,
    pub pruned_haplotypes: usize,
    pub rescued_haplotypes: usize,
    pub maximum_pruned_mass: f64,
    pub steps: usize,
}

#[expect(
    clippy::too_many_arguments,
    reason = "The search and scoring boundary share validated evidence and configuration"
)]
pub(super) fn expand(
    reference: &Reference,
    focal: Site,
    candidates: &[Site],
    indel: Option<&Indel>,
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
    require_links: bool,
) -> Result<Expanded> {
    expand_indels(
        reference,
        focal,
        candidates,
        indel.map(std::slice::from_ref).unwrap_or_default(),
        fragments,
        params,
        search,
        require_links,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Shared bounded search for SNP and competing-indel hypotheses"
)]
pub(super) fn expand_indels(
    reference: &Reference,
    focal: Site,
    candidates: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
    require_links: bool,
) -> Result<Expanded> {
    let focal_pos = focal.pos;
    let mut ordered = candidates.to_vec();
    ordered.sort_by_key(|site| {
        (
            site.pos.abs_diff(focal_pos) > search.window_radius,
            site.pos.abs_diff(focal_pos),
            site.pos,
        )
    });
    let links = if require_links {
        links(reference, &focal, &ordered, indels, fragments, params)?
    } else {
        BTreeMap::new()
    };
    let mut reachable = BTreeSet::from([focal_pos]);
    if require_links {
        loop {
            let previous = reachable.len();
            for (&(a, b), &count) in &links {
                if count >= usize::from(search.min_linking_fragments)
                    && (reachable.contains(&a) || reachable.contains(&b))
                {
                    reachable.extend([a, b]);
                }
            }
            if reachable.len() == previous {
                break;
            }
        }
    } else {
        reachable.extend(ordered.iter().map(|site| site.pos));
    }
    let disconnected_sites = ordered.iter().filter(|site| !reachable.contains(&site.pos)).count();
    ordered.retain(|site| reachable.contains(&site.pos));
    let mut result = Expanded {
        haplotypes: focal
            .alleles
            .iter()
            .map(|&base| Haplotype { alleles: vec![(focal_pos, base)], indel: None })
            .collect(),
        sites: vec![focal],
        disconnected_sites,
        resource_limited_sites: 0,
        pruned_haplotypes: 0,
        rescued_haplotypes: 0,
        maximum_pruned_mass: 0.0,
        steps: 0,
    };
    for indel in indels {
        let alternates = result
            .haplotypes
            .iter()
            .filter(|haplotype| haplotype.indel.is_none())
            .cloned()
            .map(|mut haplotype| {
                haplotype.indel = Some(indel.clone());
                haplotype
            })
            .collect::<Vec<_>>();
        result.haplotypes.extend(alternates);
    }
    let mut held_out = Vec::<Haplotype>::new();
    while !ordered.is_empty() && result.sites.len() < usize::from(search.max_sites) {
        let next = ordered.iter().position(|site| {
            !require_links
                || result.sites.iter().any(|active| {
                    links
                        .get(&(active.pos.min(site.pos), active.pos.max(site.pos)))
                        .is_some_and(|&count| count >= usize::from(search.min_linking_fragments))
                })
        });
        let Some(next) = next else { break };
        let site = ordered.remove(next);
        let mut sites = result.sites.clone();
        sites.push(site.clone());
        sites.sort_by_key(|site| site.pos);
        let retained_count = result.haplotypes.len();
        let proposals = result
            .haplotypes
            .iter()
            .chain(&held_out)
            .flat_map(|haplotype| {
                site.alleles.iter().map(|&base| {
                    let mut haplotype = haplotype.clone();
                    haplotype.alleles.push((site.pos, base));
                    haplotype.alleles.sort_by_key(|&(pos, _)| pos);
                    haplotype
                })
            })
            .collect::<Vec<_>>();
        if proposals.len() <= usize::from(search.max_haplotypes) {
            result.haplotypes = proposals;
            held_out.clear();
        } else {
            // Recompute each molecular likelihood on the expanded context. Old
            // posteriors are never used as priors for already observed evidence.
            let matrix = super::likelihood_matrix_geometry(
                reference,
                &sites,
                &proposals,
                indels.first(),
                fragments,
                params,
                search.read_local,
            )?;
            let pairs = super::posterior_pairs(
                &proposals,
                &sites,
                &matrix.values,
                params,
                !indels.is_empty(),
            )?;
            let mut support = vec![0.0; proposals.len()];
            for pair in &pairs {
                *support.get_mut(pair.first).ok_or_eyre("Missing first pruning candidate")? +=
                    *pair.posterior;
                *support.get_mut(pair.second).ok_or_eyre("Missing second pruning candidate")? +=
                    *pair.posterior;
            }
            let mut ranking = (0..proposals.len()).collect::<Vec<_>>();
            ranking.sort_by(|&a, &b| {
                support
                    .get(b)
                    .copied()
                    .unwrap_or(0.0)
                    .total_cmp(&support.get(a).copied().unwrap_or(0.0))
                    .then(a.cmp(&b))
            });
            let mut keep = BTreeSet::new();
            if let Some(index) = proposals.iter().position(|haplotype| {
                haplotype.indel.is_none()
                    && sites
                        .iter()
                        .all(|site| haplotype.allele_at(site.pos) == Some(site.reference))
            }) {
                keep.insert(index);
            }
            // Protect every focal genotype, not only the most likely chromosome.
            // This includes compound hets and fully converted latent C/G alleles.
            let mut best = [None::<usize>; 10];
            for (index, pair) in pairs.iter().enumerate() {
                let (a, b) = super::pair_bases(pair, &proposals, focal_pos)?;
                let genotype =
                    super::genotype_index(a, b).ok_or_eyre("Unknown protected genotype")?;
                let current = best.get_mut(genotype).ok_or_eyre("Invalid protected genotype")?;
                if current.and_then(|index| pairs.get(index)).is_none_or(|other| {
                    pair.log_likelihood + pair.log_prior > other.log_likelihood + other.log_prior
                }) {
                    *current = Some(index);
                }
            }
            for index in best.into_iter().flatten() {
                let pair = pairs.get(index).ok_or_eyre("Missing protected pair")?;
                keep.extend([pair.first, pair.second]);
            }
            for active in &sites {
                for &base in &active.alleles {
                    if let Some(index) = ranking.iter().copied().find(|&index| {
                        proposals
                            .get(index)
                            .is_some_and(|haplotype| haplotype.allele_at(active.pos) == Some(base))
                    }) {
                        keep.insert(index);
                    }
                }
            }
            if !indels.is_empty()
                && let Some(index) = ranking.iter().copied().find(|&index| {
                    proposals.get(index).is_some_and(|haplotype| haplotype.indel.is_some())
                })
            {
                keep.insert(index);
            }
            if indels.len() > 1 {
                let mut best = BTreeMap::new();
                for (index, pair) in pairs.iter().enumerate() {
                    let a = proposals
                        .get(pair.first)
                        .ok_or_eyre("Missing joint indel haplotype")?
                        .indel
                        .as_ref();
                    let b = proposals
                        .get(pair.second)
                        .ok_or_eyre("Missing joint indel haplotype")?
                        .indel
                        .as_ref();
                    let key = (a.min(b), a.max(b));
                    let previous = best.entry(key).or_insert(index);
                    if pairs.get(*previous).is_some_and(|other| {
                        pair.log_likelihood + pair.log_prior
                            > other.log_likelihood + other.log_prior
                    }) {
                        *previous = index;
                    }
                }
                for index in best.into_values() {
                    let pair = pairs.get(index).ok_or_eyre("Missing protected indel genotype")?;
                    keep.extend([pair.first, pair.second]);
                }
            }
            if keep.len() > usize::from(search.max_haplotypes) {
                ordered.push(site);
                break;
            }
            for &index in &ranking {
                if keep.len() == usize::from(search.max_haplotypes) {
                    break;
                }
                keep.insert(index);
            }
            let retained_mass = pairs
                .iter()
                .filter(|pair| keep.contains(&pair.first) && keep.contains(&pair.second))
                .map(|pair| *pair.posterior)
                .sum::<f64>();
            result.maximum_pruned_mass =
                result.maximum_pruned_mass.max((1.0 - retained_mass).clamp(0.0, 1.0));
            result.pruned_haplotypes += proposals.len() - keep.len();
            result.rescued_haplotypes +=
                keep.iter().filter(|&&index| index >= retained_count * site.alleles.len()).count();
            held_out = ranking
                .iter()
                .copied()
                .filter(|index| !keep.contains(index))
                .take(usize::from(search.haplotype_holdout))
                .map(|index| proposals.get(index).cloned().ok_or_eyre("Missing holdout candidate"))
                .collect::<Result<_>>()?;
            result.haplotypes = keep
                .into_iter()
                .map(|index| proposals.get(index).cloned().ok_or_eyre("Missing retained candidate"))
                .collect::<Result<_>>()?;
        }
        result.sites = sites;
        result.steps += 1;
    }
    result.resource_limited_sites = ordered.len();
    Ok(result)
}

fn links(
    reference: &Reference,
    focal: &Site,
    candidates: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
) -> Result<BTreeMap<(u32, u32), usize>> {
    let mut sites = candidates.iter().collect::<Vec<_>>();
    sites.push(focal);
    sites.sort_by_key(|site| site.pos);
    let mut context_pairs = BTreeSet::new();
    for site in &sites {
        if !site.alleles.contains(&seqair_types::Base::C) {
            continue;
        }
        for edit in std::iter::once(None).chain(indels.iter().map(Some)) {
            if edit.is_some_and(|edit| {
                !edit.inserted.is_empty() && site.pos.checked_add(1) == Some(edit.start)
            }) {
                continue;
            }
            let haplotype = Haplotype {
                alleles: vec![(site.pos, seqair_types::Base::C)],
                indel: edit.cloned(),
            };
            if let Some((partner, _)) = haplotype.adjacent_base(reference, site.pos, true)
                && sites.iter().any(|site| {
                    site.pos == partner && site.alleles.contains(&seqair_types::Base::G)
                })
            {
                context_pairs.insert((site.pos.min(partner), site.pos.max(partner)));
            }
        }
    }
    let mut links = BTreeMap::<(u32, u32), usize>::new();
    let mut context_support = BTreeMap::<(u32, u32), usize>::new();
    let mut genetic_alternates = BTreeMap::<u32, usize>::new();
    for fragment in fragments {
        let mut covered = Vec::new();
        for site in &sites {
            let pos = site.pos;
            let base = super::likelihood::consensus_at(fragment, pos)?;
            let Some(read) = base.source.and_then(|index| fragment.reads.get(index)) else {
                continue;
            };
            if !base.quality.is_some_and(|quality| quality >= 20)
                || !read.mapq.is_some_and(|mapq| mapq >= 20)
            {
                continue;
            }
            let conversion_side = (!params.non_taps)
                .then(|| {
                    [
                        crate::metrics::methylation::CpgSide::C,
                        crate::metrics::methylation::CpgSide::G,
                    ]
                    .into_iter()
                    .find(|side| {
                        if !site.alleles.contains(&side.unmod_base())
                            || (read.strand != side.strand()
                                && read.strand != seqair_types::Strand::Unknown)
                        {
                            return false;
                        }
                        std::iter::once(None).chain(indels.iter().map(Some)).any(|edit| {
                            let haplotype = Haplotype {
                                alleles: vec![(pos, side.unmod_base())],
                                indel: edit.cloned(),
                            };
                            haplotype
                                .adjacent_base(
                                    reference,
                                    pos,
                                    *side == crate::metrics::methylation::CpgSide::C,
                                )
                                .is_some_and(|(adjacent, base)| {
                                    let partner =
                                        if *side == crate::metrics::methylation::CpgSide::C {
                                            seqair_types::Base::G
                                        } else {
                                            seqair_types::Base::C
                                        };
                                    base == partner
                                        || sites.iter().any(|site| {
                                            site.pos == adjacent && site.alleles.contains(&partner)
                                        })
                                })
                        })
                    })
                })
                .flatten();
            if base.base != site.reference
                && !conversion_side.is_some_and(|side| {
                    site.reference == side.unmod_base() && base.base == side.mod_base()
                })
            {
                *genetic_alternates.entry(pos).or_default() += 1;
            }
            let confounded = conversion_side.is_some_and(|side| {
                base.base == side.mod_base()
                    && site
                        .alleles
                        .iter()
                        .all(|&base| base == side.unmod_base() || base == side.mod_base())
            });
            // A conversion-only observation cannot establish genetic phase. It
            // remains available to the rate fit and the final fragment likelihood.
            covered.push((pos, confounded));
        }
        for (index, &(a, confounded_a)) in covered.iter().enumerate() {
            for &(b, confounded_b) in covered.iter().skip(index + 1) {
                if !confounded_a && !confounded_b {
                    *links.entry((a, b)).or_default() += 1;
                }
            }
        }
        for &(a, b) in &context_pairs {
            if covered.iter().any(|&(pos, _)| pos == a || pos == b) {
                *context_support.entry((a, b)).or_default() += 1;
            }
        }
    }
    // Conversion alone must not turn every methylated CpG into a genetic
    // bridge. The focal target is always evaluated, including latent alleles.
    links.retain(|&(a, b), _| {
        // CpG partners must be inferred together: fixing an ambiguous partner
        // to reference can incorrectly preserve a CpG destroyed by a SNP.
        context_pairs.contains(&(a, b))
            || [a, b].into_iter().all(|pos| {
                pos == focal.pos || genetic_alternates.get(&pos).is_some_and(|&count| count >= 2)
            })
    });
    for (pair, count) in context_support {
        links.entry(pair).and_modify(|current| *current = (*current).max(count)).or_insert(count);
    }
    Ok(links)
}
