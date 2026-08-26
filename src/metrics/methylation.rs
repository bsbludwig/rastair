use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, DenovoAdjecent, PileupMetrics, ReadKey},
    utils::{Base::*, IntoF64, logging::ThisIsABug},
    vcf::{CpgBeta, CpgOrigin, InCpG, Methylated},
};
use color_eyre::{Result, eyre::Context};
use seqair_types::{Base, Probability, SmallVec, Strand};
use tracing::instrument;

#[cfg(test)]
mod tests;

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig(), pos = current.pos()),
    name = "methylation_call"
)]
pub fn call(current: &PileupMetrics) -> Result<Option<Methylated>> {
    let mut betas: SmallVec<CpgBeta, 2> = SmallVec::new();

    if let Some(b) = compute_beta(current, CpgSide::C).wrap_err("C-side beta")? {
        betas.push(b);
    }
    if let Some(b) = compute_beta(current, CpgSide::G).wrap_err("G-side beta")? {
        betas.push(b);
    }

    Ok((!betas.is_empty()).then_some(Methylated(betas)))
}

/// Determine whether this position has a CpG allele on the given side,
/// and if so compute the methylation beta value.
fn compute_beta(record: &PileupMetrics, side: CpgSide) -> Result<Option<CpgBeta>> {
    let Some(origin) = cpg_origin(record, side) else { return Ok(None) };

    let (raw_mod, raw_unmod) = read_counts(record, side);
    let mod_count = raw_mod.f();
    let unmod_count = raw_unmod.f();

    if mod_count + unmod_count == 0. {
        return Ok(None);
    }

    let adjustment = genotype_adjustment(record, side, origin);
    let beta = adjusted_beta(mod_count, unmod_count, adjustment);

    Ok(Some(CpgBeta {
        origin,
        beta: Probability::new(beta).this_is_a_bug()?,
        mod_count: raw_mod,
        total_count: raw_mod + raw_unmod,
    }))
}

/// Which side of the CpG dinucleotide we are looking at.
#[derive(Debug, Clone, Copy)]
enum CpgSide {
    C,
    G,
}

impl CpgSide {
    fn strand(self) -> Strand {
        match self {
            CpgSide::C => Strand::OT,
            CpgSide::G => Strand::OB,
        }
    }

    /// The base that appears when methylated (T for C-side, A for G-side).
    fn mod_base(self) -> Base {
        match self {
            CpgSide::C => T,
            CpgSide::G => A,
        }
    }

    /// The base that appears when unmethylated (C for C-side, G for G-side).
    fn unmod_base(self) -> Base {
        match self {
            CpgSide::C => C,
            CpgSide::G => G,
        }
    }

    /// The required adjacent base (G after C, C before G).
    fn adjacent_base(self) -> Base {
        match self {
            CpgSide::C => G,
            CpgSide::G => C,
        }
    }
}

/// Determine the CpG origin for this position on the given side, or None if
/// this position doesn't have a CpG allele on that side.
fn cpg_origin(record: &PileupMetrics, side: CpgSide) -> Option<CpgOrigin> {
    let is_original = match side {
        CpgSide::C => {
            record.pos_metrics.cpg == InCpG::C
                || record.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC
        }
        CpgSide::G => {
            record.pos_metrics.cpg == InCpG::G
                || record.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG
        }
    };

    if is_original {
        return Some(CpgOrigin::Original);
    }

    let cpg_base = side.unmod_base();
    let has_denovo_alt =
        record.alts.iter().any(|a| a.base == cpg_base && a.call == AltCall::RealVariant);

    let adjacent_present = match side {
        CpgSide::C => record.context.after_1 == Some(G),
        CpgSide::G => record.context.before_1 == Some(C),
    };

    if has_denovo_alt && adjacent_present {
        return Some(CpgOrigin::DeNovo);
    }

    None
}

/// Look up the mod and unmod read counts for the given CpG side.
fn read_counts(record: &PileupMetrics, side: CpgSide) -> (u32, u32) {
    let strand = side.strand();
    let adj = side.adjacent_base();

    let (counts, mod_base, unmod_base) = match side {
        CpgSide::C => (&record.after_counts, T, C),
        CpgSide::G => (&record.before_counts, A, G),
    };

    let raw_mod = counts.get(ReadKey { strand, current: mod_base, adj });
    let raw_unmod = counts.get(ReadKey { strand, current: unmod_base, adj });
    (raw_mod, raw_unmod)
}

/// How the genotype at this position affects the beta calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypeAdjustment {
    /// Normal beta: mod / (mod + unmod).
    None,
    /// The confounding base (T for C-side, A for G-side) is present as a het
    /// allele, so some mod-base reads come from the SNP, not methylation.
    HetConfounded,
    /// Original ref base is fully replaced by a homozygous variant — the
    /// original CpG no longer exists on either chromosome.
    HomAlt,
}

fn genotype_adjustment(
    record: &PileupMetrics,
    side: CpgSide,
    origin: CpgOrigin,
) -> GenotypeAdjustment {
    let Some(gt) = record.pos_metrics.genotype else {
        return GenotypeAdjustment::None;
    };

    if origin == CpgOrigin::Original && gt.genotype.is_homozygous() && !gt.genotype.is_hom_ref() {
        return GenotypeAdjustment::HomAlt;
    }

    if gt.genotype.is_heterozygous() {
        let confounding = side.mod_base();
        let confounded = match origin {
            CpgOrigin::Original => het_alt_is_base(record, &gt.genotype, confounding),
            CpgOrigin::DeNovo => {
                record.ref_base() == confounding
                    || het_alt_is_base(record, &gt.genotype, confounding)
            }
        };
        if confounded {
            return GenotypeAdjustment::HetConfounded;
        }
    }

    GenotypeAdjustment::None
}

fn het_alt_is_base(record: &PileupMetrics, gt: &GenotypeTag, base: Base) -> bool {
    match *gt {
        GenotypeTag::RefHet(idx) => {
            let i = (idx.get() as usize).saturating_sub(1);
            record.alts.get(i).map(|a| a.base) == Some(base)
        }
        GenotypeTag::AltHet(a, b) => {
            let ai = (a.get() as usize).saturating_sub(1);
            let bi = (b.get() as usize).saturating_sub(1);
            record.alts.get(ai).map(|x| x.base) == Some(base)
                || record.alts.get(bi).map(|x| x.base) == Some(base)
        }
        _ => false,
    }
}

fn adjusted_beta(mod_count: f64, unmod_count: f64, adjustment: GenotypeAdjustment) -> f64 {
    match adjustment {
        GenotypeAdjustment::HomAlt => 0.0,
        GenotypeAdjustment::HetConfounded => {
            let total = mod_count + unmod_count;
            let excess_mod = (mod_count - total / 2.).max(0.0);
            excess_mod / (unmod_count + excess_mod)
        }
        GenotypeAdjustment::None => mod_count / (mod_count + unmod_count),
    }
}
