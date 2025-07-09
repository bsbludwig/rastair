use crate::{
    bed::Rastair1BedFormat,
    utils::Base::*,
    vcf::{AlleleSpecificStrandBias, ByStrand, GenotypeConfidence, GenotypeLikelihood, Methylated},
};
use color_eyre::eyre::{Context as _, ContextCompat as _, Report, bail};
use rastair2_vcf::{
    VcfField as _,
    standard_fields::{Genotype, ReadDepth},
};
use rust_htslib::bcf::Record as HtslibRecord;
use smallvec::{SmallVec, smallvec_inline};
use smol_str::SmolStr;
use tracing::trace;

impl TryFrom<&HtslibRecord> for Rastair1BedFormat {
    type Error = Report;

    #[allow(clippy::cast_possible_truncation)]
    fn try_from(r: &HtslibRecord) -> Result<Self, Self::Error> {
        let contig = r
            .rid()
            .wrap_err("Record has not ID")
            .and_then(|id| r.header().rid2name(id).wrap_err("Header does not contain ID"))
            .and_then(|name| str::from_utf8(name).wrap_err("Contig name is not valid UTF-8"))
            .map(SmolStr::new)
            .wrap_err("Could not fetch contig name")?;
        let r#ref = r
            .alleles()
            .first()
            .and_then(|x| str::from_utf8(x).ok())
            .map(SmolStr::new)
            .wrap_err("No reference base found in record")?;
        let beta = r
            .format(Methylated::ID.as_bytes())
            .float()
            .ok()
            .and_then(|x| x.first().copied())
            .and_then(|x| x.first())
            .copied()
            .unwrap_or_default();
        let read_depth = r
            .info(ReadDepth::ID.as_bytes())
            .integer()
            .wrap_err("Could not fetch read depth from record")?
            .and_then(|x| x.first())
            .copied()
            .unwrap_or_default() as usize;
        let alleles = r
            .alleles()
            .iter()
            .map(|a| str::from_utf8(a).map(SmolStr::new))
            .collect::<Result<SmallVec<_, 4>, _>>()
            .wrap_err("Failed to parse alleles")?;
        let assb = r
            .info(AlleleSpecificStrandBias::ID.as_bytes())
            .integer()
            .wrap_err("Failed to fetch AS_SB field")?
            .wrap_err("AS_SB field not set")?;

        let count = StrandCount::from_alleles_and_assb(&alleles, &assb)
            .wrap_err("Failed to parse strand counts from record")?;
        let (unmod, r#mod, no_snp, snp) = if r#ref == "C" {
            (count.c.ot, count.t.ot, count.c.ob, count.t.ob)
        } else if r#ref == "G" {
            (count.g.ob, count.a.ob, count.g.ot, count.a.ot)
        } else {
            (0, 0, 0, 0)
        };

        let genotype = r
            .genotypes()
            .map(|g| g.get(0).iter().map(|x| (*x).into()).collect())
            .unwrap_or_default();
        let genotype_likelihood = GenotypeLikelihood(smallvec_inline![
            r.format(GenotypeLikelihood::ID.as_bytes())
                .float()
                .ok()
                .and_then(|x| x.first().copied())
                .and_then(|x| x.first())
                .copied()
                .map(f64::from)
        ]);
        let genotype_confidence = GenotypeConfidence(smallvec_inline![
            r.format(GenotypeConfidence::ID.as_bytes())
                .float()
                .ok()
                .and_then(|x| x.first().copied())
                .and_then(|x| x.first())
                .copied()
                .map(f64::from),
        ]);

        Ok(Rastair1BedFormat {
            contig,
            pos: r.pos() as usize,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: read_depth,
            genotype: Genotype(genotype),
            genotype_likelihood,
            genotype_confidence,
        })
    }
}

#[derive(Debug, Default)]
struct StrandCount {
    a: ByStrand<u32>,
    c: ByStrand<u32>,
    g: ByStrand<u32>,
    t: ByStrand<u32>,
}

impl StrandCount {
    fn from_alleles_and_assb(alleles: &[SmolStr], assb: &[i32]) -> Result<StrandCount, Report> {
        // AS_SB is encoded with two integers per allele, so we need to parse it accordingly
        if assb.len() % 2 != 0 {
            bail!("AS_SB field has an odd number of integers");
        }
        let mut counts = StrandCount::default();
        for (i, count) in assb.chunks(2).enumerate() {
            if i >= alleles.len() {
                bail!("AS_SB field has more counts than alleles");
            }
            let allele = &alleles[i];
            let count = ByStrand {
                base: match allele.as_str() {
                    "A" => A,
                    "C" => C,
                    "G" => G,
                    "T" => T,
                    _ => {
                        trace!(%allele, "Unknown allele in AS_SB field");
                        continue;
                    }
                },
                ot: count[0] as u32,
                ob: count[1] as u32,
            };
            match count.base {
                A => counts.a = count,
                C => counts.c = count,
                G => counts.g = count,
                T => counts.t = count,
                _ => {}
            }
        }

        Ok(counts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use color_eyre::eyre::Result;
    use insta::assert_debug_snapshot;

    #[test]
    fn test_strand_count_from_record() -> Result<()> {
        let alleles = &[SmolStr::new_static("G"), SmolStr::new_static("A")];
        let assb = &[19, 1, 0, 16];
        let count = StrandCount::from_alleles_and_assb(alleles, assb)?;
        assert_debug_snapshot!(count, @r"
        StrandCount {
            a: ByStrand {
                base: A,
                ot: 0,
                ob: 16,
            },
            c: ByStrand {
                base: N,
                ot: 0,
                ob: 0,
            },
            g: ByStrand {
                base: G,
                ot: 19,
                ob: 1,
            },
            t: ByStrand {
                base: N,
                ot: 0,
                ob: 0,
            },
        }
        ");

        Ok(())
    }
}
