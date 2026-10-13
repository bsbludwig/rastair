use super::LikelihoodParams;
use crate::sequence::ReaderParams;
use clap::value_parser;
use clio::ClioPath;
use color_eyre::eyre::{Result, bail};

#[derive(Debug, clap::Args)]
pub struct HaplotypeDiagnosticParams {
    #[command(flatten)]
    pub input: ReaderParams,
    #[command(flatten)]
    pub likelihood: LikelihoodParams,
    #[command(flatten)]
    pub reads: DiagnosticReadParams,
    /// Write per-SNP diagnostic features as TSV (stdout by default)
    #[arg(short, long)]
    pub output: Option<ClioPath>,
    /// Also write the competing haplotype-pair likelihoods and posteriors
    #[arg(long)]
    pub pair_output: Option<ClioPath>,
    /// Score these positions regardless of observed alternates (plain BED, one-base intervals)
    #[arg(long, value_parser = value_parser!(ClioPath).exists().is_file())]
    pub sites_bed: Option<ClioPath>,
    /// Include neighboring SNP candidates within this many bases of the focal SNP
    #[arg(long, default_value_t = 100, value_parser = value_parser!(u32).range(1..=1000))]
    pub window_radius: u32,
    /// Maximum candidate sites including the focal SNP; all focal A/C/G/T alleles are evaluated
    #[arg(long, default_value_t = 3, value_parser = value_parser!(u16).range(1..=3))]
    pub max_sites: u16,
    /// Minimum fragment observations of a base for a position to seed a candidate
    #[arg(long, default_value_t = 1, value_parser = value_parser!(u32).range(1..))]
    pub min_alt_observations: u32,
    /// Bases per input chunk; SNP windows have independent padding across chunk boundaries
    #[arg(long, default_value_t = 1000, value_parser = value_parser!(u32).range(1..=100_000))]
    pub chunk_length: u32,
}

#[derive(Debug, clap::Args, serde::Serialize)]
pub struct DiagnosticReadParams {
    #[command(flatten)]
    pub quality: crate::call::variant_calling::QualityFilterParams,
    #[command(flatten)]
    pub read_flags: crate::call::variant_calling::ReadFlags,
    #[command(flatten)]
    pub read_masking: crate::call::variant_calling::ReadMaskParams,
    /// Accept single-end reads and assign OT/OB from alignment direction
    #[arg(long)]
    pub unpaired: bool,
    /// Infer OT/OB from mismatch motifs, as in call --guess-read-orientation
    #[arg(long)]
    pub guess_read_orientation: bool,
}

pub fn diagnose(params: &HaplotypeDiagnosticParams) -> Result<()> {
    params.likelihood.validate()?;
    #[cfg(feature = "experimental-seqair")]
    return backend::run(params);
    #[cfg(not(feature = "experimental-seqair"))]
    bail!("ml haplotypes needs the `experimental-seqair` backend")
}

#[cfg(feature = "experimental-seqair")]
mod backend {
    use super::*;
    use crate::{
        call::haplotypes::{
            Fragment, Observation, ReadEvidence, Reference, Site, SnpEvidence, infer_snp,
        },
        sequence::{ChunkRegion, ReferenceWindow, Region, SeqairReaders},
    };
    use color_eyre::eyre::{OptionExt as _, WrapErr as _, ensure};
    use seqair::{
        bam::cigar::CigarOpType,
        reader::{DepthLimit, SegmentOptions},
    };
    use seqair_types::{Base, Pos0, Strand};
    use std::{
        collections::{BTreeMap, BTreeSet},
        io::{BufRead as _, BufReader, BufWriter, Write},
        num::NonZeroU32,
        sync::Arc,
    };
    use tracing::{info, warn};

    const HEADER: &str = "chrom\tpos\tref\talt\tlog_bayes_factor\tposterior\tp_het\tp_hom_alt\texpected_fraction\tgenotype_entropy\tsite_log_bayes_factor\tlinkage_log_bayes_factor\tfragments\tobservations\tcontext_sites\thaplotypes\tomitted_sites\tcpg_context_limited\talignment_indel_reads";

    struct InputChunk {
        reference: Reference,
        fragments: Vec<Fragment>,
        candidates: Vec<Site>,
        neighbors: Vec<Site>,
        indel_intervals: Vec<(u32, u32)>,
    }

    #[derive(Debug, PartialEq, Eq)]
    enum Destination {
        Stdout,
        File(std::path::PathBuf),
    }

    fn destination(path: &ClioPath) -> Result<Destination> {
        if path.is_std() {
            Ok(Destination::Stdout)
        } else {
            Ok(Destination::File(destination_path(path.path())?))
        }
    }

    pub(super) fn run(params: &HaplotypeDiagnosticParams) -> Result<()> {
        let regions = params
            .input
            .regions
            .as_ref()
            .ok_or_eyre("The experimental diagnostic requires an explicit --region")?;
        let forced_sites = params.sites_bed.as_ref().map(read_sites).transpose()?;
        let output_destination =
            params.output.as_ref().map(destination).transpose()?.unwrap_or(Destination::Stdout);
        let pair_destination = params.pair_output.as_ref().map(destination).transpose()?;
        ensure!(
            pair_destination.as_ref() != Some(&output_destination),
            "Diagnostic outputs must have different destinations"
        );
        for destination in std::iter::once(&output_destination).chain(&pair_destination) {
            if let Destination::File(destination) = destination {
                for input in [&params.input.bam_file, &params.input.fasta_file]
                    .into_iter()
                    .chain(&params.sites_bed)
                {
                    if input.is_std() {
                        continue;
                    }
                    ensure!(
                        *destination != std::fs::canonicalize(input.path())?,
                        "A diagnostic output would overwrite an input file"
                    );
                }
            }
        }
        let mut readers = params.input.open_seqair()?;
        let output: Box<dyn Write> = match &params.output {
            Some(path) => {
                Box::new(path.clone().create().wrap_err("Failed to create SNP diagnostics")?)
            }
            None => Box::new(std::io::stdout()),
        };
        let mut output = BufWriter::new(output);
        let mut pair_output = params
            .pair_output
            .as_ref()
            .map(|path| path.clone().create().map(BufWriter::new))
            .transpose()?;
        writeln!(output, "#rastair_haplotype_model=aligned_snp_v1")?;
        writeln!(output, "#likelihood_parameters={}", serde_json::to_string(&params.likelihood)?)?;
        writeln!(output, "#read_filters={}", serde_json::to_string(&params.reads)?)?;
        writeln!(
            output,
            "#window_radius={};max_sites={};min_alt_observations={}",
            params.window_radius, params.max_sites, params.min_alt_observations
        )?;
        writeln!(
            output,
            "#methylation_prior=independent_uniform_per_haplotype_cpg;maximum_observed_cpgs_per_haplotype=1;quadrature_is_approximate"
        )?;
        writeln!(output, "{HEADER}")?;
        if let Some(writer) = &mut pair_output {
            writeln!(
                writer,
                "chrom\tpos\tsites\thaplotype_a\thaplotype_b\tlog_likelihood\tlog_prior\tposterior"
            )?;
        }
        let mut completed = 0usize;
        let mut seen = BTreeSet::new();
        for region in regions.regions() {
            let contig = region.chromosome.as_str();
            let length = readers.header().resolve_contig(contig)?.len;
            let start = region
                .start
                .map_or(1, |pos| pos.as_u64())
                .checked_sub(1)
                .ok_or_eyre("Regions are 1-based")?;
            let end = region.end.map_or(length, |pos| pos.as_u64());
            ensure!(start < end && end <= length, "Invalid region {contig}:{}-{end}", start + 1);
            let mut core_start = start;
            while core_start < end {
                let core_end = core_start.saturating_add(u64::from(params.chunk_length)).min(end);
                let padding = u64::from(params.window_radius) + 2;
                let fetch_start = core_start.saturating_sub(padding);
                let fetch_end = core_end.saturating_add(padding).min(length);
                let chunk =
                    read_chunk(&mut readers, contig, fetch_start, fetch_end, length, params)?;
                let foci = match &forced_sites {
                    Some(sites) => {
                        let mut foci = Vec::new();
                        if let Some(positions) = sites.get(contig) {
                            for &pos in positions
                                .range(u32::try_from(core_start)?..u32::try_from(core_end)?)
                            {
                                match chunk
                                    .reference
                                    .base_at(pos)
                                    .filter(|base| base.known_index().is_some())
                                {
                                    Some(base) => foci.push(Site::new(pos, base, Base::KNOWN)?),
                                    None => warn!(
                                        contig,
                                        pos, "Skipping a target without a known reference base"
                                    ),
                                }
                            }
                        }
                        foci
                    }
                    None => chunk
                        .candidates
                        .iter()
                        .filter(|site| {
                            u64::from(site.pos) >= core_start && u64::from(site.pos) < core_end
                        })
                        .cloned()
                        .collect(),
                };
                for focal in &foci {
                    if !seen.insert((region.chromosome.clone(), focal.pos)) {
                        continue;
                    }
                    let nearby: Vec<_> = chunk
                        .neighbors
                        .iter()
                        .filter(|site| {
                            site.pos != focal.pos
                                && site.pos.abs_diff(focal.pos) <= params.window_radius
                        })
                        .cloned()
                        .collect();
                    let evidence = infer_snp(
                        &chunk.reference,
                        focal.pos,
                        &nearby,
                        &chunk.fragments,
                        &params.likelihood,
                        usize::from(params.max_sites),
                    )?;
                    let site_only = if evidence.sites.len() > 1 {
                        Some(infer_snp(
                            &chunk.reference,
                            focal.pos,
                            &[],
                            &chunk.fragments,
                            &params.likelihood,
                            1,
                        )?)
                    } else {
                        None
                    };
                    let first = focal.pos.saturating_sub(params.window_radius);
                    let last = focal.pos.saturating_add(params.window_radius);
                    let indels = chunk
                        .indel_intervals
                        .iter()
                        .filter(|&&(start, end)| start <= last && end >= first)
                        .count();
                    write_snp(&mut output, contig, &evidence, site_only.as_ref(), indels)?;
                    if let Some(writer) = &mut pair_output {
                        write_pairs(writer, contig, &evidence)?;
                    }
                    completed += 1;
                }
                core_start = core_end;
            }
        }
        output.flush().wrap_err("Failed to flush SNP diagnostics")?;
        if let Some(mut writer) = pair_output {
            writer.flush().wrap_err("Failed to flush pair diagnostics")?;
        }
        info!(sites = completed, "Wrote experimental haplotype likelihood diagnostics");
        Ok(())
    }

    fn read_sites(path: &ClioPath) -> Result<BTreeMap<String, BTreeSet<u32>>> {
        let mut sites: BTreeMap<String, BTreeSet<u32>> = BTreeMap::new();
        for (index, line) in BufReader::new(path.clone().open()?).lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }
            let mut fields = line.split_whitespace();
            let contig = fields.next().ok_or_eyre("BED contig is missing")?;
            let start = fields
                .next()
                .ok_or_eyre("BED start is missing")?
                .parse::<u32>()
                .wrap_err_with(|| format!("Invalid BED start on line {}", index + 1))?;
            let end = fields
                .next()
                .ok_or_eyre("BED end is missing")?
                .parse::<u64>()
                .wrap_err_with(|| format!("Invalid BED end on line {}", index + 1))?;
            ensure!(
                end == u64::from(start) + 1,
                "--sites-bed requires one-base intervals; line {} spans {start}..{end}",
                index + 1
            );
            sites.entry(contig.to_string()).or_default().insert(start);
        }
        Ok(sites)
    }

    fn destination_path(path: &std::path::Path) -> Result<std::path::PathBuf> {
        if path.exists() {
            return Ok(std::fs::canonicalize(path)?);
        }
        let name = path.file_name().ok_or_eyre("Output path has no filename")?;
        let parent = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| std::path::Path::new("."));
        Ok(std::fs::canonicalize(parent)?.join(name))
    }

    fn read_chunk(
        readers: &mut SeqairReaders,
        contig: &str,
        start: u64,
        end: u64,
        length: u64,
        params: &HaplotypeDiagnosticParams,
    ) -> Result<InputChunk> {
        let region = ChunkRegion {
            region: Region { contig: contig.into(), start, end },
            last_position: length,
            overlap_start: 0,
            overlap_end: 0,
        };
        let segment = readers.segment(&region, 0)?;
        let reference = Reference {
            start: u32::try_from(start)?,
            bases: segment.sequence.iter().copied().map(Base::from).collect(),
        };
        let custom = readers.inner_mut().customize_mut();
        custom.reference = Some(ReferenceWindow {
            start: Pos0::try_from(start)?,
            bases: Arc::from(reference.bases.clone()),
        });
        custom.read_flags = params.reads.read_flags.clone();
        custom.unpaired = params.reads.unpaired;
        custom.read_masking = params.reads.read_masking.clone();
        custom.guess_orientation = params.reads.guess_read_orientation;
        let span = (Pos0::try_from(start)?
            ..=Pos0::try_from(end.checked_sub(1).ok_or_eyre("Empty inference interval")?)?)
            .into();
        let max_len =
            NonZeroU32::new(u32::try_from(end - start)?).ok_or_eyre("Empty inference interval")?;
        let segments: Vec<_> = readers
            .inner_mut()
            .segments((contig, span), SegmentOptions::new(max_len).without_byte_budget())?
            .collect();
        let mut named: BTreeMap<(Vec<u8>, Vec<u8>), Vec<ReadEvidence>> = BTreeMap::new();
        let mut independent = Vec::new();
        let mut indel_intervals = Vec::new();
        for segment in segments {
            let guard = readers.inner_mut().pileup(&segment, DepthLimit::Unlimited).run()?;
            let store = guard.store();
            for record in store.records() {
                if record.mapq < params.reads.quality.min_mapq {
                    continue;
                }
                let extra = record.extra(store)?;
                let strand = if params.reads.unpaired {
                    if record.flags.is_reverse() { Strand::OB } else { Strand::OT }
                } else {
                    extra.strand
                };
                let Some(mask) = params.reads.read_masking.keep_range(
                    strand,
                    record.flags.is_reverse(),
                    record.seq_len,
                ) else {
                    continue;
                };
                let bases = record.seq(store)?;
                let qualities = record.qual(store)?;
                let mut observations = Vec::new();
                let mut has_indel = false;
                let mut pos = record.pos.as_u64();
                let mut query = 0u32;
                for op in record.cigar(store)? {
                    let len = op.len();
                    match op.op_type() {
                        CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                            for offset in 0..len {
                                let qpos = query
                                    .checked_add(offset)
                                    .ok_or_eyre("Query position overflow")?;
                                let genomic = pos + u64::from(offset);
                                if genomic < start || genomic >= end || !mask.contains(&qpos) {
                                    continue;
                                }
                                let index = usize::try_from(qpos)?;
                                let quality = qualities
                                    .get(index)
                                    .ok_or_eyre("CIGAR exceeds quality sequence")?
                                    .get();
                                if quality.is_some_and(|q| q < params.reads.quality.min_baseq)
                                    || (quality.is_none() && params.reads.quality.min_baseq > 0)
                                {
                                    continue;
                                }
                                let base =
                                    *bases.get(index).ok_or_eyre("CIGAR exceeds base sequence")?;
                                if base.known_index().is_some() {
                                    observations.push(Observation {
                                        pos: u32::try_from(genomic)?,
                                        base,
                                        quality,
                                    });
                                }
                            }
                        }
                        CigarOpType::Insertion | CigarOpType::Deletion | CigarOpType::RefSkip => {
                            has_indel = true;
                        }
                        CigarOpType::Unknown(code) => bail!("Unsupported BAM CIGAR opcode {code}"),
                        _ => {}
                    }
                    if op.op_type().consumes_ref() {
                        pos += u64::from(len);
                    }
                    if op.op_type().consumes_query() {
                        query = query.checked_add(len).ok_or_eyre("Query position overflow")?;
                    }
                }
                ensure!(
                    query == record.seq_len,
                    "CIGAR query length does not match sequence length"
                );
                if has_indel {
                    indel_intervals.push((
                        u32::try_from(record.pos.as_u64())?,
                        u32::try_from(record.end_pos.as_u64())?,
                    ));
                }
                if observations.is_empty() {
                    continue;
                }
                let read = ReadEvidence {
                    mapq: (record.mapq != 255).then_some(record.mapq),
                    strand,
                    observations,
                };
                let name = record.qname(store)?;
                if name.is_empty()
                    || name == b"*"
                    || params.reads.unpaired
                    || !record.flags.is_paired()
                {
                    independent.push(Fragment { reads: vec![read] });
                } else {
                    let group = match record.aux(store)?.get::<&[u8]>(*b"RG") {
                        Ok(group) => group.to_vec(),
                        Err(seqair::bam::aux::GetAuxError::TagNotFound { .. }) => Vec::new(),
                        Err(error) => return Err(error.into()),
                    };
                    named.entry((group, name.to_vec())).or_default().push(read);
                }
            }
        }
        for reads in named.into_values() {
            ensure!(
                reads.len() <= 2,
                "More than two eligible alignments share a read name and read group; the prototype cannot distinguish template collisions from additional alignments"
            );
            independent.push(Fragment { reads });
        }
        let mut counts: BTreeMap<u32, [u32; 4]> = BTreeMap::new();
        for fragment in &independent {
            let mut observed = BTreeSet::new();
            for obs in fragment.reads.iter().flat_map(|read| &read.observations) {
                if observed.insert((obs.pos, obs.base as u8))
                    && let Some(index) = obs.base.known_index()
                    && let Some(count) = counts.entry(obs.pos).or_default().get_mut(index)
                {
                    *count = count.checked_add(1).ok_or_eyre("Candidate depth overflow")?;
                }
            }
        }
        let mut candidates = Vec::new();
        let mut neighbors = Vec::new();
        for (pos, counts) in counts {
            let Some(reference_base) =
                reference.base_at(pos).filter(|base| base.known_index().is_some())
            else {
                continue;
            };
            let mut alts: Vec<_> = Base::KNOWN
                .into_iter()
                .filter(|&base| {
                    base != reference_base
                        && base
                            .known_index()
                            .and_then(|index| counts.get(index))
                            .is_some_and(|&count| count >= params.min_alt_observations)
                })
                .collect();
            let observed_alternate = !alts.is_empty();
            // Converted de-novo C/G alleles can look exactly like the reference.
            if !params.likelihood.non_taps {
                if reference_base == Base::T
                    && reference.base_at(pos.saturating_add(1)) == Some(Base::G)
                {
                    alts.push(Base::C);
                }
                if reference_base == Base::A
                    && pos.checked_sub(1).and_then(|pos| reference.base_at(pos)) == Some(Base::C)
                {
                    alts.push(Base::G);
                }
            }
            if !alts.is_empty() {
                let site = Site::new(pos, reference_base, alts)?;
                if observed_alternate {
                    neighbors.push(site.clone());
                }
                candidates.push(site);
            }
        }
        Ok(InputChunk { reference, fragments: independent, candidates, neighbors, indel_intervals })
    }

    fn write_snp(
        writer: &mut impl Write,
        contig: &str,
        evidence: &SnpEvidence,
        site_only: Option<&SnpEvidence>,
        indels: usize,
    ) -> Result<()> {
        for allele in &evidence.alleles {
            let baseline = site_only
                .and_then(|site| site.allele(allele.base))
                .unwrap_or(allele)
                .log_bayes_factor;
            writeln!(
                writer,
                "{contig}\t{}\t{}\t{}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{}\t{}\t{}\t{}\t{}\t{}\t{indels}",
                u64::from(evidence.pos) + 1,
                evidence.reference,
                allele.base,
                allele.log_bayes_factor,
                *allele.posterior,
                *allele.heterozygous,
                *allele.homozygous,
                allele.expected_fraction,
                evidence.genotype_entropy,
                baseline,
                allele.log_bayes_factor - baseline,
                evidence.focal_fragments,
                evidence.focal_observations,
                evidence.sites.len(),
                evidence.haplotypes.len(),
                evidence.omitted_sites,
                evidence.cpg_context_limited
            )?;
        }
        Ok(())
    }

    fn write_pairs(writer: &mut impl Write, contig: &str, evidence: &SnpEvidence) -> Result<()> {
        let sites = evidence
            .sites
            .iter()
            .map(|site| (u64::from(site.pos) + 1).to_string())
            .collect::<Vec<_>>()
            .join(",");
        for pair in &evidence.pairs {
            let sequence = |index: usize| -> Result<String> {
                Ok(evidence
                    .haplotypes
                    .get(index)
                    .ok_or_eyre("Missing diagnostic haplotype")?
                    .alleles
                    .iter()
                    .map(|(_, base)| char::from(*base as u8))
                    .collect())
            };
            writeln!(
                writer,
                "{contig}\t{}\t{sites}\t{}\t{}\t{:.12e}\t{:.12e}\t{:.12e}",
                u64::from(evidence.pos) + 1,
                sequence(pair.first)?,
                sequence(pair.second)?,
                pair.log_likelihood,
                pair.log_prior,
                *pair.posterior
            )?;
        }
        Ok(())
    }
}
