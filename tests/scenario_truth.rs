//! The scenario as truth oracle: generated scenarios with planted SNVs,
//! indels, CpGs and de-novo CpGs, at a depth where the evidence leaves no
//! doubt, and what the caller must get right there.
//!
//! The reads are clean: base quality 40, mapping quality 60, both haplotypes
//! and both strands at equal depth, every read pair the same length. At 32x
//! a het SNV has 16 alt and 16 reference reads; calling it `0/0` or `1/1`
//! takes 16 sequencing errors at 10⁻⁴ each. So a genotype is asserted
//! exactly. Methylation is drawn per read pair, so β̂ at one side of a CpG is
//! binomial over the reads of the strand that shows it there, and is held to
//! five standard deviations of that binomial. A het de-novo CpG is left out:
//! there an OT `T` is the reference's or a methylated `C`'s, so β rests on the
//! genotype and is no binomial of the reads.

#[path = "../src/pinned.rs"]
mod pinned;
mod utils;

use hegel::TestCase;
use utils::{
    planted::{
        self, Carrier, Event, FIRST_EVENT, IndelFlanks, Kinds, Plan, PlantedCpg, PlantedIndel,
        PlantedSnv,
    },
    *,
};

/// Plans at 32x or 64x whose indels make no CpG (see
/// [`an_insertion_that_makes_a_cpg_is_no_snv`]).
#[hegel::composite]
fn plans(tc: &TestCase, indels: bool) -> Plan {
    let kinds = Kinds { indels, max_events: 6, indel_flanks: IndelFlanks::Neutral };
    Plan { step: 25, ..tc.draw(planted::plans(kinds)) }
}

/// The columns of a VCF data line the oracle reads.
#[derive(Debug)]
struct Record {
    /// 1-based.
    pos: usize,
    reference: String,
    alts: Vec<String>,
    filter: String,
    line: String,
}

impl Record {
    fn parse(line: &str) -> Result<Self> {
        let columns: Vec<&str> = line.split('\t').collect();
        let column = |index: usize| -> Result<&str> {
            columns.get(index).copied().ok_or_else(|| eyre!("short line {line}"))
        };
        Ok(Self {
            pos: column(1)?.parse()?,
            reference: column(3)?.to_owned(),
            alts: column(4)?.split(',').map(str::to_owned).collect(),
            filter: column(6)?.to_owned(),
            line: line.to_owned(),
        })
    }

    fn field(&self, key: &str) -> Option<&str> {
        format_field(&self.line, key)
    }

    /// `GT` as sorted allele numbers; `None` for a missing allele.
    fn genotype(&self) -> Vec<Option<usize>> {
        let mut alleles: Vec<Option<usize>> = self
            .field("GT")
            .unwrap_or(".")
            .split(['/', '|'])
            .map(|allele| allele.parse().ok())
            .collect();
        alleles.sort_unstable();
        alleles
    }

    fn carries_an_alt(&self) -> bool {
        self.genotype().iter().any(|allele| allele.is_some_and(|allele| allele > 0))
    }

    fn passes(&self) -> bool {
        self.filter == "PASS"
    }

    /// The reference with this record's first alt applied.
    fn applied_to(&self, reference: &[u8]) -> Option<Vec<u8>> {
        let start = self.pos.checked_sub(1)?;
        let alt = self.alts.first()?;
        let mut sequence = reference.get(..start)?.to_vec();
        sequence.extend_from_slice(alt.as_bytes());
        sequence.extend_from_slice(reference.get(start + self.reference.len()..)?);
        Some(sequence)
    }

    fn is_indel(&self) -> bool {
        self.reference.len() != 1 || self.alts.iter().any(|alt| alt.len() != 1 && alt != ".")
    }
}

fn call(plan: &Plan) -> Result<Vec<Record>> {
    let vcf = call_scenario(&plan.scenario(), &["--vcf-all-fields"])?;
    vcf_content_lines(&vcf).map(Record::parse).collect()
}

/// The genotype a carrier gives an alt numbered `allele`.
fn genotype_of(carrier: Carrier, allele: usize) -> Vec<Option<usize>> {
    match carrier {
        Carrier::One | Carrier::Two => vec![Some(0), Some(allele)],
        Carrier::Both => vec![Some(allele), Some(allele)],
    }
}

fn indel_applied(reference: &[u8], indel: &PlantedIndel) -> Option<Vec<u8>> {
    let mut sequence = reference.get(..=indel.anchor)?.to_vec();
    sequence.extend_from_slice(&indel.inserted);
    sequence.extend_from_slice(reference.get(indel.anchor + 1 + indel.deleted..)?);
    Some(sequence)
}

fn check_snvs(plan: &Plan, records: &[Record]) -> Result<()> {
    for snv @ PlantedSnv { pos, reference, alt, carrier } in plan.snvs() {
        let alt = char::from(alt).to_string();
        let record = records
            .iter()
            .find(|record| record.pos == pos + 1 && record.passes() && record.alts.contains(&alt))
            .ok_or_else(|| eyre!("{snv:?} is not called"))?;
        ensure!(record.reference == char::from(reference).to_string(), "{snv:?}: {record:?}");
        let allele = record.alts.iter().position(|a| *a == alt).map_or(0, |index| index + 1);
        ensure!(
            record.genotype() == genotype_of(carrier, allele),
            "{snv:?} is called {:?}: {}",
            record.field("GT"),
            record.line
        );
    }
    Ok(())
}

/// A het or hom SNV on clean reads at 32x or 64x is called PASS with its
/// genotype, de-novo CpGs and transitions TAPS confounds with methylation
/// included.
#[hegel::test(test_cases = pinned::cases(4))]
fn planted_snvs_are_called_with_their_genotype(tc: TestCase) {
    let plan = tc.draw(plans(false));
    let result = call(&plan).and_then(|records| check_snvs(&plan, &records));
    assert!(result.is_ok(), "{result:?}");
}

fn check_indels(plan: &Plan, records: &[Record]) -> Result<()> {
    let reference = plan.reference();
    let reference = reference.reference();
    for indel in plan.indels() {
        let planted = indel_applied(reference, &indel).ok_or_else(|| eyre!("{indel:?}"))?;
        let record = records
            .iter()
            .filter(|record| record.passes() && record.is_indel())
            .find(|record| record.applied_to(reference).as_ref() == Some(&planted))
            .ok_or_else(|| eyre!("{indel:?} is not called"))?;
        ensure!(
            record.genotype() == genotype_of(indel.carrier, 1),
            "{indel:?} is called {:?}: {}",
            record.field("GT"),
            record.line
        );
    }
    Ok(())
}

/// A het or hom indel of up to three bases on clean reads at 32x or 64x is
/// called PASS with its genotype, wherever the caller places it in a repeat.
#[hegel::test(test_cases = pinned::cases(4))]
fn planted_indels_are_called_with_their_genotype(tc: TestCase) {
    let plan = tc.draw(plans(true));
    let result = call(&plan).and_then(|records| check_indels(&plan, &records));
    assert!(result.is_ok(), "{result:?}");
}

/// Every PASS record that carries an alt is a planted SNV or indel.
fn check_no_false_calls(plan: &Plan, records: &[Record]) -> Result<()> {
    let reference = plan.reference();
    let reference = reference.reference();
    let snvs = plan.snvs();
    let indels: Vec<Vec<u8>> =
        plan.indels().iter().filter_map(|indel| indel_applied(reference, indel)).collect();
    for record in records.iter().filter(|record| record.passes() && record.carries_an_alt()) {
        let planted = if record.is_indel() {
            record.applied_to(reference).is_some_and(|sequence| indels.contains(&sequence))
        } else {
            snvs.iter().any(|snv| {
                snv.pos + 1 == record.pos && record.alts.contains(&char::from(snv.alt).to_string())
            })
        };
        ensure!(planted, "nothing was planted at {}", record.line);
    }
    Ok(())
}

/// Nothing is called where nothing was planted, methylation included: a
/// methylated CpG is no `C>T` or `G>A`.
#[hegel::test(test_cases = pinned::cases(4))]
fn nothing_is_called_where_nothing_was_planted(tc: TestCase) {
    let plan = tc.draw(plans(true));
    let result = call(&plan).and_then(|records| check_no_false_calls(&plan, &records));
    assert!(result.is_ok(), "{result:?}");
}

/// The `M5mC` record at each side of each CpG both haplotypes carry.
fn methylation_records<'r>(
    plan: &Plan,
    records: &'r [Record],
) -> Result<Vec<(PlantedCpg, &'r Record)>> {
    let mut sides = Vec::new();
    for cpg in plan.cpgs() {
        if !matches!(cpg.carrier, None | Some(Carrier::Both)) {
            continue;
        }
        for pos in [cpg.c + 1, cpg.c + 2] {
            let record = records
                .iter()
                .find(|record| record.pos == pos && record.field("M5mC").is_some())
                .ok_or_else(|| eyre!("{cpg:?}: no M5mC at {pos}"))?;
            sides.push((cpg, record));
        }
    }
    Ok(sides)
}

/// β̂ at each side of a CpG both haplotypes carry is within five standard
/// deviations of the planted level, over the `DPM5mC` reads it counts. Those
/// are the reads of the strand that shows methylation there, the two
/// haplotypes alike, less any that end on the CpG (which β leaves out); so at
/// least half the planted depth.
fn check_methylation(plan: &Plan, records: &[Record]) -> Result<()> {
    let planted = 2 * plan.depth_per_hap_and_strand();
    let [one, two] = plan.methylation;
    let beta = f64::midpoint(one, two);
    for (cpg, record) in methylation_records(plan, records)? {
        let reads: usize = record.field("DPM5mC").unwrap_or(".").parse()?;
        ensure!(
            (planted / 2..=planted).contains(&reads),
            "{cpg:?}: DPM5mC {reads} of {planted} planted: {}",
            record.line
        );
        let variance = (one * (1.0 - one) + two * (1.0 - two)) * reads as f64 / 2.0;
        let tolerance = 5.0 * variance.sqrt() / reads as f64 + 1e-5;
        let observed: f64 = record.field("M5mC").unwrap_or(".").parse()?;
        ensure!(
            (observed - beta).abs() <= tolerance,
            "{cpg:?}: β̂ {observed}, planted {beta} ± {tolerance}: {}",
            record.line
        );
    }
    Ok(())
}

/// Planted methylation is recovered at reference CpGs and at de-novo CpGs on
/// both haplotypes.
#[hegel::test(test_cases = pinned::cases(4))]
fn planted_methylation_is_recovered(tc: TestCase) {
    let plan = tc.draw(plans(false));
    let result = call(&plan).and_then(|records| check_methylation(&plan, &records));
    assert!(result.is_ok(), "{result:?}");
}

/// `DPM5mC` counts the reads `M5mC_Strands` counts as unmodified or modified.
/// Smallest failure: a CpG where 4 of the 16 OT reads end on the `C`;
/// `M5mC_Strands=16,0,16,0`, `DPM5mC=12`.
#[hegel::test(test_cases = pinned::cases(4))]
#[hegel::explicit_test_case(plan = Plan {
    seed: 1,
    events: vec![(FIRST_EVENT + 24, Event::Cpg)],
    step: 25,
    copies: 1,
    methylation: [0.0, 0.0],
    indel_flanks: IndelFlanks::Neutral,
})]
fn methylation_depth_is_the_strand_count(tc: TestCase) {
    let plan = tc.draw(plans(false));
    let check = |records: Vec<Record>| -> Result<()> {
        for (cpg, record) in methylation_records(&plan, &records)? {
            let info = record.line.split('\t').nth(7).unwrap_or_default();
            let strands: Vec<usize> = info
                .split(';')
                .find_map(|field| field.strip_prefix("M5mC_Strands="))
                .ok_or_else(|| eyre!("no M5mC_Strands: {}", record.line))?
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let depth: usize = record.field("DPM5mC").unwrap_or(".").parse()?;
            let counted: usize = strands.iter().take(2).sum();
            ensure!(depth == counted, "{cpg:?}: DPM5mC {depth}, M5mC_Strands {strands:?}");
        }
        Ok(())
    };
    let result = call(&plan).and_then(check);
    assert!(result.is_ok(), "{result:?}");
}

/// An insertion of `C` before a reference `G`, on both haplotypes, makes a
/// CpG, 50 % methylated on haplotype one and fully on two, with a het SNV
/// 8 bp before it on haplotype one. The OB reads show `A` at the `G` (11 of
/// 16) and the OT reads none (0 of 16): methylation, since a `G>A` would show
/// on OT as well. The `G` is a de-novo CpG's, its β (¾) from the OB reads.
/// Before the insertion made a CpG, the SNV engine called the `G` a PASS
/// `G>A`, phased with the SNV.
#[test]
fn an_insertion_that_makes_a_cpg_is_no_snv() -> Result<()> {
    let plan = Plan {
        seed: 5,
        events: vec![
            (FIRST_EVENT + 80, Event::Snv { alt: 0, carrier: Carrier::One }),
            (FIRST_EVENT + 88, Event::Insertion { bases: b"C".to_vec(), carrier: Carrier::Both }),
        ],
        step: 25,
        copies: 1,
        methylation: [0.5, 1.0],
        indel_flanks: IndelFlanks::Random,
    };
    let reference = plan.reference();
    ensure!(reference.reference().get(FIRST_EVENT + 89) == Some(&b'G'), "no G after the insertion");
    let records = call(&plan)?;
    check_no_false_calls(&plan, &records)?;
    let g = records
        .iter()
        .find(|record| record.pos == FIRST_EVENT + 90 && !record.is_indel())
        .ok_or_else(|| eyre!("no record at the G"))?;
    let beta: f64 = g.field("M5mC").unwrap_or(".").parse().map_err(|_| eyre!("{}", g.line))?;
    ensure!(g.line.contains("CPGnovo") && (beta - 0.75).abs() < 0.25, "{}", g.line);
    Ok(())
}

/// An insertion of `G` after a reference `C`, on both haplotypes, makes the
/// `C` a de-novo CpG's, its β (¾) from the OT reads, whose next base is the
/// inserted `G`.
#[test]
fn an_insertion_after_a_c_makes_a_cpg_of_it() -> Result<()> {
    let insertion = |pos| Plan {
        seed: 5,
        events: vec![(pos, Event::Insertion { bases: b"G".to_vec(), carrier: Carrier::Both })],
        step: 25,
        copies: 1,
        methylation: [0.5, 1.0],
        indel_flanks: IndelFlanks::Random,
    };
    let reference = insertion(FIRST_EVENT).reference();
    let base = |pos: usize| reference.reference().get(pos).copied();
    let anchor = (FIRST_EVENT..FIRST_EVENT + 100)
        .find(|&pos| base(pos) == Some(b'C') && base(pos + 1) != Some(b'G'))
        .ok_or_else(|| eyre!("no C to anchor on"))?;
    let plan = insertion(anchor);
    let records = call(&plan)?;
    check_no_false_calls(&plan, &records)?;
    let c = records
        .iter()
        .find(|record| record.pos == anchor + 1 && !record.is_indel())
        .ok_or_else(|| eyre!("no record at the C"))?;
    let beta: f64 = c.field("M5mC").unwrap_or(".").parse().map_err(|_| eyre!("{}", c.line))?;
    ensure!(c.line.contains("CPGnovo") && (beta - 0.75).abs() < 0.25, "{}", c.line);
    Ok(())
}
