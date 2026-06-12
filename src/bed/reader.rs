use color_eyre::{
    Result, Section as _,
    eyre::{Context as _, ContextCompat as _, bail},
};
use noodles::{core::Region, tabix};
use seqair_types::{Base, SmallVec, SmolStr};
use std::{fs::File, path::Path};
use tracing::{instrument, warn};

type TabixReader = noodles::csi::io::IndexedReader<
    noodles::bgzf::io::Reader<File>,
    noodles::csi::binning_index::Index<Vec<noodles::bgzf::VirtualPosition>>,
>;

pub struct RastairBedReader {
    reader: TabixReader,
}

impl RastairBedReader {
    #[instrument(level = "debug")]
    pub fn new(bed_path: &Path) -> Result<Self> {
        let reader = tabix::io::indexed_reader::Builder::default()
            .build_from_path(bed_path)
            .wrap_err("Failed to open tabix file");

        let reader = if "-" == bed_path.to_string_lossy() {
            reader.note("Rastair can only read calls from local .bed.gz files as it also needs to load the tabix index")?
        } else if bed_path.ends_with("gz") {
            reader.with_suggestion(|| {
                let path = bed_path.display();
                format!(
                    "You can create a tabix index with:\n\
                     \n\
                     tabix {path}"
                )
            })?
        } else {
            reader.with_suggestion(|| {
                let path = bed_path.display();
                format!(
                    "You can compress and index the calls file with:\n\
                     \n\
                     bgzip {path}\n\
                     tabix {path}.gz"
                )
            })?
        };

        // Warn if the tabix index is older than the BED file
        let tbi_path = {
            let mut p = bed_path.as_os_str().to_owned();
            p.push(".tbi");
            std::path::PathBuf::from(p)
        };
        if let (Ok(bed_meta), Ok(tbi_meta)) =
            (std::fs::metadata(bed_path), std::fs::metadata(&tbi_path))
            && let (Ok(bed_modified), Ok(tbi_modified)) = (bed_meta.modified(), tbi_meta.modified())
            && bed_modified > tbi_modified
        {
            let path = bed_path.display();
            warn!(
                "Tabix index for {path} appears to be older than the BED file. \
                         This can cause query errors. Regenerate it with: tabix {path}"
            );
        }

        Ok(RastairBedReader { reader })
    }

    #[instrument(level = "debug", skip(self))]
    pub fn query(&mut self, region: &Region) -> Result<Vec<SimpleRastairBedRecord>> {
        let query = match self.reader.query(region) {
            Ok(query) => query,
            Err(e) if e.kind() == std::io::ErrorKind::InvalidInput => {
                // Region's chromosome not in the tabix index — no calls for this region
                return Ok(Vec::new());
            }
            Err(e) => return Err(e).wrap_err("Failed to query tabix file"),
        };
        let mut records = Vec::new();
        for result in query {
            let record =
                result.wrap_err("Failed to read record from tabix file").and_then(|record| {
                    SimpleRastairBedRecord::from_bed_line(record.as_ref())
                        .wrap_err("Failed to parse rastair BED record")
                });

            match record {
                Ok(record) => {
                    records.push(record);
                }
                Err(error) => {
                    warn!(%error, "Failed to read record");
                }
            }
        }
        Ok(records)
    }
}

#[derive(Debug, Clone)]
#[expect(unused, reason = "for completeness")]
pub struct SimpleRastairBedRecord {
    pub chrom: SmolStr,
    pub pos: u32,
    pub call: RastairCall,
}

#[derive(Debug, Clone)]
#[allow(unused, reason = "we list all fields for clarity")]
struct RefRastairBedRecord<'src> {
    chr: &'src str,
    start: &'src str,
    end: &'src str,
    name: &'src str,
    beta: &'src str,
    strand: &'src str,
    unmod: &'src str,
    r#mod: &'src str,
    no_snp: &'src str,
    snp: &'src str,
    coverage: &'src str,
    genotype: &'src str,
    genotype_likelihood: &'src str,
    genotype_confidence: &'src str,
    de_novo: &'src str,
}

#[derive(Debug, Clone)]
#[expect(unused, reason = "for completeness")]
pub enum RastairCall {
    /// A CpG site on the reference genome
    Cpg {
        /// The base on the reference strand (C or G)
        base: Base,
        /// Whether the CpG is methylated
        methylated: bool,
    },
    /// A de novo CpG site not present in the reference genome
    DeNovoCpg {
        /// The base on the reference strand (C or G)
        base: Base,
        /// Whether the CpG is methylated
        methylated: bool,
    },
    /// A SNP (single nucleotide polymorphism)
    Snp {
        /// Reference base
        from: Base,
        /// Variant base
        to: Base,
    },
}

impl SimpleRastairBedRecord {
    fn from_bed_line(record: &str) -> Result<SimpleRastairBedRecord> {
        let mut columns = record.split('\t');
        let rec = RefRastairBedRecord {
            chr: columns.next().wrap_err("Missing chr column")?,
            start: columns.next().wrap_err("Missing start column")?,
            end: columns.next().wrap_err("Missing end column")?,
            name: columns.next().wrap_err("Missing name column")?,
            beta: columns.next().wrap_err("Missing beta column")?,
            strand: columns.next().wrap_err("Missing strand column")?,
            unmod: columns.next().wrap_err("Missing unmod column")?,
            r#mod: columns.next().wrap_err("Missing mod column")?,
            no_snp: columns.next().wrap_err("Missing no_snp column")?,
            snp: columns.next().wrap_err("Missing snp column")?,
            coverage: columns.next().wrap_err("Missing coverage column")?,
            genotype: columns.next().wrap_err("Missing genotype column")?,
            genotype_likelihood: columns.next().wrap_err("Missing genotype_likelihood column")?,
            genotype_confidence: columns.next().wrap_err("Missing genotype_confidence column")?,
            de_novo: columns.next().wrap_err("Missing de_novo column")?,
        };

        let genotype = rec.genotype.split('/').collect::<SmallVec<&str, 2>>();
        let [from, to] = genotype.as_slice() else {
            bail!(
                "Failed to parse genotype column, expected format 'N/N' or similar, but got {genotype:?}"
            );
        };
        let from = from.parse().wrap_err("Failed to parse from base")?;
        let to = to.parse().wrap_err("Failed to parse to base")?;
        let beta = rec.r#beta.parse::<f32>().wrap_err("Failed to parse beta value")?;

        let call = if beta == 0.0 && rec.no_snp != "0" {
            RastairCall::Snp { from, to }
        } else if rec.de_novo == "NEW" {
            RastairCall::DeNovoCpg { base: from, methylated: beta > 0.5 }
        } else {
            RastairCall::Cpg { base: from, methylated: beta > 0.5 }
        };

        Ok(SimpleRastairBedRecord {
            chrom: SmolStr::from(rec.chr),
            pos: rec.start.parse().wrap_err("Failed to parse position")?,
            call,
        })
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;
    use std::io::Write;

    use super::*;

    /// Creates a small bgzf-compressed BED file with records only on chr19,
    /// runs `tabix` to index it, and returns the path to the .bed.gz file.
    fn create_test_bed_gz(dir: &std::path::Path) -> Result<std::path::PathBuf> {
        let bed_path = dir.join("test.bed.gz");
        let compression_level =
            bgzf::CompressionLevel::try_from(6).map_err(|e| color_eyre::eyre::eyre!("{e}"))?;
        let mut writer =
            bgzf::Writer::new(std::io::BufWriter::new(File::create(&bed_path)?), compression_level);
        writeln!(
            writer,
            "#chr\tstart\tend\tname\tbeta_est\tstrand\tunmod\tmod\tno_snp\tsnp\tcoverage\tgenotype\tgt_p_score\tgt_conf_score\tcpg"
        )?;
        writeln!(
            writer,
            "chr19\t6107663\t6107664\t.\t0.05\t+\t19\t1\t20\t0\t41\tC/C\t99\t60\tNEW"
        )?;
        writeln!(writer, "chr19\t6107750\t6107751\t.\t0.73\t+\t3\t8\t7\t0\t18\tC/C\t99\t21\tREF")?;
        writer.finish().map_err(|e| color_eyre::eyre::eyre!("Failed to finish bgzf: {e}"))?;

        let output = std::process::Command::new("tabix")
            .args(["-p", "bed"])
            .arg(&bed_path)
            .output()
            .wrap_err("Failed to run tabix")?;
        color_eyre::eyre::ensure!(
            output.status.success(),
            "tabix failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        Ok(bed_path)
    }

    #[test]
    fn query_missing_chromosome_returns_empty() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = create_test_bed_gz(dir.path())?;

        let mut reader = RastairBedReader::new(&bed_path)?;

        // chr19 exists in the file
        let region: Region = "chr19:6107663-6107751".parse()?;
        let records = reader.query(&region)?;
        assert_eq!(records.len(), 2);

        // chrX does NOT exist in the tabix index — should return empty, not error
        let region: Region = "chrX:1-1000".parse()?;
        let records = reader.query(&region)?;
        assert!(records.is_empty());

        Ok(())
    }

    #[test]
    fn query_region_with_no_records_returns_empty() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = create_test_bed_gz(dir.path())?;

        let mut reader = RastairBedReader::new(&bed_path)?;

        // chr19 exists but this region has no records
        let region: Region = "chr19:1-100".parse()?;
        let records = reader.query(&region)?;
        assert!(records.is_empty());

        Ok(())
    }

    #[test]
    fn parse_rastair_bed_file() -> Result<()> {
        let example_bed = "\
            #chr	start	end	name	beta_est	strand	unmod	mod	no_snp	snp	coverage	genotype	gt_p_score	gt_conf_score	cpg
            chr19	6107663	6107664	.	0.05	+	19	1	20	0	41	C/C	99	60	NEW
            chr19	6107750	6107751	.	0.73	+	3	8	7	0	18	C/C	99	21	REF
            chr19	6107751	6107752	.	0.71	-	2	5	11	0	18	G/G	99	33	REF
        ";

        let parsed = example_bed
            .lines()
            .map(|l| l.trim())
            .enumerate()
            .filter(|(_i, x)| !x.is_empty())
            .filter(|(_i, x)| !x.starts_with('#'))
            .map(|(i, l)| {
                SimpleRastairBedRecord::from_bed_line(l)
                    .wrap_err_with(|| format!("Failed to parse line {i}: {l}"))
            })
            .collect::<Result<Vec<_>>>()?;

        assert_debug_snapshot!(parsed, @r#"
        [
            SimpleRastairBedRecord {
                chrom: "chr19",
                pos: 6107663,
                call: DeNovoCpg {
                    base: C,
                    methylated: false,
                },
            },
            SimpleRastairBedRecord {
                chrom: "chr19",
                pos: 6107750,
                call: Cpg {
                    base: C,
                    methylated: true,
                },
            },
            SimpleRastairBedRecord {
                chrom: "chr19",
                pos: 6107751,
                call: Cpg {
                    base: G,
                    methylated: true,
                },
            },
        ]
        "#);

        Ok(())
    }
}
