use crate::{
    CallParams,
    call::{call, record_filters::RecordFilters},
    utils::default,
};
use clio::ClioPath;
use color_eyre::Result;

#[test]
// Un-ignore, and press "debug this test" in your IDE. Better than trying to
// deal with a debug config that first has do a build and then pass the right
// args :)
#[ignore = "used for debugging"]
fn call_default_bed() -> Result<()> {
    let _test_writer =
        tracing_subscriber::fmt().with_test_writer().with_env_filter("rastair=debug").try_init();

    call(CallParams {
        segments: crate::sequence::ReaderParams {
            bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
            fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
            region: Some("chr19:6103000-6103100".parse()?),
        },
        record_filters: RecordFilters { vcf_all: false, cpgs_only: true },
        segmentation: default(),
        require_tags: default(),
        variant_calling: default(),
        indel: default(),
        denovo_cpg: default(),
        methylation: default(),
        ml: default(),
        vcf: default(),
        bed: default(),
        total_threads: 2,
    })
}
