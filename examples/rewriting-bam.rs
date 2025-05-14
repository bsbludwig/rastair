use color_eyre::eyre::{Context, ContextCompat, Result};
use rastair2::utils::{Base, MethylatedPositions};
use rust_htslib::bam::{self, Header, Read, Record, Writer, header::HeaderRecord};
use smallvec::SmallVec;
use tracing_subscriber::layer::SubscriberExt as _;

fn main() -> Result<()> {
    color_eyre::install()?;
    let subscriber =
        tracing_subscriber::Registry::default().with(tracing_error::ErrorLayer::default()).with(
            tracing_subscriber::fmt::Layer::default()
                .with_target(true)
                .with_writer(std::io::stderr),
        );
    tracing::subscriber::set_global_default(subscriber)?;

    let bam_file = "./tests/data/test.bam";
    let mut bam = bam::IndexedReader::from_path(bam_file).wrap_err("failed to open bam file")?;
    bam.fetch("chr19:6105400-6105410").wrap_err("error fetching range")?;

    rewrite_bam(&mut bam, "enhanced.bam").wrap_err("rewrite")?;

    Ok(())
}

#[tracing::instrument(skip_all)]
fn rewrite_bam(bam: &mut bam::IndexedReader, output_file: &str) -> Result<()> {
    let header = {
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"HD").push_tag(b"VN", "1.0").push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ").push_tag(b"SN", "chr19").push_tag(b"LN", 6105400),
        );
        header
    };
    let mut writer = Writer::from_path(output_file, &header, bam::Format::Bam)
        .wrap_err("failed to create writer")?;

    let mut record = Record::new();
    bam.read(&mut record).wrap_err("failed to read record")??;

    record.set_mapq(67); // marker
    // record.push_aux(b"Mm", Aux::String("C+m,1;C-m,1;")).wrap_err("failed to push aux")?; // this is not validated at all

    let mods = MethylatedPositions { base: Base::C, positions: SmallVec::from([6, 19, 20]) };
    mods.apply_to_record(&mut record)?;

    writer.write(&record).wrap_err("failed to write record")?;
    drop(writer);

    Ok(())
}
