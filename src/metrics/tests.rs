use crate::{call::process, sequence::ReaderParams};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _},
};

#[test]
fn test_cpg_detection() -> Result<()> {
    let mut readers = ReaderParams::test_data().around("chr19", 6105711).pileup_readers()?;
    let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

    let pileup_mapping_params = process::PileupMappingParams::default();
    let (_segment, pileups) = process::get_pileups(&mut readers, &chunk, &pileup_mapping_params)
        .wrap_err("failed to process region")?;

    let pileups: Vec<_> = pileups.collect();

    // get pileups for a CpG site
    let ref_c =
        pileups.iter().find(|p| p.pos.as_u32() == 6105711).wrap_err("Could not find C pileup")?;
    let ref_g =
        pileups.iter().find(|p| p.pos.as_u32() == 6105712).wrap_err("Could not find G pileup")?;

    assert!(ref_c.reference_cpg().is_some());
    assert!(ref_g.reference_cpg().is_some());

    Ok(())
}
