mod debug;
pub mod utils;

mod vcf_tests;

#[cfg(test)]
mod figure_out_outputs_tests {
    use crate::{CallParams, call::record_filters::RecordFilters, utils::default};
    use clio::ClioPath;
    use color_eyre::Result;

    fn test_params() -> CallParams {
        CallParams {
            segments: crate::sequence::ReaderParams {
                bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
                fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
                region: None,
            },
            segmentation: default(),
            variant_calling: default(),
            denovo_cpg: default(),
            methylation: default(),
            ml: default(),
            vcf: default(),
            bed: default(),
            record_filters: RecordFilters { vcf_all: false, cpgs_only: false },
            total_threads: 2,
        }
    }

    #[test]
    fn test_default_vcf_output() -> Result<()> {
        let mut params = test_params();
        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_none());
        assert!(!params.record_filters.cpgs_only);

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_some());
        assert!(params.bed.bed.is_none());
        assert!(!params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_default_bed_output_with_cpgs_only() -> Result<()> {
        let mut params = test_params();
        params.record_filters.cpgs_only = true;

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_some());
        assert!(params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_vcf_output_specified() -> Result<()> {
        let mut params = test_params();
        params.vcf.vcf = Some(ClioPath::new("output.vcf").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_some());
        assert!(params.bed.bed.is_none());
        assert!(!params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_bed_output_specified() -> Result<()> {
        let mut params = test_params();
        params.bed.bed = Some(ClioPath::new("output.bed").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_some());
        assert!(params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_both_outputs_different_files() -> Result<()> {
        let mut params = test_params();
        params.vcf.vcf = Some(ClioPath::new("output.vcf").unwrap());
        params.bed.bed = Some(ClioPath::new("output.bed").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_some());
        assert!(params.bed.bed.is_some());

        Ok(())
    }

    #[test]
    fn test_both_outputs_same_file_errors() {
        let mut params = test_params();
        let same_path = ClioPath::new("output.txt").unwrap();
        params.vcf.vcf = Some(same_path.clone());
        params.bed.bed = Some(same_path);

        let result = params.figure_out_outputs();

        assert!(result.is_err());
        let err_msg = format!("{:#}", result.unwrap_err());
        assert!(err_msg.contains("Can't write both VCF and BED output to the same file"));
    }

    #[test]
    fn test_vcf_with_bed_extension_switches_to_bed() -> Result<()> {
        let mut params = test_params();
        params.vcf.vcf = Some(ClioPath::new("output.bed").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_some());
        assert!(params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_vcf_with_bed_gz_extension_switches() -> Result<()> {
        let mut params = test_params();
        params.vcf.vcf = Some(ClioPath::new("output.bed.gz").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_some());
        assert!(params.record_filters.cpgs_only);

        Ok(())
    }

    #[test]
    fn test_bed_output_sets_cpgs_only() -> Result<()> {
        let mut params = test_params();
        params.bed.bed = Some(ClioPath::new("output.bed").unwrap());
        assert!(!params.record_filters.cpgs_only);

        params.figure_out_outputs()?;

        assert!(params.record_filters.cpgs_only);
        assert!(params.vcf.vcf.is_none());
        assert!(params.bed.bed.is_some());

        Ok(())
    }

    #[test]
    fn test_vcf_and_bed_both_specified_no_switch() -> Result<()> {
        let mut params = test_params();
        params.vcf.vcf = Some(ClioPath::new("output.bed").unwrap());
        params.bed.bed = Some(ClioPath::new("methylation.bed").unwrap());

        params.figure_out_outputs()?;

        assert!(params.vcf.vcf.is_some());
        assert!(params.bed.bed.is_some());

        Ok(())
    }
}
