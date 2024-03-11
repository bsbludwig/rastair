use log::warn;
use rust_htslib::bam::Record;
use bio::bio_types::sequence::SequenceReadPairOrientation::{F1R2, F2R1, R1F2, R2F1, self};

pub fn read_pair_orientation(record: &Record, exclude_ambiguous: bool) -> SequenceReadPairOrientation
{
    let mut read_pair_orientation = record.read_pair_orientation();
    if ! exclude_ambiguous
    {
        read_pair_orientation = match read_pair_orientation
        {
            F1R2 | R2F1 => F1R2,
            F2R1 | R1F2 => F2R1,
            SequenceReadPairOrientation::None => {
                warn!("Orientation of {} cannot be unambiguously determined", String::from_utf8(Vec::from(record.qname())).unwrap_or_default());

                if record.is_first_in_template() && record.is_mate_reverse() ||
                record.is_last_in_template() && record.is_reverse()
                {
                    F1R2
                }
                // F2R1
                else if record.is_first_in_template() && record.is_reverse() ||
                        record.is_last_in_template() && record.is_mate_reverse()
                {
                    F2R1
                }
                else {
                    SequenceReadPairOrientation::None
                }
            },
            _   =>  SequenceReadPairOrientation::None // This should be impossible?
        };
    }
    read_pair_orientation
}