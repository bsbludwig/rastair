use crate::{
    sequence::Segment,
    utils::{Base, logging::ThisIsABug as _},
};
use better_default::Default;
use color_eyre::eyre::{Context as _, ContextCompat as _, Result};
use rastair_vcf::{HeaderField, InfoField, InfoFieldNumber, VcfField};
use rust_htslib::bcf::Record;
use smol_str::{SmolStr, SmolStrBuilder};

/// 5-base sequence context centered on the variant position
///
/// Printed in VCF as string with up to 5 characters.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SequenceContext {
    pub before_2: Option<Base>,
    pub before_1: Option<Base>,
    #[default(Base::Unknown)]
    pub me: Base,
    pub after_1: Option<Base>,
    pub after_2: Option<Base>,
}

impl SequenceContext {
    pub fn new(idx: usize, segment: &Segment) -> Result<Self> {
        const N: usize = 2;
        let me =
            segment.sequence.get(idx).wrap_err("Failed to get self base!").this_is_a_bug()?.into();
        let (before_2, before_1) =
            match segment.sequence_slice::<2>(idx.saturating_sub(N), idx)?.as_slice() {
                [b2, b1] => (Some(*b2), Some(*b1)),
                [b1] => (None, Some(*b1)),
                _ => (None, None),
            };
        let (after_1, after_2) = match segment.sequence_slice::<2>(idx + 1, idx + N + 1)?.as_slice()
        {
            [a1, a2] => (Some(*a1), Some(*a2)),
            [a1] => (None, Some(*a1)),
            _ => (None, None),
        };
        Ok(SequenceContext { before_2, before_1, me, after_1, after_2 })
    }

    fn to_smol_str(&self) -> SmolStr {
        let mut res = SmolStrBuilder::new();
        if let Some(base) = self.before_2 {
            res.push_str(base.into());
        };
        if let Some(base) = self.before_1 {
            res.push_str(base.into());
        };
        res.push_str(self.me.into());
        if let Some(base) = self.after_1 {
            res.push_str(base.into());
        };
        if let Some(base) = self.after_2 {
            res.push_str(base.into());
        };
        res.finish()
    }
}

impl VcfField for SequenceContext {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("SC5");
}

impl HeaderField for SequenceContext {
    const DESCRIPTION: &'static str = "5-base sequence context centered on the variant position";
}

impl InfoField for SequenceContext {
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Num(1);
    type Type = SmolStr;

    fn write(&self, record: &mut Record) -> Result<()> {
        record
            .push_info_string(Self::ID, &[self.to_smol_str().as_bytes()])
            .wrap_err("Failed to set field")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        sequence::{ChunkRegion, Region},
        utils::Base::*,
    };

    #[test]
    fn test_new_sequence_context() -> Result<()> {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr_test".into(), start: 100, end: 105 },
                last_position: 105,
            },
            sequence: b"ACGTA".to_vec(),
        };

        let context = SequenceContext::new(2, &segment)?;
        assert_eq!(
            SequenceContext {
                before_2: Some(A),
                before_1: Some(C),
                me: G,
                after_1: Some(T),
                after_2: Some(A)
            },
            context
        );
        Ok(())
    }

    #[test]
    fn test_new_sequence_context_at_start() -> Result<()> {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr_test".into(), start: 100, end: 105 },
                last_position: 105,
            },
            sequence: b"ACGTA".to_vec(),
        };

        let context = SequenceContext::new(0, &segment)?;
        assert_eq!(
            SequenceContext {
                before_2: None,
                before_1: None,
                me: A,
                after_1: Some(C),
                after_2: Some(G)
            },
            context
        );
        Ok(())
    }

    #[test]
    fn test_new_sequence_context_at_end() -> Result<()> {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr_test".into(), start: 100, end: 105 },
                last_position: 105,
            },
            sequence: b"ACGTA".to_vec(),
        };

        let context = SequenceContext::new(4, &segment)?;
        assert_eq!(
            SequenceContext {
                before_2: Some(G),
                before_1: Some(T),
                me: A,
                after_1: None,
                after_2: None,
            },
            context
        );
        Ok(())
    }

    #[test]
    fn test_new_sequence_context_oor() -> Result<()> {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr_test".into(), start: 100, end: 105 },
                last_position: 105,
            },
            sequence: b"ACGTA".to_vec(),
        };

        let res = SequenceContext::new(12, &segment);
        assert!(res.is_err());
        Ok(())
    }

    #[test]
    fn test_to_smol_str_complete_context() {
        let context = SequenceContext {
            before_2: Some(Base::A),
            before_1: Some(Base::T),
            me: Base::G,
            after_1: Some(Base::C),
            after_2: Some(Base::T),
        };
        assert_eq!(context.to_smol_str(), "ATGCT");
    }

    #[test]
    fn test_to_smol_str_partial_context() {
        let context = SequenceContext {
            before_2: Some(Base::C),
            before_1: Some(Base::C),
            me: Base::A,
            after_1: Some(Base::G),
            after_2: None,
        };
        assert_eq!(context.to_smol_str(), "CCAG");
    }

    #[test]
    fn test_to_smol_str_only_center() {
        let context = SequenceContext {
            before_2: None,
            before_1: None,
            me: Base::T,
            after_1: None,
            after_2: None,
        };
        assert_eq!(context.to_smol_str(), "T");
    }
}
