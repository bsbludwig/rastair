//! Feature calculation implementations for ML models
//!
//! Trait-based abstraction for calculating features from variant metrics.
//! Different implementations can be swapped to support various feature sets and
//! model variants.

use super::types::MlFeatureSet;
use crate::metrics::{MetricsForAlt, MetricsForIndel, PileupMetrics};
use color_eyre::{Result, eyre::Context as _};
use ndarray::Array2;
use std::fmt;

pub mod shared;
pub mod standard;
pub mod utils;

/// Define a feature struct whose memory layout *is* its ML feature vector.
///
/// Every field is `f64` or `[f64; N]`, so a `#[repr(C)]` struct has no padding
/// and is [`bytemuck::Pod`]. The field declaration order is the feature order,
/// which removes the hand-counted `buf[start..end]` index arithmetic the old
/// code relied on. `as_row()` reinterprets the struct as `&[f64]` with no copy.
///
/// Names are derived from the same declaration, so they cannot drift from the
/// values: a `scalar` field yields one name (its identifier), an `array` field
/// yields its explicit per-slot names, and a `flatten`ed field delegates to the
/// nested type's names. `names()` is a cold path (TSV headers, importance
/// export), so it allocates a `Vec` rather than fighting const concatenation.
///
/// ```ignore
/// define_features! {
///     pub struct InsertionFeatures {
///         flatten common: CommonIndelFeatures;
///         /// RMS base quality of the inserted bases.
///         scalar insertion_baseq_rms;
///     }
/// }
/// ```
macro_rules! define_features {
    (
        $(#[$meta:meta])*
        $vis:vis struct $name:ident { $($body:tt)* }
    ) => {
        define_features!(@munch
            meta { $(#[$meta])* } vis { $vis } name { $name }
            fields { } names { }
            rest { $($body)* }
        );
    };

    // scalar f64 field
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($nm:tt)* }
        rest { $(#[$attr:meta])* scalar $field:ident ; $($rest:tt)* }
    ) => {
        define_features!(@munch meta {$($m)*} vis {$v} name {$n}
            fields { $($f)* $(#[$attr])* pub $field: f64, }
            names { $($nm)* (strs stringify!($field)) }
            rest { $($rest)* }
        );
    };

    // [f64; N] field with one explicit name per slot
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($nm:tt)* }
        rest {
            $(#[$attr:meta])* array $field:ident : $len:literal = [ $($label:literal),+ $(,)? ] ;
            $($rest:tt)*
        }
    ) => {
        define_features!(@munch meta {$($m)*} vis {$v} name {$n}
            fields { $($f)* $(#[$attr])* pub $field: [f64; $len], }
            names { $($nm)* (strs $($label),+) }
            rest { $($rest)* }
        );
    };

    // flattened nested feature struct
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($nm:tt)* }
        rest { $(#[$attr:meta])* flatten $field:ident : $ty:ty ; $($rest:tt)* }
    ) => {
        define_features!(@munch meta {$($m)*} vis {$v} name {$n}
            fields { $($f)* $(#[$attr])* pub $field: $ty, }
            names { $($nm)* (flat $ty) }
            rest { $($rest)* }
        );
    };

    // emit one accumulated name descriptor into `$out` (kept in the same
    // expansion as the `extend_names` binding so hygiene unifies)
    (@emit $out:ident (strs $($name:expr),+)) => { $( $out.push($name); )+ };
    (@emit $out:ident (flat $ty:ty)) => { <$ty>::extend_names($out); };

    // done: emit the struct + impl
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($grp:tt)* } rest { }
    ) => {
        $($m)*
        #[repr(C)]
        #[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
        $v struct $n { $($f)* }

        impl $n {
            /// Number of `f64` features in this struct's flat layout.
            pub const FEATURES: usize =
                ::core::mem::size_of::<Self>() / ::core::mem::size_of::<f64>();

            /// Append this struct's feature names, in layout order, to `out`.
            pub fn extend_names(out: &mut ::std::vec::Vec<&'static str>) {
                $( define_features!(@emit out $grp); )*
            }

            /// Feature names in layout order; one entry per `f64` slot.
            pub fn names() -> ::std::vec::Vec<&'static str> {
                let mut out = ::std::vec::Vec::with_capacity(Self::FEATURES);
                Self::extend_names(&mut out);
                out
            }

            /// Reinterpret the struct as its flat feature row, with no copy.
            #[inline]
            pub fn as_row(&self) -> &[f64] {
                bytemuck::cast_slice(::core::slice::from_ref(self))
            }
        }
    };
}

pub(crate) use define_features;

#[derive(Debug, Clone, Copy)]
pub struct FeatureNum {
    pub cpg: usize,
    pub denovo_cpg: usize,
    pub others: usize,
    pub insertion: usize,
    pub deletion: usize,
}

/// Feature names per model, in the same layout order as the feature vectors.
#[derive(Debug, Clone)]
pub struct FeatureNames {
    pub cpg: Vec<&'static str>,
    pub denovo_cpg: Vec<&'static str>,
    pub others: Vec<&'static str>,
    pub insertion: Vec<&'static str>,
    pub deletion: Vec<&'static str>,
}

pub type FeatureCalculatorBox = Box<dyn FeatureCalculator>;

/// Calculate ML features from variant metrics
pub trait FeatureCalculator: fmt::Debug + Send + Sync {
    fn feature_num(&self) -> FeatureNum;

    /// Feature names per model, aligned with the feature vectors' layout order.
    fn feature_names(&self) -> FeatureNames;

    /// Calculate features for a CpG methylation candidate
    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;

    /// Calculate features for a denovo CpG candidate
    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;

    /// Calculate features for other variants
    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;

    /// Calculate features for insertion
    fn calculate_insertion(&self, current: &MetricsForIndel) -> Result<Array2<f64>>;

    /// Calculate features for deletion
    fn calculate_deletion(&self, current: &MetricsForIndel) -> Result<Array2<f64>>;
}

impl MlFeatureSet {
    pub fn get_calculator(&self) -> FeatureCalculatorBox {
        match self {
            MlFeatureSet::Standard => Box::new(StandardFeatures),
            MlFeatureSet::Simple => Box::new(SimpleFeatures),
        }
    }
}

/// Standard implementation of feature calculation using all features
#[derive(Debug, Clone, Copy)]
pub struct StandardFeatures;

/// Wrap a flat feature row into a `(1, N)` `Array2` for the forest predictor.
fn row_to_array(row: &[f64]) -> Result<Array2<f64>> {
    Array2::from_shape_vec((1, row.len()), row.to_vec())
        .wrap_err("Failed to build feature array from row")
}

impl FeatureCalculator for StandardFeatures {
    fn feature_num(&self) -> FeatureNum {
        FeatureNum {
            cpg: standard::CpgFeatures::FEATURES,
            denovo_cpg: standard::DenovoCpgFeatures::FEATURES,
            others: standard::OthersFeatures::FEATURES,
            insertion: standard::InsertionFeatures::FEATURES,
            deletion: standard::DeletionFeatures::FEATURES,
        }
    }

    fn feature_names(&self) -> FeatureNames {
        FeatureNames {
            cpg: standard::CpgFeatures::names(),
            denovo_cpg: standard::DenovoCpgFeatures::names(),
            others: standard::OthersFeatures::names(),
            insertion: standard::InsertionFeatures::names(),
            deletion: standard::DeletionFeatures::names(),
        }
    }

    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        row_to_array(standard::CpgFeatures::extract(current, before, after)?.as_row())
    }

    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        row_to_array(standard::DenovoCpgFeatures::extract(current, before, after)?.as_row())
    }

    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        row_to_array(standard::OthersFeatures::extract(current, before, after)?.as_row())
    }

    fn calculate_insertion(&self, current: &MetricsForIndel) -> Result<Array2<f64>> {
        row_to_array(standard::InsertionFeatures::extract(current).as_row())
    }

    fn calculate_deletion(&self, current: &MetricsForIndel) -> Result<Array2<f64>> {
        row_to_array(standard::DeletionFeatures::extract(current).as_row())
    }
}

/// Very basic feature calculation using small subset of features
#[derive(Debug, Clone, Copy)]
pub struct SimpleFeatures;

impl SimpleFeatures {
    /// `SimpleFeatures` uses all of [`CommonFeatures`] except `region_entropy`.
    const FEATURES: usize = shared::CommonSectionA::FEATURES - 1 + shared::CommonSectionB::FEATURES;

    fn calculate_basic(&self, current: &MetricsForAlt) -> Result<Array2<f64>> {
        let common = shared::CommonFeatures::extract(current);
        let mut features = Vec::with_capacity(Self::FEATURES);
        features.extend_from_slice(&common.base_encoding);
        features.extend_from_slice(&common.position_metrics);
        features.extend_from_slice(&common.context_encoding);
        features.extend_from_slice(&common.depth_ratios);
        features.extend_from_slice(&common.base_quality_metrics);
        features.extend_from_slice(&common.mapping_quality_metrics);
        features.extend_from_slice(&common.read_metrics);
        Array2::from_shape_vec((1, features.len()), features)
            .wrap_err("Failed to create basic feature array")
    }

    /// Names mirroring [`calculate_basic`](Self::calculate_basic): all of
    /// [`CommonFeatures`] except `region_entropy`.
    fn basic_names() -> Vec<&'static str> {
        shared::CommonSectionA::names()
            .into_iter()
            .filter(|name| *name != "region_entropy")
            .chain(shared::CommonSectionB::names())
            .collect()
    }
}

impl FeatureCalculator for SimpleFeatures {
    fn feature_num(&self) -> FeatureNum {
        let n = Self::FEATURES;
        FeatureNum { cpg: n, denovo_cpg: n, others: n, insertion: 0, deletion: 0 }
    }

    fn feature_names(&self) -> FeatureNames {
        FeatureNames {
            cpg: Self::basic_names(),
            denovo_cpg: Self::basic_names(),
            others: Self::basic_names(),
            insertion: Vec::new(),
            deletion: Vec::new(),
        }
    }

    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }

    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }

    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }

    fn calculate_insertion(&self, _current: &MetricsForIndel) -> Result<Array2<f64>> {
        todo!()
    }

    fn calculate_deletion(&self, _current: &MetricsForIndel) -> Result<Array2<f64>> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::shared::{CommonSectionA, CommonSectionB};
    use super::standard::indel::CommonIndelFeatures;
    use super::standard::{
        CpgFeatures, DeletionFeatures, DenovoCpgFeatures, InsertionFeatures, OthersFeatures,
    };

    /// One name per `f64` slot, for every feature struct. A drift here means
    /// the `names()` output no longer aligns with the feature vectors.
    #[test]
    fn names_align_with_feature_count() {
        assert_eq!(CommonSectionA::names().len(), CommonSectionA::FEATURES);
        assert_eq!(CommonSectionB::names().len(), CommonSectionB::FEATURES);
        assert_eq!(CommonIndelFeatures::names().len(), CommonIndelFeatures::FEATURES);
        assert_eq!(CpgFeatures::names().len(), CpgFeatures::FEATURES);
        assert_eq!(DenovoCpgFeatures::names().len(), DenovoCpgFeatures::FEATURES);
        assert_eq!(OthersFeatures::names().len(), OthersFeatures::FEATURES);
        assert_eq!(InsertionFeatures::names().len(), InsertionFeatures::FEATURES);
        assert_eq!(DeletionFeatures::names().len(), DeletionFeatures::FEATURES);
    }

    /// The feature vector layout is frozen by every trained model. These counts
    /// are the historical layout sizes; changing one silently invalidates the
    /// corresponding model, so they are pinned here.
    #[test]
    fn feature_counts_are_stable() {
        assert_eq!(CommonSectionA::FEATURES, 33);
        assert_eq!(CommonSectionB::FEATURES, 18);
        assert_eq!(CpgFeatures::FEATURES, 55);
        assert_eq!(DenovoCpgFeatures::FEATURES, 56);
        assert_eq!(OthersFeatures::FEATURES, 54);
        assert_eq!(CommonIndelFeatures::FEATURES, 33);
        assert_eq!(InsertionFeatures::FEATURES, 34);
        assert_eq!(DeletionFeatures::FEATURES, 38);
    }

    /// Snapshot the full name→index layout of every model. This is the
    /// human-checkable replacement for the old "never change the order"
    /// comments: reordering a field changes the snapshot and fails the test.
    #[test]
    fn feature_name_layout() {
        let render = |name: &str, names: Vec<&str>| {
            let rows: String =
                names.iter().enumerate().map(|(i, n)| format!("  {i:>2}  {n}\n")).collect();
            format!("{name} ({} features)\n{rows}", names.len())
        };

        let report = [
            render("cpg", CpgFeatures::names()),
            render("denovo_cpg", DenovoCpgFeatures::names()),
            render("others", OthersFeatures::names()),
            render("insertion", InsertionFeatures::names()),
            render("deletion", DeletionFeatures::names()),
        ]
        .join("\n");

        insta::assert_snapshot!(report);
    }
}
