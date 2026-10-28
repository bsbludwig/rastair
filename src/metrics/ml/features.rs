//! Feature vectors of the ML models.
//!
//! Each model's vector is a `#[repr(C)]` struct of `f32` / `[f32; N]` fields
//! built with `define_features!`: **field order is the feature order**. The
//! macro generates `FEATURES`, `names()`/`extend_names()` and a zero-copy
//! `as_row()`. Field kinds: `scalar name;`, `array name: N = ["..", ..];` and
//! `flatten name: Type;` (a nested feature struct). The shared features
//! are split into `CommonSectionA` and `CommonSectionB` because the alt-based models put
//! their own scalars (e.g. `alt_score`) between the two halves.
//!
//! **Every trained model freezes the order**: reordering a field silently
//! corrupts predictions. `feature_counts_are_stable` pins the counts and
//! `feature_name_layout` snapshots every `name→index` mapping; update the
//! snapshot only for an intended layout change, which needs a retrained
//! model. The same names head the `--feature-analytics` and
//! `--export-features` outputs of `ml train`.
//!
//! - **One column pipeline.** `call` and `ml train` both get their columns from
//!   `call::segment_columns`; a per-column quantity a feature reads belongs
//!   there, or training will not see it.
//!   `train::tests::training_and_calling_see_the_same_common_features` checks
//!   that the two agree.
//! - **Layout fingerprint.** `RastairFlatModel.feature_layout` holds an FNV-1a
//!   hash of each forest's `feature_names()`, one for both forest sets;
//!   `load_combined_model` refuses a model whose fingerprint differs from this
//!   build's, or a file of another format, and says to retrain. A layout
//!   change ships a retrained bundled model.

use super::types::MlModel;
use crate::metrics::{MetricsForAlt, PileupMetrics};
use color_eyre::{Result, eyre::Context as _};
use ndarray::Array2;

pub mod shared;
pub mod standard;
pub mod utils;

/// Define a feature struct whose memory layout *is* its ML feature vector.
///
/// Every field is `f32` or `[f32; N]`, so a `#[repr(C)]` struct has no padding
/// and is [`bytemuck::Pod`]. The field declaration order is the feature order,
/// so no index arithmetic places a feature. `as_row()` reinterprets the
/// struct as `&[f32]` with no copy, feeding the f32 inference forests
/// directly.
///
/// Names are derived from the same declaration, so they cannot drift from the
/// values: a `scalar` field yields one name (its identifier), an `array` field
/// yields its explicit per-slot names, and a `flatten`ed field delegates to the
/// nested type's names. `names()` is a cold path (TSV headers, importance
/// export), so it allocates a `Vec` rather than fighting const concatenation.
///
/// ```ignore
/// define_features! {
///     pub struct OthersFeatures {
///         flatten section_a: CommonSectionA;
///         /// Generic (strand-combined) alt score.
///         scalar alt_score;
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

    // scalar f32 field
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($nm:tt)* }
        rest { $(#[$attr:meta])* scalar $field:ident ; $($rest:tt)* }
    ) => {
        define_features!(@munch meta {$($m)*} vis {$v} name {$n}
            fields { $($f)* $(#[$attr])* pub $field: f32, }
            names { $($nm)* (strs stringify!($field)) }
            rest { $($rest)* }
        );
    };

    // [f32; N] field with one explicit name per slot
    (@munch meta {$($m:tt)*} vis {$v:vis} name {$n:ident}
        fields { $($f:tt)* } names { $($nm:tt)* }
        rest {
            $(#[$attr:meta])* array $field:ident : $len:literal = [ $($label:literal),+ $(,)? ] ;
            $($rest:tt)*
        }
    ) => {
        define_features!(@munch meta {$($m)*} vis {$v} name {$n}
            fields { $($f)* $(#[$attr])* pub $field: [f32; $len], }
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
            /// Number of `f32` features in this struct's flat layout.
            pub const FEATURES: usize =
                ::core::mem::size_of::<Self>() / ::core::mem::size_of::<f32>();

            /// Append this struct's feature names, in layout order, to `out`.
            pub fn extend_names(out: &mut ::std::vec::Vec<&'static str>) {
                $( define_features!(@emit out $grp); )*
            }

            /// Feature names in layout order; one entry per `f32` slot.
            pub fn names() -> ::std::vec::Vec<&'static str> {
                let mut out = ::std::vec::Vec::with_capacity(Self::FEATURES);
                Self::extend_names(&mut out);
                out
            }

            /// Reinterpret the struct as its flat feature row, with no copy.
            #[inline]
            pub fn as_row(&self) -> &[f32] {
                bytemuck::cast_slice(::core::slice::from_ref(self))
            }
        }
    };
}

pub(crate) use define_features;

/// How many features `model` reads.
pub fn feature_count(model: MlModel) -> usize {
    match model {
        MlModel::Others => standard::OthersFeatures::FEATURES,
        MlModel::Cpg => standard::CpgFeatures::FEATURES,
        MlModel::DenovoCpg => standard::DenovoCpgFeatures::FEATURES,
    }
}

/// `model`'s feature names, in the layout order of its feature vector.
pub fn feature_names(model: MlModel) -> Vec<&'static str> {
    match model {
        MlModel::Others => standard::OthersFeatures::names(),
        MlModel::Cpg => standard::CpgFeatures::names(),
        MlModel::DenovoCpg => standard::DenovoCpgFeatures::names(),
    }
}

/// The model that scores `current`, and its feature row for that model.
pub fn features(
    current: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> (MlModel, Result<Array2<f32>>) {
    if current.is_conversion() {
        let row = standard::CpgFeatures::extract(current, before, after);
        (MlModel::Cpg, row.and_then(|f| row_to_array(f.as_row())))
    } else if current.denovo_side().is_some() {
        let row = standard::DenovoCpgFeatures::extract(current, before, after);
        (MlModel::DenovoCpg, row.and_then(|f| row_to_array(f.as_row())))
    } else {
        let row = standard::OthersFeatures::extract(current, before, after);
        (MlModel::Others, row.and_then(|f| row_to_array(f.as_row())))
    }
}

/// Wrap a flat feature row into a `(1, N)` `Array2` for the forest predictor.
fn row_to_array(row: &[f32]) -> Result<Array2<f32>> {
    Array2::from_shape_vec((1, row.len()), row.to_vec())
        .wrap_err("Failed to build feature array from row")
}

#[cfg(test)]
mod tests {
    use super::shared::{CommonSectionA, CommonSectionB};
    use super::standard::{CpgFeatures, DenovoCpgFeatures, OthersFeatures};

    /// One name per `f32` slot, for every feature struct. A drift here means
    /// the `names()` output no longer aligns with the feature vectors.
    #[test]
    fn names_align_with_feature_count() {
        assert_eq!(CommonSectionA::names().len(), CommonSectionA::FEATURES);
        assert_eq!(CommonSectionB::names().len(), CommonSectionB::FEATURES);
        assert_eq!(CpgFeatures::names().len(), CpgFeatures::FEATURES);
        assert_eq!(DenovoCpgFeatures::names().len(), DenovoCpgFeatures::FEATURES);
        assert_eq!(OthersFeatures::names().len(), OthersFeatures::FEATURES);
    }

    /// The feature vector layout is frozen by every trained model. These counts
    /// are the layout the bundled model was fitted on; changing one
    /// invalidates the corresponding model, so they are pinned here.
    #[test]
    fn feature_counts_are_stable() {
        assert_eq!(CommonSectionA::FEATURES, 33);
        assert_eq!(CommonSectionB::FEATURES, 23);
        assert_eq!(CpgFeatures::FEATURES, 60);
        assert_eq!(DenovoCpgFeatures::FEATURES, 61);
        assert_eq!(OthersFeatures::FEATURES, 59);
    }

    /// Snapshot the full name→index layout of every model: reordering a field
    /// changes the snapshot and fails the test.
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
        ]
        .join("\n");

        insta::assert_snapshot!(report);
    }
}
