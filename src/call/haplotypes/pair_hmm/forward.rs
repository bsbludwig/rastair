//! The scaled forward algorithm of the local `PairHMM`, run on several
//! problems at a time.
//!
//! Every lane is an independent problem and performs exactly the operations of
//! the scalar recurrence, in the same order, so a batch scores its members
//! bit-for-bit as scoring them one at a time would. That is what lets the
//! batching, the phasing memo and the target-feature dispatch below leave
//! results unchanged.

use color_eyre::eyre::{OptionExt as _, Result, ensure};
use rayon::prelude::*;
use std::collections::BTreeMap;

/// Observation terms of one profile row: the sequencing emission under each
/// known template base, then the insertion emission.
pub(super) const TERMS: usize = 5;
/// Slot of the insertion emission. A column standing for an unknown template
/// base weights this slot alone, so it emits the insertion probability.
pub(super) const INSERTION: usize = 4;

pub(super) type Row = [f64; TERMS];
pub(super) type Column = [f64; TERMS];

/// A profile against a template. The emission at `(i, j)` is the slot-wise
/// product of `profile[i]` and `template[j]`, summed in slot order.
#[derive(Clone, Copy)]
pub(super) struct Problem<'a> {
    pub profile: &'a [Row],
    pub template: &'a [Column],
}

#[derive(Clone, Copy, PartialEq)]
pub(super) struct Gaps {
    pub open: f64,
    pub extend: f64,
}

/// The widest batch: one AVX-512 register of `f64` per lane-vector.
const LANES: usize = 8;

/// The log-likelihood of every problem, in order.
///
/// Problems of one shape share a batch. Batches run on the enclosing rayon
/// pool when there is one, so a worker phasing a dense segment borrows the
/// idle workers; the lane arithmetic makes the split invisible in the result.
pub(super) fn forward_all(problems: &[Problem<'_>], gaps: Gaps) -> Result<Vec<f64>> {
    let mut shapes: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (index, problem) in problems.iter().enumerate() {
        shapes.entry((problem.profile.len(), problem.template.len())).or_default().push(index);
    }
    let batches: Vec<Vec<usize>> = shapes
        .into_values()
        .flat_map(|members| members.chunks(LANES).map(<[usize]>::to_vec).collect::<Vec<_>>())
        .collect();
    let run = |batch: &Vec<usize>| -> Result<Vec<f64>> {
        let selected = batch
            .iter()
            .map(|&index| problems.get(index).copied().ok_or_eyre("Missing forward problem"))
            .collect::<Result<Vec<_>>>()?;
        run_slice(&selected, gaps)
    };
    // A few batches per task keep rayon's per-job cost below the DP's.
    let scored: Vec<Result<Vec<f64>>> =
        if rayon::current_thread_index().is_some() && batches.len() > TASK_BATCHES {
            batches
                .par_chunks(TASK_BATCHES)
                .flat_map_iter(|chunk| chunk.iter().map(run).collect::<Vec<_>>())
                .collect()
        } else {
            batches.iter().map(run).collect()
        };
    let mut results = vec![f64::NAN; problems.len()];
    for (batch, scores) in batches.iter().zip(scored) {
        for (&index, score) in batch.iter().zip(scores?) {
            *results.get_mut(index).ok_or_eyre("Forward result outside its batch")? = score;
        }
    }
    Ok(results)
}

/// Batches one rayon task scores before another thread can steal.
const TASK_BATCHES: usize = 4;

/// Splits a batch into power-of-two lane groups, widest first.
fn run_slice(problems: &[Problem<'_>], gaps: Gaps) -> Result<Vec<f64>> {
    let mut results = Vec::with_capacity(problems.len());
    let mut rest = problems;
    while !rest.is_empty() {
        let (head, tail) = if let Ok(head) = <&[Problem<'_>; 8]>::try_from(rest) {
            (run(head, gaps)?.to_vec(), rest.get(8..))
        } else if let Some(head) =
            rest.get(..4).and_then(|head| <&[Problem<'_>; 4]>::try_from(head).ok())
        {
            (run(head, gaps)?.to_vec(), rest.get(4..))
        } else if let Some(head) =
            rest.get(..2).and_then(|head| <&[Problem<'_>; 2]>::try_from(head).ok())
        {
            (run(head, gaps)?.to_vec(), rest.get(2..))
        } else {
            let head = rest.first().copied().ok_or_eyre("Empty forward batch")?;
            (run(&[head], gaps)?.to_vec(), rest.get(1..))
        };
        results.extend(head);
        rest = tail.ok_or_eyre("Forward batch split outside its slice")?;
    }
    Ok(results)
}

/// Dispatches to the widest vector unit the CPU has. Every variant is the same
/// generic function and no fused operations are enabled, so they agree bit for
/// bit.
fn run<const L: usize>(problems: &[Problem<'_>; L], gaps: Gaps) -> Result<[f64; L]> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f") {
            // SAFETY: the feature was detected on the running CPU just above.
            return unsafe { run_avx512(problems, gaps) };
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the feature was detected on the running CPU just above.
            return unsafe { run_avx2(problems, gaps) };
        }
    }
    run_lanes(problems, gaps)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
fn run_avx512<const L: usize>(problems: &[Problem<'_>; L], gaps: Gaps) -> Result<[f64; L]> {
    run_lanes(problems, gaps)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn run_avx2<const L: usize>(problems: &[Problem<'_>; L], gaps: Gaps) -> Result<[f64; L]> {
    run_lanes(problems, gaps)
}

type Lanes<const L: usize> = [f64; L];

#[inline(always)]
fn mul<const L: usize>(a: &Lanes<L>, b: &Lanes<L>) -> Lanes<L> {
    let mut out = [0.0; L];
    for ((out, a), b) in out.iter_mut().zip(a).zip(b) {
        *out = a * b;
    }
    out
}

#[inline(always)]
fn add<const L: usize>(a: &Lanes<L>, b: &Lanes<L>) -> Lanes<L> {
    let mut out = [0.0; L];
    for ((out, a), b) in out.iter_mut().zip(a).zip(b) {
        *out = a + b;
    }
    out
}

#[inline(always)]
fn scaled<const L: usize>(a: &Lanes<L>, factor: f64) -> Lanes<L> {
    let mut out = [0.0; L];
    for (out, a) in out.iter_mut().zip(a) {
        *out = a * factor;
    }
    out
}

/// `f64::max` folded from zero skips a NaN; so does this comparison.
#[inline(always)]
fn raise<const L: usize>(maximum: &mut Lanes<L>, values: &Lanes<L>) {
    for (maximum, &value) in maximum.iter_mut().zip(values) {
        if value > *maximum {
            *maximum = value;
        }
    }
}

/// Transposes one row or column of every lane's problem into lane vectors.
#[inline(always)]
fn gather<const L: usize>(terms: &[Option<&Row>; L]) -> Result<[Lanes<L>; TERMS]> {
    let mut vectors = [[0.0; L]; TERMS];
    for (lane, terms) in terms.iter().enumerate() {
        let terms = terms.ok_or_eyre("Forward problem shorter than its batch")?;
        for (vector, &value) in vectors.iter_mut().zip(terms) {
            *vector.get_mut(lane).ok_or_eyre("Forward lane outside its batch")? = value;
        }
    }
    Ok(vectors)
}

#[inline(always)]
fn run_lanes<const L: usize>(problems: &[Problem<'_>; L], gaps: Gaps) -> Result<[f64; L]> {
    let first = problems.first().ok_or_eyre("Empty forward batch")?;
    let (rows, columns) = (first.profile.len(), first.template.len());
    ensure!(
        problems
            .iter()
            .all(|problem| problem.profile.len() == rows && problem.template.len() == columns),
        "A forward batch needs problems of one shape"
    );
    if rows == 0 {
        return Ok([0.0; L]);
    }
    ensure!(columns > 0, "Empty PairHMM haplotype");
    let width = columns.checked_add(1).ok_or_eyre("PairHMM width overflow")?;
    let mut weights = Vec::with_capacity(columns);
    for column in 0..columns {
        weights.push(gather(&problems.map(|problem| problem.template.get(column)))?);
    }
    let zero = [0.0; L];
    let mut previous_m = vec![zero; width];
    let mut previous_i = vec![zero; width];
    let mut previous_d = vec![[1.0 / columns as f64; L]; width];
    *previous_d.last_mut().ok_or_eyre("Missing final HMM cell")? = zero;
    let mut current_m = vec![zero; width];
    let mut current_i = vec![zero; width];
    let mut current_d = vec![zero; width];
    let (open, extend) = (gaps.open, gaps.extend);
    let (matched_to_matched, gap_to_matched) = (1.0 - 2.0 * open, 1.0 - extend);
    let mut scale = zero;
    for row in 0..rows {
        let [term_a, term_c, term_g, term_t, insertion] =
            gather(&problems.map(|problem| problem.profile.get(row)))?;
        let (first_m, rest_m) = current_m.split_first_mut().ok_or_eyre("Missing HMM cell")?;
        let (first_i, rest_i) = current_i.split_first_mut().ok_or_eyre("Missing HMM cell")?;
        let (first_d, rest_d) = current_d.split_first_mut().ok_or_eyre("Missing HMM cell")?;
        *first_m = zero;
        *first_i = zero;
        *first_d = zero;
        let (diagonal_m, above_m) = previous_m.split_first().ok_or_eyre("Missing HMM cell")?;
        let (diagonal_i, above_i) = previous_i.split_first().ok_or_eyre("Missing HMM cell")?;
        let (diagonal_d, above_d) = previous_d.split_first().ok_or_eyre("Missing HMM cell")?;
        let (mut diagonal_m, mut diagonal_i, mut diagonal_d) =
            (*diagonal_m, *diagonal_i, *diagonal_d);
        let (mut left_m, mut left_d) = (zero, zero);
        let mut maximum = zero;
        let cells = rest_m
            .iter_mut()
            .zip(rest_i.iter_mut())
            .zip(rest_d.iter_mut())
            .zip(above_m)
            .zip(above_i)
            .zip(above_d)
            .zip(&weights);
        for ((((((cell_m, cell_i), cell_d), above_m), above_i), above_d), weight) in cells {
            let [weight_a, weight_c, weight_g, weight_t, weight_insertion] = weight;
            let emission = add(
                &add(
                    &add(
                        &add(&mul(&term_a, weight_a), &mul(&term_c, weight_c)),
                        &mul(&term_g, weight_g),
                    ),
                    &mul(&term_t, weight_t),
                ),
                &mul(&insertion, weight_insertion),
            );
            let matched = mul(
                &emission,
                &add(
                    &add(
                        &scaled(&diagonal_m, matched_to_matched),
                        &scaled(&diagonal_i, gap_to_matched),
                    ),
                    &scaled(&diagonal_d, gap_to_matched),
                ),
            );
            let inserted = mul(&insertion, &add(&scaled(above_m, open), &scaled(above_i, extend)));
            let deleted = add(&scaled(&left_m, open), &scaled(&left_d, extend));
            raise(&mut maximum, &matched);
            raise(&mut maximum, &inserted);
            raise(&mut maximum, &deleted);
            *cell_m = matched;
            *cell_i = inserted;
            *cell_d = deleted;
            (diagonal_m, diagonal_i, diagonal_d) = (*above_m, *above_i, *above_d);
            (left_m, left_d) = (matched, deleted);
        }
        // A common row scale leaves every relative path weight intact and prevents underflow on long reads.
        ensure!(
            maximum.iter().all(|maximum| maximum.is_finite() && *maximum > 0.0),
            "PairHMM has no finite alignment path"
        );
        for cell in current_m.iter_mut().chain(&mut current_i).chain(&mut current_d) {
            for (value, maximum) in cell.iter_mut().zip(&maximum) {
                *value /= maximum;
            }
        }
        for (scale, maximum) in scale.iter_mut().zip(&maximum) {
            *scale += maximum.ln();
        }
        std::mem::swap(&mut previous_m, &mut current_m);
        std::mem::swap(&mut previous_i, &mut current_i);
        std::mem::swap(&mut previous_d, &mut current_d);
    }
    let mut mass = zero;
    for (cell_m, cell_i) in previous_m.iter().zip(&previous_i).skip(1) {
        mass = add(&mass, cell_m);
        mass = add(&mass, cell_i);
    }
    let mut values = zero;
    for ((value, scale), mass) in values.iter_mut().zip(&scale).zip(&mass) {
        *value = scale + mass.ln();
        ensure!(value.is_finite(), "PairHMM forward likelihood is not finite");
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    #![expect(clippy::cast_possible_truncation, reason = "Test code")]
    use super::*;

    fn problems(seed: u64, count: usize) -> Vec<(Vec<Row>, Vec<Column>)> {
        let mut state = seed;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..count)
            .map(|index| {
                let rows = 1 + (index * 37) % 23;
                let columns = 1 + (index * 53) % 29;
                let profile = (0..rows)
                    .map(|_| {
                        let error = 10.0_f64.powf(-(next() * 40.0) / 10.0);
                        let mut row = [error / 3.0; TERMS];
                        *row.get_mut(((next() * 4.0) as usize).min(3)).expect("slot") = 1.0 - error;
                        *row.get_mut(INSERTION).expect("slot") =
                            if next() < 0.5 { 0.25 } else { 1.0 };
                        row
                    })
                    .collect();
                let template = (0..columns)
                    .map(|_| {
                        let mut column = [0.0; TERMS];
                        match (next() * 4.0) as usize {
                            0 => *column.get_mut(INSERTION).expect("slot") = 1.0,
                            1 => {
                                *column.get_mut(1).expect("slot") = 0.5;
                                *column.get_mut(3).expect("slot") = 0.5;
                            }
                            2 => {
                                *column.get_mut(((next() * 4.0) as usize).min(3)).expect("slot") =
                                    1.0
                            }
                            _ => {
                                let rate = next();
                                *column.get_mut(2).expect("slot") = 1.0 - rate;
                                *column.get_mut(0).expect("slot") = rate;
                            }
                        }
                        column
                    })
                    .collect();
                (profile, template)
            })
            .collect()
    }

    fn gaps() -> Gaps {
        Gaps { open: 10.0_f64.powf(-4.5), extend: 0.1 }
    }

    #[test]
    fn batches_and_vector_units_agree_with_single_lanes_bit_for_bit() -> Result<()> {
        let owned = problems(0x9e37_79b9_7f4a_7c15, 96);
        let problems: Vec<_> =
            owned.iter().map(|(profile, template)| Problem { profile, template }).collect();
        let singles = problems
            .iter()
            .map(|&problem| {
                let [single] = run_lanes(&[problem], gaps())?;
                Ok(single)
            })
            .collect::<Result<Vec<_>>>()?;
        let batched = forward_all(&problems, gaps())?;
        assert_eq!(singles.len(), batched.len());
        for (single, batched) in singles.iter().zip(&batched) {
            assert_eq!(single.to_bits(), batched.to_bits(), "{single} versus {batched}");
        }
        for (index, &problem) in problems.iter().enumerate() {
            let [dispatched] = run(&[problem], gaps())?;
            assert_eq!(
                dispatched.to_bits(),
                singles.get(index).copied().expect("single").to_bits()
            );
        }
        let mut same_shape = Vec::new();
        for (profile, template) in &owned {
            if profile.len() == 1 + 37 % 23 && template.len() == 1 + 53 % 29 {
                same_shape.push(Problem { profile, template });
            }
        }
        if let Ok(eight) = <[Problem<'_>; 8]>::try_from(
            same_shape.iter().copied().cycle().take(8).collect::<Vec<_>>(),
        ) {
            let wide = run(&eight, gaps())?;
            for (lane, problem) in eight.iter().enumerate() {
                let [single] = run_lanes(&[*problem], gaps())?;
                assert_eq!(wide.get(lane).copied().expect("lane").to_bits(), single.to_bits());
            }
        }
        Ok(())
    }

    #[test]
    fn mismatched_shapes_and_empty_templates_are_refused() -> Result<()> {
        let (profile, template) = problems(7, 1).pop().expect("problem");
        let short = template.get(..template.len() - 1).unwrap_or_default();
        assert!(
            run_lanes(
                &[
                    Problem { profile: &profile, template: &template },
                    Problem { profile: &profile, template: short }
                ],
                gaps()
            )
            .is_err()
        );
        assert!(run_lanes(&[Problem { profile: &profile, template: &[] }], gaps()).is_err());
        assert_eq!(run_lanes(&[Problem { profile: &[], template: &template }], gaps())?, [0.0]);
        Ok(())
    }
}
