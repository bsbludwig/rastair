---
name: write-vcf-tests
description: Write comprehensive VCF integration tests for Rastair's methylation and variant calling. Creates tests using the pileups! macro and vcf_assert! to verify expected VCF output for complex CpG and denovo scenarios.
---

# Write VCF Tests for Rastair

Guide for writing VCF integration tests that verify Rastair's methylation calling and variant detection logic.

## Test Structure Overview

VCF tests verify the complete pipeline from pileup reads to VCF records:
1. Create synthetic pileup data with `pileups!` macro
2. Run through the calling pipeline with `test_call()`
3. Optionally modify ML scores with `set_pass()`/`set_fail()`
4. Reprocess with `reprocess()` to recalculate genotypes/methylation
5. Assert expected VCF output with `vcf_assert!`

## The pileups! Macro

Creates synthetic read data with clear, readable syntax:

```rust
let (segment, pileups) = pileups!(
    [ C G ] Ref,        // Reference bases
    [ T G ] OT,         // Read 1: C→T on OT strand
    [ T G ] OT,         // Read 2: C→T on OT strand
    [ C A ] OB,         // Read 3: G→A on OB strand
);
```

### Key Patterns

**Strand semantics in TAPS:**
- `OT` (Original Top): C→T indicates methylated C
- `OB` (Original Bottom): G→A indicates methylated G

**Common scenarios:**
```rust
// Original CpG with methylation
[ C G ] Ref,
[ T G ] OT,  // Methylation on C-side
[ C A ] OB,  // Methylation on G-side

// Denovo CpG (A→G creates CpG with previous C)
[ C A ] Ref,
[ C G ] OB,  // G variant creates denovo CpG

// Het variant
[ C G ] Ref,
[ C G ] OT,  // Ref allele
[ A G ] OT,  // Alt allele
```

## The vcf_assert! Macro

Defines expected VCF records with optional field assertions:

```rust
let expected_vcf = vcf_assert![
    (C .) PASS M5mC=0.5,           // Ref with no alt, methylation 50%
    (C T) FAIL,                     // C→T fails ML threshold
    (G A) PASS M5mC=1.0 GT="0/1",  // Het G/A, fully methylated
];
```

### Field Assertions

- `M5mC=<float>` - Single beta value
- `M5mC=vec![<f1>, <f2>]` - Dual context (original + denovo)
- `M5mC=None` - No methylation info
- `GT="0/1"` - Genotype (0/0, 0/1, 1/1, 1/2, etc.)
- `ML=<float>` - ML score

## Test Workflow Patterns

### Basic Test

```rust
#[test]
fn test_name() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ C G ] OT,
    );

    let records = test_call(segment, pileups, RecordFilters::all())?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.5,
        (G .) PASS M5mC=0.0,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
```

### With ML Score Manipulation

```rust
#[test]
fn test_with_ml_control() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ A G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // T is methylation evidence
    set_pass(&mut records[0], A); // A is real variant
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C A) PASS M5mC=0.5 GT="0/1",
        (C T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
```

## Edge Case Patterns

### Dual-Role Positions

When a position is both an original CpG side AND creates a denovo CpG:

```rust
// CGG with G→C at middle position
// Position 1 is G-side of original (0-1) AND C-side of denovo (1-2)
let (segment, pileups) = pileups!(
    [ C G G ] Ref,
    [ C C G ] OT,  // C variant at pos 1
    [ C T G ] OT,  // T = methylation on denovo C-side
    [ C A G ] OB,  // A = methylation on original G-side
);

// Expected: dual methylation values [original, denovo]
let expected_vcf = vcf_assert![
    (C .) PASS,
    (G C) PASS M5mC=vec![1.0, 0.5],  // Both contexts
    (G T) FAIL,
    (G A) FAIL,
    (G .) PASS,
];
```

### Compound Het

When both alts create different denovo CpGs:

```rust
// CAG with C/G het - each alt creates different denovo CpG
let (segment, pileups) = pileups!(
    [ C A G ] Ref,
    [ C C G ] OT,  // C creates denovo at 1-2
    [ C G G ] OB,  // G creates denovo at 0-1
);

let expected_vcf = vcf_assert![
    (C .) PASS,              // C-side of denovo 0-1
    (A C) PASS GT="1/2",     // Compound het
    (A G) PASS GT="1/2",
    (G .) PASS,              // G-side of denovo 1-2
];
```

### HomAlt Scenarios

HomAlt on original vs denovo CpGs behaves differently:

```rust
// Original CpG + HomAlt → beta=0.0 (ref base gone)
let (segment, pileups) = pileups!(
    [ C G ] Ref,
    [ T G ] OT,  // All T reads
);
set_pass(&mut records[0], T);
// Expected: M5mC=0.0

// Denovo CpG + HomAlt → normal beta (CpG on both chromosomes)
let (segment, pileups) = pileups!(
    [ A G ] Ref,
    [ C G ] OT,  // All C reads (denovo)
    [ T G ] OT,  // T = methylation
);
// Expected: M5mC=0.333 (normal beta calculation)
```

## Confounding Base Scenarios

When the confounding base (T for C-side, A for G-side) is also the ref:

```rust
// CAG with A→G het: A is ref AND confounding for G-side methylation
let (segment, pileups) = pileups!(
    [ C A G ] Ref,
    [ C A G ] OB,  // A reads: ref allele OR methylation? Ambiguous!
    [ C G G ] OB,  // G reads: alt allele
);

// Cannot distinguish ref A from methylation A
// Test should focus on unambiguous genotype call
```

## Test Organization

Organize tests by scenario type in separate modules:
- `basic.rs` - Simple variants outside CpG context
- `cpgs.rs` - Original CpG methylation scenarios
- `denovo.rs` - Denovo CpG creation
- `edge_cases.rs` - Subtle corner cases (Het non-confounding, HomAlt, etc.)
- `genotyping.rs` - Genotype estimation accuracy
- `ben_edge_cases.rs` - Complex dual-role and compound het scenarios

## Common Pitfalls

1. **Don't set ref base as pass/fail** - Only alt bases can be marked
2. **Test both strands** - OT and OB show different methylation patterns
3. **Document expected values** - Explain why beta=X is expected
4. **Consider genotype impact** - Het/HomAlt affects methylation calculation
5. **Match read counts** - Ensure read distribution supports expected genotype

## Test Naming Convention

Use descriptive names that indicate:
- Context (cgg, cag, cpg)
- Variant type (g_to_c, het_a_g)
- Key feature (methylation_on_both, dual_denovo)

Examples:
- `cgg_middle_g_to_c_methylation_on_both`
- `het_with_non_confounding_alt_original_cpg`
- `cag_middle_het_c_g_dual_denovo`

## Verifying Tests

After writing tests, check:
1. Compilation: `cargo test --lib <test_module> --no-run`
2. Run tests: `cargo test --lib <test_module>`
3. Review failures for insights about current implementation
4. Document expected vs actual behavior in test comments

Tests are documentation - failing tests show what needs implementation.
