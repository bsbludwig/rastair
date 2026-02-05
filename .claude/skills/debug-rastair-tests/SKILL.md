---
name: debug-rastair-tests
description: Debug failing VCF tests in src/call/tests/vcf_tests/. Use when methylation or genotype test failures occur, or when investigating beta calculation issues.
---

# Debug Rastair Methylation/Genotype Tests

Systematic approach to debugging test failures in `src/call/tests/vcf_tests/`.

## Understanding Test Structure

### Pileup Format
```rust
pileups!(
    [ C G ] Ref,      // Reference bases
    [ T A ] OT,       // Read 1: pos0=T, pos1=A on OT strand
    [ C G ] OB,       // Read 2: pos0=C, pos1=G on OB strand
)
```

**Key:**
- Each row = one read across all positions
- `Ref` = reference sequence
- `OT` = Original Top (forward strand)
- `OB` = Original Bottom (reverse strand)

### TAPS Biology Essentials
- **C position methylation:** C→T on OT strand
- **G position methylation:** G→A on OB strand
- Strand matters: C uses OT, G uses OB

### Test Modifications
- `set_pass(records[i], Base)` → RealVariant (ml=1.0)
- `set_fail(records[i], Base)` → MethylationEvidenceOnly (ml=0.0)
- `reprocess(records)` → recalculates genotype, methylation, strand info

## Debugging Process

### Step 1: Add Debug Output

Insert before the failing assertion:

```rust
eprintln!("\n=== DEBUG: $0 ===");
for (i, rec) in records.iter().enumerate() {
    eprintln!("\nPos {}: ref={:?}", i, rec.ref_base());
    eprintln!("  Genotype: {:?}", rec.pos_metrics.genotype);
    eprintln!("  Methylated: {:?}", rec.pos_metrics.methylated);
    eprintln!("  Ref strand: {:?}", rec.ref_metrics.strand_count);
    for alt in &rec.alts {
        eprintln!("  Alt {:?}: call={:?}, ml={:?}, strand={:?}",
            alt.base, alt.call, alt.filters.ml, alt.metrics.strand_count);
    }
}
```

### Step 2: Run Test

```bash
cargo test $0 -- --nocapture 2>&1 | grep -A 50 "=== DEBUG"
```

### Step 3: Interpret Output

**Genotype values:**
- `HomRef` = 0/0 (reference)
- `HomAlt(n)` = n/n (e.g., HomAlt(1) = 1/1 for first alt)
- `RefHet(n)` = 0/n (heterozygous)
- `AltHet(a,b)` = a/b (two alts)

**Alt call types:**
- `RealVariant` → counts for genotype
- `MethylationEvidenceOnly { for_base }` → should count as ref base
- `ReadError` → ignored

**Methylated variants:**
- `OriginalCpG(beta)` → beta for original CpG
- `DeNovoCpG(beta)` → beta for de-novo CpG
- `Both { original_beta, denovo_beta }` → both values
- `NoEvidence` / `Unknown` → no methylation data

### Step 4: Check Expectations

Compare debug output to test assertions:

```rust
vcf_assert![
    (C A) PASS M5mC=0.5 GT="0/1",
    //   ^alt   ^filter  ^beta  ^genotype
];
```

**Common mismatches:**
- Wrong genotype → MethylationEvidence not counted
- Beta=0.0 (expect >0) → likely hom alt instead of het
- Beta=1.0 (expect <1) → missing unmethylated reads
- Wrong strand counts → reads on unexpected strand

### Step 5: Trace Logic

**For ref C:**
1. OT strand: C (unmod) vs T (mod)
2. If het C/X: T = methylation of C allele
3. If hom X/X: beta = 0 (no C)

**For ref G:**
1. OB strand: G (unmod) vs A (mod)
2. If het G/X: A = methylation of G allele
3. If hom X/X: beta = 0 (no G)

### Step 6: Verify Calculations

**Beta formula:**
```
beta = mod_reads / (mod_reads + unmod_reads)
```

**Het SNP adjustment:** `calculate_het_snp_beta()` accounts for 50% expected alt reads

**Genotype should count:**
- RealVariant alts
- MethylationEvidence as ref: T→C, A→G

## Quick Reference

### Strand-Specific Patterns
| Ref | Methylated | Unmethylated | Strand |
|-----|------------|--------------|--------|
| C   | T          | C            | OT     |
| G   | A          | G            | OB     |

### Common Test Patterns
- **Simple methylation:** All T (at C) or all A (at G)
- **Het + methylation:** Real variant + methylation evidence
- **De-novo CpG:** X→C or X→G creates new CpG
- **Dual role:** Position in original + de-novo CpG

## Common Fixes

1. **Wrong genotype:** Fix `estimate_genotype()` to count MethylationEvidence
2. **Wrong beta:** Check genotype first, then strand counting
3. **Test expectation:** Update expected values for correct behavior
4. **Strand asymmetry:** Reads on unexpected strand

## Cleanup

Remove debug `eprintln!` statements after fixing.
