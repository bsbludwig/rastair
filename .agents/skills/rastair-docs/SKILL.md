---
name: rastair-docs
description: >
  Look up Rastair project documentation. Use when you need explanations of domain concepts (TAPS, methylation, CpG, genotyping, de-novo CpGs), output formats (VCF, BED, BAM), calling logic (variants, methylation algorithm, ML, de-novo detection), CLI usage, or any other Rastair-specific topic. Directs you to search docs/src/ for answers.
---

When you need to understand something about Rastair — a concept, format, algorithm, or CLI flag — look in `docs/src/` before guessing or relying on memory.

## Finding the right doc

**Glossary first**: For any unfamiliar term (TAPS, OT/OB, VAF, CpG, de-novo, methylation, SNV…), search in `docs/src/glossary.yaml` — it defines all domain terminology used across the project in 3-10 lines each.

If you need more Rastair-specific information, do this:

1. **Orientation**: `docs/src/SUMMARY.md` lists all docs with their topics — read it first to identify which file covers your question.

2. **Then read the relevant file**

3. **When in doubt, search**: Use Grep to search across all docs for a keyword rather than guessing.

## Keeping docs current

When a session reveals knowledge that isn't captured in `docs/src/` — a new behaviour, edge case, algorithm detail, field meaning, or corrected understanding — proactively prompt the user:

> "I learned something about `<topic>` that isn't in the docs yet. Would you like me to add it to `<suggested file>`?"

Use your judgement about what's worth documenting:
- **Do suggest**: newly understood algorithm details, edge cases in calling logic, corrected misconceptions in existing docs
- **Skip**: session-specific debugging steps, temporary workarounds, things already well covered

Run `cargo xtask docs` before updating the docs to make sure the generated files are up-to-date.
