# Migration from the preserved original to 2.1.0

The refactor keeps Trajectory → mutation → evaluation and the separate manuscript pipeline. It consolidates numerical code behind public and compatibility APIs. The historical implementation remains in Git at `backup/package-before-refactor-2026-10-14` (commit `4c06417c2df7a12dfac7858540b2947054d29830`). Its source snapshot is also stored in the workspace's `codex_work/sinabro_review_2026-10-14/sinabro_before_refactor.tar.gz`.

## Deliberate behavior corrections

| Topic | Before | 2.1.0 |
|---|---|---|
| Public import | README import lacked Trajectory export | `import sinabro as snbr` works |
| Coordinates | Mixed raw index, HGVS, start offset | Internal zero-based target, explicit CDS slice |
| CDS-plus defaults | Legacy default start=0 could evaluate wrong codon | CDS starts at 1 and ends at len(sequence)-1 exclusive |
| `end` | Backend-dependent inclusive/exclusive semantics | Always exclusive; negative coordinates rejected |
| Plain CDS | Implicit/inconsistent frame | Explicit `cds_start=0, cds_end=len(sequence)` |
| RNG | Time-based reseeding and global RNG | Per-trajectory Generator; explicit seed |
| Motifs | Non-overlapping regex matches | All overlapping opportunities |
| Strand bias | Same multiplier on both strands | Forward bias and reverse 1-bias |
| Empty probabilities | NaNs or invalid codon evaluation | Explicit no_available_event |
| Counts | Blind length-1 even for no-motif/capped path | accepted_steps=None if incomplete |
| BLOSUM | Distributed development code and mixed conventions | Original-vs-current score; threshold >=; explicit stop policy |
| Legacy BLASTP mode | Fixed temporary files and unsupported calls | Explicit error; no silent BLOSUM substitution |
| `max_length` | Appended up to that many extra steps | Total target number of substitutions |
| `add_mutated_sequence` | Legacy returned sequence list | Returns MutInfo; trajectory data still available |
| Manual append | Parallel arrays could diverge | Exactly one in-CDS substitution; synchronized record arrays |
| Fixed-n v2 | Serial argument and missing-method bugs | Explicit endpoint phenotype evaluation |

Corrected results may differ from the original. Do not substitute the new package into a historical manuscript calculation without an old/new impact assessment. The original mutation-profile interpretation (per-position weights followed by global normalization) is retained; no genome-opportunity correction or biological null model has been silently introduced.

## v2 interface

```python
import sinabro_v2 as snbr
record = snbr.MutationRecord("TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT")
trajectory = snbr.Trajectory(0, record, seed=42)
trajectory.autofill("signature", "blosum", mutational_signature="SBS2", threshold=20)
```

RobustnessComputer and its generated trajectories use the shared engine. Their getters, records, lengths, and result dictionary remain available. Directly generated engine trajectories use the public `autofill(condition=..., method=...)` form, rather than the facade's positional v2 form. This distinction is documented rather than hidden.

Pairwise dN/dS remains an optional exploratory metric. It needs complete, equal-length CDS with a reference terminal stop and is not dNdScv or a test of observed tumor selection. Hamming in the v2 metric module is translated protein Hamming; DNA Hamming is not inferred from that name. Excluding a terminal record requires at least two records.

## Safe restoration

After the refactor is committed and the working tree is clean, an independent checkout is preferable to moving the working branch:

```bash
git worktree add ../sinabro_original backup/package-before-refactor-2026-10-14
```

This leaves the new checkout intact and opens the original source beside it. Do not use `git reset --hard` to switch between scientific versions. Run historical code in its own environment and preserve output/input hashes.
