# Reproducibility and scientific provenance

## Manuscript and package

The associated manuscript is *Human and bats genome robustness under COSMIC mutational signatures*. Sinabro supplies the sequential mutation simulator. The [separate analysis repository](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds) supplies preprocessing and biological analyses.

Version 3.0.0 is the consolidated, refactored implementation for ongoing research and revision. The mutation → evaluation → trajectory workflow is retained, but coordinate, sampling, termination, and random-state defects have been corrected. Therefore this edition must not be described as the exact code that produced historical manuscript outputs. Full old/new biological-output validation remains a separate research task.

The historical repository snapshot is commit `4c06417c2df7a12dfac7858540b2947054d29830`, available on [backup/before-cleanup-2026-10-14](https://github.com/StudyingAnt/sinabro/tree/backup/before-cleanup-2026-10-14). This preserves the source as received; the correspondence between every manuscript result and a specific historical commit still needs confirmation from the analysis records.

## Recorded environment

The [original conda package listing](reproducibility/original-packages.txt) records:

| Component | Recorded version |
|---|---|
| Python | 3.9.15 |
| NumPy | 1.21.6 |
| Biopython | 1.79 |
| SciPy | 1.8.1 |

The root `environment.yml` specifies these four scientific/runtime versions. The original listing contains additional tools from the broader research environment and is a provenance record, not the package's minimal dependency list. No conda channel release of Sinabro is claimed; conda installs the dependencies and pip installs this checkout.

The preceding refactor's complete package suite passed in the recorded Python/NumPy/Biopython environment. Separate Biopython 1.79, 1.80, and 1.81 checks found identical numeric payloads for 500 synthetic trajectories per environment, codon translation, BLOSUM80 scores, and selected pairwise dN/dS examples. This evidence concerns the tested cases, not all biological inputs or future library releases. The consolidation is retested as documented in [testing](testing.md).

## Criteria and output interpretation

Protein identity is the strict criterion used for synonymous robustness. BLOSUM distance is a relaxed criterion that accepts some protein changes; its stopping time must be reported separately from synonymous RUMS. Both are sequence-level modeling criteria rather than experimental measurements of fitness. See [algorithm](algorithm.md) for thresholds and stopping behavior.

Before using the refactor to revise manuscript results, compare a fixed biological subset against the historical implementation with identical input sequences, signature tables, preprocessing, run counts, and evaluation definitions. Record differences caused by package corrections separately from those caused by a relaxed criterion.

## Archive each scientific run

Record the source commit and package version, Python/NumPy/Biopython/SciPy versions as applicable, input FASTA hash, profile hash, CDS bounds, seed, trajectory count, evaluation threshold, stop policy, and incomplete-path counts. The CLI generates a manifest alongside its CSV. Store the analysis-script commit and preprocessing parameters with downstream results as well.

## Restore historical code separately

An independent checkout keeps the current package intact:

```bash
git fetch origin
git worktree add ../sinabro_original origin/backup/before-cleanup-2026-10-14
```

Run it in a separate environment and retain the original input/output files. Historical scripts and notebooks are preserved in that Git history, rather than included in the current installable package.
