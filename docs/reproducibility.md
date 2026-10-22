# Scientific provenance and reproducibility

## Package and manuscript analyses

Sinabro 3.0.0 implements the mutation → evaluation → trajectory workflow used by the associated [human/bat manuscript analysis repository](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds). The scientific source and signature tables are unchanged by the publication-oriented repository organization.

The analysis workflows have been checked on the full human COSMIC campaign (73,214 transcripts × 79 signatures × 1,000 paths), two bat COSMIC cohorts, and six experimental enzyme profiles on all three species. Independent sampler checks, downstream numerical comparisons and the BRCA1 example are documented in the [analysis validation reference](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/blob/main/docs/validation_reference.md). These analyses use the corrected implementation; they do not establish identity with unrecovered historical random streams.

## Recorded environment

| Component | Version |
| --- | --- |
| Python | 3.9.15 |
| NumPy | 1.21.6 |
| Biopython | 1.79 |
| SciPy | 1.8.1 |

The package `environment.yml` records these dependencies. The analysis repository additionally pins Matplotlib 3.6.2. The [historical environment listing](reproducibility/original-packages.txt) is retained as source provenance; it contains tools beyond the package's minimal dependencies. Conda installs dependencies and pip installs this checkout; a published Conda channel package is not assumed.

Biopython 1.79, 1.80 and 1.81 checks gave identical numeric payloads for the tested 500 synthetic trajectories per environment, translations, BLOSUM80 scores and selected pairwise dN/dS examples. This evidence is limited to those cases and versions. Use the pinned analysis environment for manuscript reproduction.

## Scientific criteria

Protein identity defines synonymous RUMS. BLOSUM distance accepts some protein changes below an explicit alignment-score threshold and is reported as a separate criterion-neutral stopping time. Neither is an experimental assay of fitness. Incomplete paths have no complete robustness estimate and must not be silently averaged. See [the algorithm](algorithm.md).

## Historical source

The preserved earlier implementation is commit `4c06417c2df7a12dfac7858540b2947054d29830`, available on [backup/before-cleanup-2026-10-14](https://github.com/StudyingAnt/sinabro/tree/backup/before-cleanup-2026-10-14). Coordinate, opportunity-sampling, termination and random-state corrections distinguish the present edition; [migration notes](migration.md) describe them. The current package must not be labelled as the exact executable that generated all historical figures.

The publication cleanup snapshot is available on [backup/before-paper-release-2026-10-21](https://github.com/StudyingAnt/sinabro/tree/backup/before-paper-release-2026-10-21). Historical code can be inspected separately without changing the current checkout:

```bash
git fetch origin
git worktree add ../sinabro_original origin/backup/before-cleanup-2026-10-14
```

## Record each scientific run

Preserve the package and analysis commits, source-file fingerprints, Python/scientific-library versions, input FASTA and signature hashes, CDS bounds, seeds, path counts, evaluation settings, stop policy and incomplete-path counts. The CLI writes a manifest beside its output. Use a new output directory for each distinct scientific configuration and retain the resolved environment with the campaign.
