# Change log

## 3.0.0 — 2026-10-14

This edition consolidates the refactored code into one installable `sinabro` package. Removing historical import paths is a breaking interface change; the scientific engine and packaged SBS tables are retained.

- Use `import sinabro as snbr` for trajectories and robustness computations.
- Move protein distance and optional pairwise dN/dS to `sinabro.metrics`.
- Remove `sinabro_v2`, `sinabro.sinabro`, and the historical helper, mutation, and evaluation facades.
- Keep one copy of COSMIC SBS data under `sinabro/data`.
- Preserve historical development scripts, notebooks, unused assets, and the earlier checkout outside the active package tree and in backup branches.
- Provide concise installation, examples, API, testing, and scientific-provenance documentation.
- Include a conda specification using the recorded Python, NumPy, and Biopython versions.
- Permit Biopython 1.79 installation, following compatibility checks of 1.79, 1.80, and 1.81.

The six tests specific to removed compatibility interfaces were retired. Mutation, coordinate, evaluation, statistical, metric, and CLI behavior checks remain.

## 2.1.0 — 2026-10-14

Initial refactored implementation: shared engine, explicit coding coordinates, validated weighted mutation opportunities, seeded random streams, protein-identity and BLOSUM evaluations, complete-path aggregation, portable CLI, and automated tests. Historical import facades were temporarily retained.

This implementation corrects behavior in the preserved historical source. Its output is not assumed to reproduce earlier manuscript results without an impact assessment. See [reproducibility](docs/reproducibility.md).
