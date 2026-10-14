# Moving to the consolidated package

Version 3.0.0 exposes one package:

```python
import sinabro as snbr

trajectory = snbr.Trajectory("example", "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT", seed=42)
trajectory.autofill(
    condition="nonsynonymous", method="signature", mutational_signature="SBS2",
)
```

Historical `sinabro_v2`, `sinabro.sinabro`, `_helper`, `_mutate`, and evaluation facade imports have been removed. Earlier implementations remain in backup branches. Use the public package API and [API guide](api.md) for new code.

## Useful metrics retained

```python
from sinabro.metrics import hamming, compute_distance, compute_dnds
```

`hamming` compares translated proteins from CDS strings. `compute_distance` uses trajectory records and explicit CDS bounds. Pairwise `compute_dnds` is optional, requires suitable complete CDS with a reference terminal stop, and is an exploratory sequence metric rather than tumor-selection inference. Install `sinabro[analysis]` to supply SciPy.

## Deliberate corrections relative to the historical source

| Topic | Consolidated behavior |
|---|---|
| Coordinates | Zero-based internal targets; coding ends are exclusive |
| CDS-plus | Default coding slice starts at 1 and ends at sequence length minus 1 |
| Plain CDS | Set `cds_start=0, cds_end=len(sequence)` explicitly |
| Randomness | Per-trajectory NumPy Generator; explicit seeds |
| Motifs | Overlapping opportunities are eligible |
| Strand weights | Forward bias and reverse 1-bias |
| No event | Explicit `no_available_event` rather than invalid probability values |
| Stopping time | Incomplete paths have `accepted_steps=None` and cannot be silently averaged |
| BLOSUM | Original/current comparison; distance >= positive threshold; explicit stop policy |
| `max_length` | Total target substitutions rather than extra appended steps |
| Mutation return | `add_mutated_sequence` returns `MutInfo`; records remain accessible |
| Fixed-n measure | Protein identity after exactly n steps, including possible restoration |

The current scientific engine is unchanged by the folder consolidation. These corrections were introduced in the preceding refactor and can change outputs relative to the preserved historical implementation. See [scientific provenance](reproducibility.md) before replacing manuscript calculations.
