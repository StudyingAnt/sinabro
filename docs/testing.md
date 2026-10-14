# Package testing

Install the optional analysis dependencies to run the complete suite:

```bash
python -m pip install -e ".[analysis]"
python -m unittest discover -s tests -v
```

| Test file | Main checks |
|---|---|
| `test_coordinates.py` | All 64 codons × nine substitutions, frame offsets, flank exclusion, record invariants |
| `test_mutation.py` | Motifs, overlapping opportunities, strand weights, SBS profiles, invalid inputs, deterministic sampling |
| `test_trajectory_evaluation.py` | Protein identity, BLOSUM thresholds, custom policies, stopping reasons, incomplete paths |
| `test_robustness.py` | Stopping-time and endpoint aggregates, matched criterion prefixes, serial/parallel equality |
| `test_metrics.py` | Protein distances, stops, pairwise dN/dS, undefined/saturated cases |
| `test_cli.py` | CSV and manifest contents, reproducibility, validation, overwrite protection |

There are 89 automated tests after removing six checks for interfaces that are no longer shipped. The exhaustive coordinate oracle additionally evaluates 576 codon/substitution combinations inside one test.

GitHub Actions is configured for Linux and Windows with Python 3.9 and 3.12, plus a pinned scientific-library job for NumPy 1.21.6 and Biopython 1.79. The pinned job uses Python 3.9 as provided by the runner; it does not claim the precise 3.9.15 patch version. Check the actual Actions result for the relevant commit before citing remote CI as passed.

Package validation also includes building a wheel, checking that it contains only the `sinabro` namespace and five SBS tables, and running the suite and CLI against an isolated wheel installation outside the source checkout. Synthetic checks verify software behavior; they do not reproduce manuscript figures.
