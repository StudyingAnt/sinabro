# Python API and command-line usage

Use `import sinabro as snbr`. All examples below use the single installed `sinabro` package.

## Protein-identity neutrality: the original RUMS criterion

CDS-plus input contains one context nucleotide before and after the CDS. The default coding slice is `sequence[1:-1]`. Coordinates are zero based and ends are exclusive.

```python
import sinabro as snbr

sequence = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
trajectory = snbr.Trajectory(id="example", data=sequence, seed=42)
trajectory.autofill(
    condition="nonsynonymous",
    method="signature",
    mutational_signature="SBS2",
)
print(trajectory.terminal_reason)
print(trajectory.accepted_steps)  # 0 when the first substitution is nonsynonymous
trajectory.show()
```

The synonymous count excludes the terminal nonsynonymous event. RUMS is its mean across complete independent trajectories. This is a sequence-level proxy, not a direct measurement of fitness or function.

```python
computer = snbr.RobustnessComputer("example", sequence, seed=42)
rums = computer.compute_l_robustness(
    n_sim=100, method="signature", eval_method="nonsynonymous",
    mutational_signature="SBS2",
)
print(rums, computer.trajs["lrho"]["mc_se"])
```

## BLOSUM neutrality: tolerate changes below a score threshold

```python
trajectory = snbr.Trajectory("example", sequence, seed=42)
trajectory.autofill(
    condition="evaluation", method="signature", mutational_signature="SBS2",
    measure="blosum", threshold=20, matrix_name="BLOSUM80",
    stop_policy="score",
)
print(trajectory.evaluation, trajectory.score, trajectory.accepted_steps)
```

The score distance is `alignment(original protein, original protein) - alignment(original protein, current protein)`. Global alignment, gap-open -10, gap-extension -0.5, and BLOSUM80 are the defaults. Termination occurs at distance **greater than or equal to** the positive threshold, measured against the original CDS at every step. An intermediate missense substitution may be tolerated. These accepted steps are **criterion-neutral steps**, not necessarily synonymous substitutions; do not label their mean as the original RUMS.

`stop_policy="score"` preserves score-based treatment of stop-gain/loss; `"terminate"` stops immediately when a stop position changes. This is an explicit modeling choice. BLOSUM similarity is not experimental evidence that a mutation is functionally neutral.

## Other mutation methods

```python
trajectory = snbr.Trajectory("example", sequence, seed=42)
trajectory.autofill(condition="max_length", max_length=10, method="random")

trajectory = snbr.Trajectory("example", sequence, seed=42)
trajectory.autofill(condition="nonsynonymous", method="mut_type",
                    mut_type="T[C>T]", strand="both")

trajectory = snbr.Trajectory("example", sequence, seed=42)
trajectory.autofill(condition="nonsynonymous", method="mut_types",
                    mut_types=["AA[C>T]", "AC[C>T]"], mut_type_probs=[0.6, 0.4])
```

IUPAC motifs and overlapping motif opportunities are supported. Weighted motifs and signatures use per-opportunity weights on the current sequence; their realized channel fractions depend on available sequence contexts. The algorithm permits repeated hits and back mutations. `strand_bias` splits weight between the supplied orientation and its reverse complement; default 0.5 is symmetric.

## Custom evaluations and plain CDS

```python
from sinabro import EvaluationResult

def custom(original_cds, current_cds, mutation):
    # Replace this demonstration rule with a scientifically justified criterion.
    return EvaluationResult(stop=current_cds != original_cds, reason="custom")

trajectory = snbr.Trajectory("plain", "ATGGCCTAA", seed=42,
                             cds_start=0, cds_end=9)
trajectory.autofill(condition="evaluation", method="random", evaluator=custom)
```

The callback receives CDS strings without flanks and the `MutInfo` event. It returns a bool or an `EvaluationResult`. Signature sampling needs a complete trinucleotide context and therefore cannot target the first/last base of a plain CDS lacking flanks. Use CDS-plus when all CDS positions must be eligible.

## Termination and reproducibility

`terminal_reason` distinguishes biological criterion termination, fixed length, unavailable events, and the safety cap. `accepted_steps` is `None` for incomplete or fixed-length paths. `compute_l_robustness` refuses a stopping-time mean when any path is incomplete and retains the paths in `computer.trajs["lrho"]` for inspection. It does not silently drop censored paths.

A default `max_steps=10000` bounds execution; increasing it is a scientific/run configuration choice. `max_length` is a total target trajectory length. `max_length=0` records the original state only. Use `compute_n_robustness` for fixed-n endpoint protein identity; endpoint identity includes trajectories that changed and later returned to the original phenotype.

Supply `seed` to `Trajectory` or `RobustnessComputer`, or supply a NumPy `Generator` to `Trajectory`. Mutation functions never reseed from time. Independent trajectory streams are derived from `SeedSequence`. The same seed and configuration produce identical serial/parallel paths; also record NumPy/Biopython versions. A completed stopping-time trajectory cannot be continued without removing its terminal record or constructing a fresh trajectory.

## Portable CLI and synthetic demonstration

```bash
sinabro --input examples/cdsplus.fa --output results/identity.csv \
  --signatures SBS2 SBS13 --n-trajectories 100 --seed 42
sinabro --input examples/cdsplus.fa --output results/blosum.csv \
  --signatures SBS2 SBS13 --n-trajectories 100 --seed 42 \
  --evaluation blosum --threshold 20
python examples/compare_criteria.py
```

The CLI writes per-trajectory data and a sibling `.manifest.json` with input/profile/output hashes, versions, configuration, and completion counts. IDs determine stable seed streams independently of FASTA ordering. Output files are never overwritten. These examples are synthetic package checks, not independent replications of the human/bat manuscript.

