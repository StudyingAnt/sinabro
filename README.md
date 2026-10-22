# Sinabro

![Sinabro](docs/assets/sinabro.png)

**Sequential mutation simulation for coding-sequence robustness.**

Sinabro applies one nucleotide substitution at a time, evaluates the translated coding sequence, and records the mutation trajectory. It supports random substitutions, motif preferences and COSMIC SBS signatures, with protein-identity and BLOSUM80 evaluation criteria.

This is the **Sinabro 3.0.0** package accompanying *Human and bats genome robustness under COSMIC mutational signatures*. The separate [manuscript analysis repository](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds) contains reference downloads, human/bat simulations, statistical analyses and figure workflows.

## Installation

Python 3.9 or later, NumPy and Biopython are required. Install from this repository:

```bash
git clone https://github.com/StudyingAnt/sinabro.git
cd sinabro
python -m pip install .
```

For the recorded scientific-library environment:

```bash
conda env create --file environment.yml
conda activate sinabro
python -m pip install --no-deps --no-build-isolation .
```

The Conda specification uses Python 3.9.15, NumPy 1.21.6, Biopython 1.79 and SciPy 1.8.1. Conda supplies the dependencies and pip installs the package. For manuscript reproduction, use the [analysis environment and setup guide](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/blob/main/CONDA_SETUP.md), which also pins plotting dependencies.

The optional `analysis` extra supplies SciPy for pairwise dN/dS calculations: `python -m pip install ".[analysis]"`. Editable installation is available for package development; a regular installation is recommended when the checkout shares a parent directory with an analysis repository.

## Quick start

```python
from sinabro import Trajectory, RobustnessComputer

# One context nucleotide on either side of the CDS.
sequence = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
trajectory = Trajectory("example", sequence, seed=42)
trajectory.autofill(
    condition="nonsynonymous",
    method="signature",
    mutational_signature="SBS2",
)
print(trajectory.accepted_steps)

computer = RobustnessComputer("example", sequence, seed=42)
rums = computer.compute_l_robustness(
    n_sim=100,
    method="signature",
    eval_method="nonsynonymous",
    mutational_signature="SBS2",
)
print(rums)
```

Under protein identity, `accepted_steps` counts neutral substitutions before the first protein-sequence change; its mean across complete trajectories is RUMS. The default CDS is `sequence[1:-1]`. Plain CDS inputs require explicit bounds. BLOSUM80 supplies a relaxed sequence criterion and must be labelled separately from synonymous robustness. These criteria are not experimental measurements of fitness or function.

## Command-line use

```bash
sinabro --input examples/cdsplus.fa --output results/identity.csv \
  --signatures SBS2 SBS13 --n-trajectories 100 --seed 42
python examples/compare_criteria.py
```

The CLI writes trajectory results and a manifest with input hashes, settings, dependency versions and completion counts. Example sequences are synthetic. [API and CLI documentation](docs/api.md) describes BLOSUM settings, coding coordinates and custom evaluation policies.

## Documentation and tests

| Resource | Contents |
| --- | --- |
| [API guide](docs/api.md) | Python and command-line examples |
| [Algorithm](docs/algorithm.md) | Mutation opportunities, evaluation and stopping rules |
| [Scientific provenance](docs/reproducibility.md) | Recorded environment, validated analysis scope and historical source |
| [Testing](docs/testing.md) | Package tests, wheel checks and CI configuration |
| [Change log](CHANGELOG.md) | Version history and interface changes |

```bash
python -m unittest discover -s tests -v
```

The 89 package tests cover coding coordinates, all 64 codons and nine possible substitutions per codon, mutation probabilities, stopping rules, metrics, reproducibility and CLI output. Package behavior tests and the full human/bat analyses have separate validation records. GitHub Actions defines Linux/Windows checks; consult the run for the relevant commit for its actual CI status.

The installable code is under `sinabro/`, packaged SBS tables under `sinabro/data/`, executable examples under `examples/`, and tests under `tests/`. Historical interfaces are documented in [migration notes](docs/migration.md); previous implementations remain recoverable in Git history.

## Citation and license

Use [CITATION.cff](CITATION.cff) for software attribution and cite the associated manuscript when using its analyses. Record the package commit and analysis commit with each scientific run. The current edition is distinguished from the preserved historical implementation in the provenance guide.

BSD-3-Clause. Copyright Joon-Hyun Song. See [LICENSE](LICENSE).
