# Sinabro - Sequential Mutations Simulator
**Sinabro** (Korean: 시나브로) means "little by little unknowingly." The program will mutate a given coding sequence one-by-one until it mets given condition, and generates a sequence of sequences which I'll call a *trajectory*.

![My Image](images/Sinabro_white.png)

## How to use
### 1. Import Sinabro as:
```python
import sinabro as snbr
```

### 2. Initialize trajectory: 
You need to provide id (`str`) and data (`str`, `Seq`, or `MutableSeq`)
```python
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
```

### 3. Display trajectory:
```python
traj.show()
```
__result__
```console
ID: test
Trajectory length: 0

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
```

### 4. Generate trajectory: 
You can generate trajectory automatically by giving a condition (`max_length` or `nonsynonymous`) and a method (`random`, `mut_type`, or `signature`)
__Example 1__
```python
# max_length with random
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="random")

traj.show()
```
__result__
```console
ID: Test
Trajectory length: 10

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCCGACTCGGTCATCGATCGGTTCTCATTGAT   c.5T>C      p.2Leu>Pro     [T>C]           .
TATGCCGACTCGGTCATCGTTCGGTTCTCATTGAT   c.19A>T     p.7Ile>Phe     [A>T]           .
TATGCAGACTCGGTCATCGTTCGGTTCTCATTGAT   c.5C>A      p.2Pro>Gln     [C>A]           .
TATGCAGACTCGGGCATCGTTCGGTTCTCATTGAT   c.13T>G     p.5Ser>Ala     [T>G]           .
TATGCAGACTCGGGCCTCGTTCGGTTCTCATTGAT   c.15A>C     .              [A>C]           .
TATGCACACTCGGGCCTCGTTCGGTTCTCATTGAT   c.6G>C      p.2Gln>His     [G>C]           .
TATGCACACTCGGGCCTCGTTCTGTTCTCATTGAT   c.22G>T     p.8Gly>Cys     [G>T]           .
TATGCACACTCGGGCCTCGTTCTTTTCTCATTGAT   c.23G>T     p.8Cys>Phe     [G>T]           .
TATGCACACTCGGGCCTCCTTCTTTTCTCATTGAT   c.18G>C     .              [G>C]           .
TATGCACACTCGGGCTTCCTTCTTTTCTCATTGAT   c.15C>T     .              [C>T]           .
```
__Example 2__
```python
# max_length with mut_type
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="mut_type", mut_type="T[C>T]")
```
__result__
```console
ID: Test
Trajectory length: 6

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCTGACTTGGTCATCGATCGGTTCTCATTGAT   c.10C>T     p.4Arg>Trp     T[C>T]          .
TATGCTGACTTGGTCATCGATCGGTTTTCATTGAT   c.26C>T     p.9Ser>Phe     T[C>T]          .
TATGCTGACTTGGTTATCGATCGGTTTTCATTGAT   c.14C>T     p.5Ser>Leu     T[C>T]          .
TATGCTGACTTGGTTATCGATTGGTTTTCATTGAT   c.21C>T     .              T[C>T]          .
TATGCTGACTTGGTTATCGATTGGTTTTTATTGAT   c.28C>T     p.10His>Tyr    T[C>T]          .
TATGCTGACTTGGTTATTGATTGGTTTTTATTGAT   c.17C>T     p.6Ser>Leu     T[C>T]          .
```
__Example 3__
```python
# max_length with signature
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="signature", mutational_signature="SBS2")
```
__result__
```console
ID: Test
Trajectory length: 10

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCTGACTCGGTCATCGATCGGTTCTCATTAAT   c.32G>A     .              T[C>T]A         .
TATGCTGACTCGGTCATCGATCGGTTTTCATTAAT   c.26C>T     p.9Ser>Phe     T[C>T]T         .
TATGCTAACTCGGTCATCGATCGGTTTTCATTAAT   c.6G>A      .              T[C>T]A         .
TATGCTAACTCGGTTATCGATCGGTTTTCATTAAT   c.14C>T     p.5Ser>Leu     T[C>T]A         .
TATGCTAACTCGGTTATCAATCGGTTTTCATTAAT   c.18G>A     .              T[C>T]G         .
TATGCTAACTCGGTTATTAATCGGTTTTCATTAAT   c.17C>T     p.6Ser>Leu     T[C>T]A         .
TATGCTAACTTGGTTATTAATCGGTTTTCATTAAT   c.10C>T     p.4Arg>Trp     T[C>T]G         .
TATGCTAACTTGGTTATTAATCGGTTTTTATTAAT   c.28C>T     p.10His>Tyr    T[C>T]A         .
TATGCTAACTTGGTTATTAATCAGTTTTTATTAAT   c.22G>A     p.8Gly>Ser     C[C>T]G         .
TATGCTAACTTGGTTATTAATTAGTTTTTATTAAT   c.21C>T     .              T[C>T]A         .
```
__Example 4__
```python
# nonsynonymous with random
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="random")
```
__result__
```console
ID: Test
Trajectory length: 1

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCTGACTCGGTTATCGATCGGTTCTCATTGAT   c.14C>T     p.5Ser>Leu     [C>T]           .
```
__Example 5__
```python
# nonsynonymous with mut_type
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="mut_type", mut_type="T[C>T]")
```
__result__
```console
ID: Test
Trajectory length: 1

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCTGACTTGGTCATCGATCGGTTCTCATTGAT   c.10C>T     p.4Arg>Trp     T[C>T]          .
```
__Example 6__
```python
# nonsynonymous with signature
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="signature", mutational_signature="SBS2")
```
__result__
```console
ID: Test
Trajectory length: 1

Sequence                              HGVS_mRNA   HGVS_Protein   Mutation_Type   Note
TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT   .           .              .               .
TATGCTGACTCGGTCATCGATCGGTTCTTATTGAT   c.28C>T     p.10His>Tyr    T[C>T]A         .
```

