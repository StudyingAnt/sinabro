# Sinabro - Sequential Mutations Simulator
**Sinabro** (Korean: 시나브로) means "little by little unknowingly." The program will mutate a given coding sequence one-by-one until it mets given condition, and generates a sequence of sequences which I'll call a *trajectory*.

![My Image](images/Sinabro_white.png)

## How to use
Import Sinabro as:
```python
import sinabro as snbr
```

Initialize trajectory: you need to provide id(`str`{:.python}) and data(str, Seq, MutableSeq)
```python
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
```

Display trajectory:
```python
traj.show()


```

Generate trajectory: you can generate trajectory automatically by giving condition(max_length or )`x = 4`{:.ruby}