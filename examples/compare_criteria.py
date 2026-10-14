"""Synthetic package demonstration, not a replication of manuscript results."""
from sinabro import RobustnessComputer

def main():
    seq = 'TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT'
    for (criterion, settings) in [('nonsynonymous', {}), ('blosum', {'threshold': 20})]:
        computer = RobustnessComputer('readme_example', seq, seed=42)
        value = computer.compute_l_robustness(100, 'signature', criterion, mutational_signature='SBS2', **settings)
        print(criterion, 'mean pre-terminal neutral steps:', value)
if __name__ == '__main__':
    main()
