from sinabro._mutate._mutate_seq_with_mutational_signature import make_full_mutational_signature, compute_mut_prob_matrix, mutate_seq_with_mutational_signature, dna_bases, trinucleotides
from sinabro.core import dna, TRIPLET_INDEX
trint_dict = TRIPLET_INDEX

def get_trinucleotides(seq):
    text = dna(seq)
    return [text[i - 1:i + 2] for i in range(1, len(text) - 1)]
