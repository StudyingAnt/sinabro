from ..core import motif_mutation

def mutate_seq_with_mut_type(seq, mut_type, strand='both', start=1, end=None, *, rng=None, seed=None, strand_bias=0.5):
    return motif_mutation(seq, [mut_type], strand=strand, start=start, end=end, rng=rng, seed=seed, strand_bias=strand_bias)

def mutate_seq_with_mut_types(seq, mut_types, mut_type_probs, strand='both', start=1, end=None, *, rng=None, seed=None, strand_bias=0.5):
    return motif_mutation(seq, mut_types, mut_type_probs, strand=strand, start=start, end=end, rng=rng, seed=seed, strand_bias=strand_bias)
