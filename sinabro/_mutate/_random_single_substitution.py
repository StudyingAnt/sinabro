from ..core import random_mutation

def random_single_substitution(seq, start=1, end=None, *, rng=None, seed=None):
    return random_mutation(seq, start=start, end=end, rng=rng, seed=seed)
