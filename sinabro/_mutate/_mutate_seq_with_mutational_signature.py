"""Compatibility names for validated SBS profile primitives."""
from Bio.Seq import Seq
from ..core import signature_mutation, full_signature, probability_matrix, BASES, TRIPLETS

dna_bases=list(BASES)
trinucleotides=[Seq(t) for t in TRIPLETS]

def make_full_mutational_signature(input_file_path,column=None,strand_bias=0.5):
    return full_signature(input_file_path,column,strand_bias)

def compute_mut_prob_matrix(seq,full_mutational_signature,start=1,end=None):
    return probability_matrix(seq,full_mutational_signature,start,end)

def mutate_seq_with_mutational_signature(seq,start=None,end=None,mutational_signature='SBS1',
    cosmic_version=3.3,genome_ref='GRCh37',custom_signature_path=None,column=None,strand_bias=0.5,
    *,rng=None,seed=None):
    return signature_mutation(seq,start=1 if start is None else start,end=end,
        mutational_signature=mutational_signature,cosmic_version=cosmic_version,
        genome_ref=genome_ref,custom_signature_path=custom_signature_path,column=column,
        strand_bias=strand_bias,rng=rng,seed=seed)
