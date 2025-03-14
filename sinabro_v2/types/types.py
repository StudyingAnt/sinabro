from dataclasses import dataclass
from typing import Union
from Bio.Seq import Seq, MutableSeq

@dataclass
class MutInfo:
    """
    Data class to store mutation information.
    
    Attributes:
        new_seq (Union[str, Seq, MutableSeq]): The mutated sequence.
        idx_target (int): The index of the mutated nucleotide.
        hgvs_mrna (str): HGVS notation for the mutation.
        mut_type (str): Mutation type description.
        e (int): Error code (0 for success).
    """
    new_seq: Union[str, Seq, MutableSeq]
    idx_target: int
    hgvs_mrna: str
    mut_type: str
    e: int

@dataclass
class MutationRecord:
    """
    Data class to store mutation record information.
    
    Attributes:
        sequence (Seq): The nucleotide sequence.
        hgvs_mrna (str): HGVS notation for mRNA.
        hgvs_aa (str): HGVS notation for amino acids.
        mut_type (str): Mutation type.
        note (str): Additional note.
    """
    sequence: Seq
    hgvs_mrna: str = "."
    hgvs_aa: str = "."
    mut_type: str = "."
    note: str = "."