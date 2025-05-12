from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord

from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

def compute_alignment(seq1, seq2, mode="global", gap_open=-10, gap_extend=-0.5, matrix_name="BLOSUM80"):
    """
    Perform an alignment between two amino acid sequences.
    
    Parameters:
        seq1 (str): The first amino acid sequence.
        seq2 (str): The second amino acid sequence.
        mode (str): Alignment mode ('global', 'local', 'semiglobal', etc.). Default is "global".
        gap_open (float): Gap opening penalty. Default is -10.
        gap_extend (float): Gap extension penalty. Default is -0.5.
        matrix_name (str): Name of the substitution matrix to use. Default is "BLOSUM62".
        
    Returns:
        best_alignment: The alignment result with the highest score.
    """
    # Create and configure the PairwiseAligner object
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend

    # Load and apply the substitution matrix
    substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.substitution_matrix = substitution_matrix

    # Perform the alignment and return the best (highest scoring) alignment
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    return best_alignment


def compute_align_distance(seq1, seq2):
    """
    Compute the distance between two sequences based on their alignment scores.
    
    This function calculates the maximum possible alignment score for seq1 (self-alignment)
    and the best alignment score between seq1 and seq2. The distance is defined as the
    difference between the maximum self-alignment score and the alignment score of seq1 and seq2.
    
    Parameters:
        seq1 (str): The first amino acid sequence.
        seq2 (str): The second amino acid sequence.
        
    Returns:
        float: The computed distance.
    """
    # Compute the maximum possible alignment score (self-alignment of seq1)
    self_alignment = compute_alignment(seq1, seq1)
    max_score = self_alignment.score

    # Compute the best alignment score between seq1 and seq2
    best_alignment = compute_alignment(seq1, seq2)
    best_score = best_alignment.score

    # Calculate the distance as the difference between the maximum score and the best alignment score
    distance = max_score - best_score
    return distance

def eval_blosum(records, mutinfo, **kwargs):
    threshold = kwargs.get('threshold', None)
    if threshold is None:
        raise ValueError("threshold must be provided")

    ori_seq = records[0].sequence[1:-1]
    new_seq = records[-1].sequence[1:-1]

    dist = compute_align_distance(ori_seq, new_seq)

    if dist < threshold:
        return False, dist
    else:
        return True, dist