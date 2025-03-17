import os
import time
import math

import numpy as np
import pandas as pd

from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord
from ._helper import preserve_seq_type

# Generate a list of all possible trinucleotides as BioPython Seq objects using list comprehension
dna_bases = ["A", "C", "G", "T"]
trinucleotides = [Seq(a + b + c) for a in dna_bases for b in dna_bases for c in dna_bases]

# Create a dictionary for fast lookup of trinucleotide indices
trint_dict = {str(trinuc): i for i, trinuc in enumerate(trinucleotides)}

def make_full_mutational_signature(
    input_file_path: str,
    column: str = None,
    strand_bias: float = 0.5
) -> np.ndarray:
    """
    Generate a full mutational signature matrix from an SBS (single base substitution) file.
    
    This function reads an SBS file where each mutation type is indexed, and constructs a 4x64 
    matrix representing the mutational signature based on trinucleotide contexts. For each mutation type,
    the function extracts the forward trinucleotide context and its reverse complement, then assigns
    the mutation percentage (adjusted by strand bias) to the corresponding positions in the matrix.
    
    Parameters:
        input_file_path (str): The file path to the SBS file (tab-separated, with headers and index).
        column (str, optional): The column name from the SBS file to use. If None, the first column is used.
        strand_bias (float, optional): The factor to apply to the mutation percentage (default is 0.5).
        
    Returns:
        np.ndarray: A 2D numpy array (4x64) representing the full mutational signature.
    """
    # Read the SBS file into a pandas DataFrame
    sbs = pd.read_csv(input_file_path, sep="\t", index_col=0, header=0)

    # Select the column to use for mutational percentages; default to the first column if not provided
    if column is None:
        column = sbs.columns[0]

    # Initialize a 4x64 numpy array with zeros for the mutational signature matrix
    full_mutational_signature = np.zeros((4, 64))
    
    # Get the list of mutation types from the DataFrame index
    mut_types = list(sbs.index)
    
    for mut_type in mut_types:
        # Extract the trinucleotide context for the forward strand based on specific indices
        trinucleotide_f = f"{mut_type[0]}{mut_type[2]}{mut_type[6]}"
        # Calculate the reverse complement of the forward trinucleotide context
        trinucleotide_r = str(Seq(trinucleotide_f).reverse_complement())
        
        # Get the mutated nucleotide for the forward mutation
        changed_nt_f = mut_type[4]
        # Calculate the complement of the mutated nucleotide for the reverse mutation
        changed_nt_r = str(Seq(changed_nt_f).complement())

        # Calculate the mutation percentage adjusted by strand bias
        mut_type_percentage = strand_bias * sbs.loc[mut_type, column]

        # Determine the index for the mutated nucleotide in the list of DNA bases
        i1 = dna_bases.index(changed_nt_f)
        # Determine the index for the forward trinucleotide in the list of trinucleotides
        j1 = trinucleotides.index(Seq(trinucleotide_f))
        # Determine the index for the complement nucleotide in the list of DNA bases
        i2 = dna_bases.index(changed_nt_r)
        # Determine the index for the reverse complement trinucleotide in the list of trinucleotides
        j2 = trinucleotides.index(Seq(trinucleotide_r))
        
        # Assign the adjusted mutation percentage to both forward and reverse contexts in the matrix
        full_mutational_signature[i1, j1] = mut_type_percentage
        full_mutational_signature[i2, j2] = mut_type_percentage

    return full_mutational_signature

def get_trinucleotides(seq: str) -> list:
    """
    Extract trinucleotides from a sequence using list comprehension.
    
    Parameters:
        seq (str): Input nucleotide sequence.
        
    Returns:
        List[str]: List of trinucleotides (each of length 3).
    """
    return [seq[i - 1: i + 2] for i in range(1, len(seq) - 1)]

def compute_mut_prob_matrix(seq, full_mutational_signature):
    """
    Compute the mutation probability matrix for a given sequence based on a full mutational signature.
    
    This function constructs a trinucleotide context vector from the input sequence and maps each trinucleotide
    (for positions with both a left and right neighbor) to its corresponding index in the global 'trinucleotides' list.
    It then builds a probability matrix (of shape 4 x len(seq)) using the provided full_mutational_signature matrix.
    The matrix is normalized by the total sum of its elements.
    
    Parameters:
        seq (str): The input nucleotide sequence.
        full_mutational_signature (np.ndarray): A 2D numpy array representing the full mutational signature (4 x 64).
        
    Returns:
        np.ndarray: A normalized mutation probability matrix with shape (4, len(seq)).
    """
    # Extract trinucleotides using list comprehension (ignoring the first and last nucleotide)
    trint_list = [str(seq[i - 1: i + 2]) for i in range(1, len(seq) - 1)]
    
    # Create the trinucleotide vector with placeholders for the first and last positions
    trint_vector = ["-"] + trint_list + ["-"]
    
    # Create the trinucleotide type vector using the precomputed dictionary,
    # with -1 as placeholder values for positions without a trinucleotide context
    trint_type_vector = [-1] + [trint_dict[trint] for trint in trint_list] + [-1]
    
    # Initialize the probability matrix (4 rows for nucleotides, columns equal to sequence length)
    P_before_norm = np.zeros((4, len(seq)))
    
    # Vectorized assignment for positions where a trinucleotide context exists (positions 1 to len(seq)-2)
    valid_indices = np.array(trint_type_vector[1:-1])
    P_before_norm[:, 1:-1] = full_mutational_signature[:, valid_indices]
    
    # Normalize the matrix by the total sum of its elements
    S = np.sum(P_before_norm)
    P = P_before_norm / S
    
    return P

@preserve_seq_type
def mutate_seq_with_mutational_signature(seq, **kwargs):
    """
    Mutate a nucleotide sequence based on a specified mutational signature.
    
    This function mutates the input sequence using a full mutational signature matrix derived from either a custom file
    or pre-defined COSMIC signatures. It computes a mutation probability matrix for the sequence, selects a mutation 
    site probabilistically, constructs HGVS notation for the mutation, and applies the mutation to the sequence.
    
    Parameters:
        seq (Union[str, Seq, MutableSeq]): The input nucleotide sequence.
        mutational_signature (str): Mutational signature to use ('SBS1', 'SBS2', ... or 'custom').
        cosmic_version (float): COSMIC version number to use for pre-defined signatures.
        genome_ref (str): Genome reference to use for pre-defined signatures.
        custom_signature_path (str, optional): File path for custom signature if 'custom' is selected.
        column (str, optional): Column name in the signature file to use.
        strand_bias (float): Factor to apply for strand bias in mutational signature.
        
    Returns:
        tuple[Union[str, Seq, MutableSeq], int, str, int]: A tuple containing the mutated sequence, target index, 
        HGVS notation, mutation type, and error code (0 for success).
    """
    # Set default parameter values
    defaults = {
        'mutational_signature': "SBS1",
        'cosmic_version': 3.3,
        'genome_ref': "GRCh37",
        'custom_signature_path': None,
        'column': None,
        'strand_bias': 0.5
    }
    # Update defaults with values provided in kwargs (if any)
    defaults.update(kwargs)
    
    # Extract parameters from the defaults dictionary
    mutational_signature = defaults['mutational_signature']
    cosmic_version = defaults['cosmic_version']
    genome_ref = defaults['genome_ref']
    custom_signature_path = defaults['custom_signature_path']
    column = defaults['column']
    strand_bias = defaults['strand_bias']

    start = 0
    end = len(seq) - 1

    # Define available mutational signatures
    signatures = [
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
        "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
        "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
        "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
        "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
        "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
        "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
        "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
        "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
        "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
        "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
        "SBS95"
    ]

    # Load full mutational signature matrix
    if mutational_signature == "custom":
        if custom_signature_path is None:
            raise ValueError("custom_signature_path must be provided")
        else:
            full_mutational_signature = make_full_mutational_signature(
                input_file_path=custom_signature_path,
                column=column,
                strand_bias=strand_bias
            )
    elif mutational_signature in signatures:
        dir_name = os.path.dirname(__file__)
        tokens = [
            "../../cosmic_signatures/COSMIC_v",
            str(cosmic_version), ".1_SBS_",
            genome_ref, ".txt"
        ]
        signature_file = "".join(tokens)
        signature_file_path = os.path.join(dir_name, signature_file)
        full_mutational_signature = make_full_mutational_signature(
            input_file_path=signature_file_path,
            column=mutational_signature,
            strand_bias=strand_bias
        )
    else:
        raise ValueError("signature must be SBSXX or custom")

    # Compute mutation probability matrix for the sequence
    P = compute_mut_prob_matrix(seq, full_mutational_signature)

    # Check if there are any valid mutation probabilities
    num_nonzero = np.count_nonzero(P)
    if not num_nonzero:
        return MutInfo(MutableSeq(""), -1, "", "", 1)
    
    # Seed the random number generator using current time
    seed1 = time.time()
    seed2 = (seed1 - int(seed1)) * (10**7)
    np.random.seed(int(seed1 + seed2))
    
    # Randomly select a mutation site based on the probability matrix P (flattened in Fortran order)
    idx_flat = np.random.choice(len(seq) * 4, p=P.flatten("F"))
    idx_base = idx_flat % 4
    idx_target = idx_flat // 4  # Using integer division instead of math.floor

    # Note: 'dna_bases' should be defined globally (e.g., dna_bases = ["A", "C", "G", "T"])
    # Construct HGVS notation and mutation type based on the target nucleotide
    if str(seq[idx_target]) not in ["C", "T"]:
        # For nucleotides other than C or T, use complementary representation
        hgvs_mrna_tokens = [
            "c.", str(idx_target),
            str(seq[idx_target]), ">",
            dna_bases[idx_base]
        ]
        mut_type_tokens = [
            str(Seq(seq[idx_target+1]).complement()), "[", 
            str(Seq(seq[idx_target]).complement()), ">", 
            str(Seq(dna_bases[idx_base]).complement()), "]",
            str(Seq(seq[idx_target-1]).complement())
        ]
        hgvs_mrna = "".join(hgvs_mrna_tokens)
        mut_type = "".join(mut_type_tokens)
    else:
        # For C or T, use the original bases
        hgvs_mrna_tokens = [
            "c.", str(idx_target),
            str(seq[idx_target]), ">",
            dna_bases[idx_base]
        ]
        mut_type_tokens = [
            str(seq[idx_target-1]), "[",
            str(seq[idx_target]), ">",
            dna_bases[idx_base], "]",
            str(seq[idx_target+1])    
        ]
        hgvs_mrna = "".join(hgvs_mrna_tokens)
        mut_type = "".join(mut_type_tokens)
        
    # Apply the mutation to the sequence
    seq[idx_target] = dna_bases[idx_base]

    # Always work with MutableSeq type inside the function.
    return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)