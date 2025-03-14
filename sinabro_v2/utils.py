import numpy as np

from Bio.Data import CodonTable
from Bio.Seq import Seq, MutableSeq

def is_codon_synonymous(codon1, codon2):
    """
    Determine if two DNA codons are synonymous, i.e., they code for the same amino acid.

    Parameters:
        codon1 (str): The first codon (e.g., "ATG").
        codon2 (str): The second codon (e.g., "ATA").

    Returns:
        bool: True if both codons encode the same amino acid, False otherwise.

    Raises:
        ValueError: If either codon is invalid or not present in the standard codon table.
    """
    # Create a copy of the standard codon table mapping (without stop codons).
    codon_table = CodonTable.standard_dna_table.forward_table.copy()
    # Manually add stop codons mapping to "*" (stop signal).
    codon_table.update({"TAA": "*", "TAG": "*", "TGA": "*"})

    # Convert codons to uppercase to ensure consistency.
    codon1 = codon1.upper()
    codon2 = codon2.upper()

    try:
        # Check if both codons translate to the same amino acid.
        return codon_table[codon1] == codon_table[codon2]
    except KeyError as e:
        # Raise a clearer error if a codon is not valid.
        raise ValueError(f"Invalid codon provided: {e.args[0]}")

def generate_random_sequence(length, gc=0.5):
    """
    Generate a random DNA sequence with a specified length and GC content.

    Parameters:
    length (int): The length of the sequence to generate.
    gc (float): The desired GC content (default is 0.5), affecting nucleotide selection probability.

    Returns:
    Seq: A Bio.Seq object representing the generated DNA sequence.
    """
    # List of possible nucleotides for DNA.
    dna = ["A", "C", "G", "T"]
    
    # Generate a random sequence by selecting nucleotides according to defined probabilities.
    # The probability for 'A' and 'T' is proportional to (1 - gc) and for 'C' and 'G' is proportional to gc.
    # Multiplying by 0.5 ensures the total probability sums to 1.
    seq = Seq("".join(
        np.random.choice(
            dna, 
            length, 
            p=0.5 * np.array([1 - gc, gc, gc, 1 - gc])
        )
    ))
    
    # Return the generated sequence as a Bio.Seq object.
    return seq

def get_amino_acid_from_codon(codon, three=False):
    """
    Return the corresponding amino acid for a given codon.
    
    Parameters:
    codon (str): A DNA codon (e.g., "ATG").
    three (bool): If True, return the three-letter amino acid code; otherwise, return the one-letter code.
    
    Returns:
    str: The amino acid represented by the codon in either one-letter or three-letter format.
    """
    # Retrieve the standard DNA codon table for forward translation.
    codon_table = CodonTable.standard_dna_table.forward_table
    
    # Manually add stop codons with '*' symbol.
    codon_table["TAA"] = "*"
    codon_table["TAG"] = "*"
    codon_table["TGA"] = "*"
    
    # Mapping from one-letter amino acid codes to three-letter codes.
    one_to_three = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
        "*": "*"  # Stop codon remains as '*'
    }
    
    # If three-letter format is requested, convert the one-letter code.
    if three:
        return one_to_three[codon_table[codon]]
    else:
        return codon_table[codon]