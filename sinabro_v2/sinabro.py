from typing import Union
from Bio.Seq import Seq, MutableSeq

from .types.types import MutInfo, MutationRecord
from . import mutate
from . import evaluate

class Trajectory:
    """
    Class to represent a trajectory of mutations.
    
    Attributes:
        _traj_id: Identifier for the trajectory.
        _original_seq (Seq): The original nucleotide sequence.
        _data (list): List of sequences (Seq objects) representing the trajectory.
        _seq_length (int): Length of the original sequence.
        _length (int): The number of mutations added.
    """
    def __init__(self, traj_id, record: MutationRecord, note="."):
        """
        Initialize the Trajectory with an initial MutationRecord.
        
        Args:
            traj_id: Identifier for the trajectory.
            record (MutationRecord): The initial mutation record.
            note (str): Additional note (unused in this context).
        
        Raises:
            ValueError: If the sequence in the record is None.
            TypeError: If the sequence type is invalid.
        """
        # Validate and convert the input mutation record
        valid_record = self._validate_and_convert_record(record)
        
        # Save the original sequence from the validated record
        self._original_seq = valid_record.sequence
        
        # Initialize trajectory attributes
        self._traj_id = traj_id
        self._records = [valid_record]  # List to store sequences over the trajectory
        self._seq_length = len(self._original_seq)
        self._length = 0

    def __str__(self) -> str:
        """
        Return a formatted string representing the trajectory details.
        """
        output = []
        output.append(f"Trajectory ID: {self._traj_id}")
        output.append(f"Original Sequence: {self._original_seq}")
        output.append(f"Sequence Length: {self._seq_length}")
        output.append(f"Number of Mutations: {self._length}")
        output.append("Mutation Steps:")
        for i, record in enumerate(self._records):
            output.append(f"  Step {i}: {record.sequence}")
        return "\n".join(output)

    def show(self, verbose: bool = True):
        """
        Display the trajectory details.
        
        Args:
            verbose (bool): If True, display detailed mutation steps.
        """
        # Print the trajectory header information
        print(f"Trajectory ID: {self._traj_id}")
        print(f"Original Sequence: {self._original_seq}")
        print(f"Sequence Length: {self._seq_length}")
        print(f"Number of Mutations: {self._length}")

        # If verbose is True, print each mutation step in detail
        if verbose:
            print("\nMutation Steps:")
            for i, record in enumerate(self._records):
                # Here you could add more detailed info for each step if available
                print(f"  Step {i}: {record.sequence}")

    def add_record(self, record: MutationRecord) -> list:
        """
        Append a new mutation record to the trajectory.
        
        Args:
            record (MutationRecord): A mutation record object containing the mutation details.
        
        Returns:
            list: The updated list of MutationRecord objects.
            
        Raises:
            ValueError: If record.sequence is None.
            TypeError: If record.sequence is not a valid type.
        """
        if record.sequence is None:
            raise ValueError("Sequence in the mutation record must not be None")
        
        # Validate and convert the mutation record.
        valid_record = self._validate_and_convert_record(record)
        self._records.append(valid_record)
        self._length += 1
        return self._records

    def remove_last_record(self) -> MutationRecord:
        """
        Remove the last mutation record from the trajectory.
        
        Returns:
            MutationRecord: The removed mutation record.
        
        Raises:
            IndexError: If attempting to remove the initial record.
        """
        if len(self._records) <= 1:
            raise IndexError("Cannot remove the initial record.")
        removed_record = self._records.pop()
        self._length -= 1
        return removed_record

    
    """
    Autofill related methods
    """
    def add_mutated_sequence(self, method, **kwargs):
        """
        Append a mutated sequence to the trajectory based on the specified mutation method.
        
        Parameters:
            **kwargs: Keyword arguments containing mutation parameters.
                Expected keys include:
                    - method (str): The mutation method to apply ("random", "mut_type", "mut_types", "signature").
                    - note (str, optional): Additional note for the mutation.
                    - ... (other method-specific parameters)
        
        Returns:
            list: The updated list of sequence records.
        
        Raises:
            ValueError: If an invalid mutation method is provided.
        """
        # Mapping of mutation method names to their corresponding functions.
        mutation_methods = {
            "random": mutate.random_single_substitution,
            "mut_type": mutate.mutate_seq_with_mut_type,
            #"mut_types": self._apply_mut_types_mutation,
            "signature": mutate.mutate_seq_with_mutational_signature
        }
        
        # Retrieve the mutation method name from keyword arguments.
        method_name = method
        if method_name not in mutation_methods:
            raise ValueError("Invalid mutation method")
        
        # Execute the corresponding mutation function with the provided parameters.
        seq = self._records[-1].sequence
        mutinfo = mutation_methods[method_name](seq, **kwargs)
        
        # Check if the mutation process was successful (error code is 0).
        if not mutinfo.e:
            # Append the new mutated sequence to the trajectory records.
            record = MutationRecord(
                mutinfo.new_seq,
                hgvs_mrna=mutinfo.hgvs_mrna,
                mut_type=mutinfo.mut_type,
                note=kwargs.get("note", ".")
            )
            self.add_record(record)
        else:
            # Add error handling logic here if necessary.
            pass
    
        return self._records

    def autofill(self, method, eval_method, **kwargs):
        eval_methods = {
            "max_length": evaluate.eval_maxlen,
            "nonsynonymous": ,
            "blosum":
        }

        if eval_method not in eval_methods:
            raise ValueError("Invalid evaluation method")

        run_flag = True
        while run_flag:
            eval_methods[eval_method](self._records)


        print("autofill")
        

    """
    Internal helper methods
    """
    def _validate_and_convert_record(self, record: MutationRecord) -> MutationRecord:
        """
        Validate and convert the sequence in the mutation record.
        
        Args:
            record (MutationRecord): A mutation record to validate.
        
        Returns:
            MutationRecord: The validated and converted mutation record.
            
        Raises:
            ValueError: If the sequence is None.
            TypeError: If the sequence is not a string, Seq, or MutableSeq.
        """
        seq = record.sequence
        if seq is None:
            raise ValueError("data must not be None")
        # Convert the sequence to a Seq object if it's a string or MutableSeq
        if isinstance(seq, (str, MutableSeq)):
            record.sequence = Seq(seq)
            return record
        elif isinstance(seq, Seq):
            return record
        else:
            raise TypeError("data should be a string, Seq, or MutableSeq object")