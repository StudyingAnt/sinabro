import numpy as np
from typing import Union
from Bio.Seq import Seq, MutableSeq

from multiprocessing import Pool
import os
import time

from .types.types import MutInfo, MutationRecord
from . import mutate
from . import evaluate

from .utils import get_codon, get_amino_acid_from_codon
from .evaluate import compute_alignment, compute_align_distance

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
            output.append(f"  Step {i}: {record.sequence}\t{record.hgvs_mrna}\t{record.hgvs_aa}\t{record.mut_type}\t{record.note}")
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
                print(f"  Step {i}: {record.sequence}\t{record.hgvs_mrna}\t{record.hgvs_aa}\t{record.mut_type}\t{record.note}")

    def get_records(self):
        return self._records
    
    def get_length(self):
        return self._length

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

        old_aa = get_amino_acid_from_codon(get_codon(seq, mutinfo.idx_target))
        new_aa = get_amino_acid_from_codon(get_codon(mutinfo.new_seq, mutinfo.idx_target))

        if old_aa != new_aa:
            pos = mutinfo.idx_target//3+1
            hgvs_aa = f"p.{old_aa}{pos}{new_aa}"
        else:
            hgvs_aa = f"."

        
        # Check if the mutation process was successful (error code is 0).
        if not mutinfo.e:
            # Append the new mutated sequence to the trajectory records.
            record = MutationRecord(
                mutinfo.new_seq,
                hgvs_mrna=mutinfo.hgvs_mrna,
                hgvs_aa=hgvs_aa,
                mut_type=mutinfo.mut_type,
                note=kwargs.get("note", ".")
            )
            self.add_record(record)
        else:
            # Add error handling logic here if necessary.
            pass
    
        return mutinfo

    def autofill(self, method, eval_method, **kwargs):
        eval_methods = {
            "max_length": evaluate.eval_maxlen,
            "nonsynonymous": evaluate.eval_nonsym,
            "blosum": evaluate.eval_blosum
        }

        if eval_method not in eval_methods:
            raise ValueError("Invalid evaluation method")

        iter_num = 0
        run_flag = True
        while run_flag:
            # Add mutation
            mutinfo = self.add_mutated_sequence(method, **kwargs)
            stop_flag, score = eval_methods[eval_method](self._records, mutinfo, **kwargs)

            if stop_flag:
                run_flag = False
            
            if eval_method == 'blosum' and iter_num > kwargs.get('max_iter', 20):
                run_flag = False

            if score is not None:
                self._records[-1].note = f"score: {score}"

            iter_num += 1

        

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
        
# Module-level helper function for generating a single trajectory.
def _gen_trajectory_helper(args):
    self, traj_id, method, eval_method, kwargs = args
    # Set a different seed for each process to ensure randomness.
    seed = (int(time.time() * 1000) + traj_id + os.getpid()) % (2**32)
    np.random.seed(seed)
    return self.generate_trajectory(traj_id, method, eval_method, **kwargs)


def _traj_contribution(traj, phi_eval, by_score, threshold):
    records = traj.get_records()
    orig_seq = records[0].sequence[1:-1]
    last_seq = records[-1].sequence[1:-1]

    if phi_eval == 'nonsym':
        return 1 if traj._compare_protein_sequences(orig_seq, last_seq) else 0
    elif phi_eval == 'blosum':
        if by_score:
            # alignment = compute_alignment(orig_seq, last_seq)
            dist = compute_align_distance(orig_seq, last_seq)
            return dist
            # return alignment.score
        else:
            if threshold is None:
                raise ValueError("threshold must be provided for blosum distance")
            dist = compute_align_distance(orig_seq, last_seq)
            return 1 if dist <= threshold else 0
    else:
        raise ValueError(f"Unknown phi_eval: {phi_eval}")

class RobustnessComputer:
    def __init__(self, gene_name, gene_seq):
        self.gene_name = gene_name
        self.gene_seq = gene_seq
        self.trajs = {}

    def generate_trajectory(self, traj_id, method, eval_method, **kwargs):
        record = MutationRecord(sequence=self.gene_seq)
        traj = Trajectory(traj_id, record)

        traj.autofill(method, eval_method, **kwargs)

        return traj
    
    def generate_trajectories(self, n_traj, method, eval_method, **kwargs):
        # If 'multiprocessing' is included in kwargs, use its value and remove it.
        use_multiprocessing = kwargs.pop("multiprocessing", False)
        
        if use_multiprocessing:
            # Determine the number of CPUs to use: default is available CPUs minus 2 (minimum 1)
            default_n_cpu = max(1, (os.cpu_count() or 1) - 2)
            n_cpu = kwargs.pop("n_cpu", default_n_cpu)

            # Prepare argument list for each trajectory.
            arg_list = [(self, traj_id, method, eval_method, kwargs) for traj_id in range(n_traj)]
            
            with Pool(n_cpu) as pool:
                trajs = pool.map(_gen_trajectory_helper, arg_list)
            return trajs
        else:
            trajs = []
            for traj_id in range(n_traj):
                traj = self.generate_trajectory(traj_id, method, eval_method, **kwargs)
                trajs.append(traj)
            return trajs

    # def generate_trajectories(self, n_traj, method, eval_method, **kwargs):
    #     trajs = []
    #     for traj_id in range(n_traj):
    #         traj = self.generate_trajectory(traj_id, method, eval_method, **kwargs)
    #         trajs.append(traj)

    #     return trajs
    
    def compute_l_robustness(self, n_sim, method, eval_method, **kwargs):
        trajs = self.generate_trajectories(n_sim, 
                                           method=method, 
                                           eval_method=eval_method, 
                                           **kwargs)
        
        l = []
        for traj in trajs:
            l.append(traj.get_length()-1)

        self.trajs["lrho"] = {
            "robustness": np.array(l).mean(),
            "trajs": trajs,
            "parameters": {
                "n_sim": n_sim,
                "method": method,
                "eval_method": eval_method,
                **kwargs
            }
        }

        return np.array(l).mean()
    
    def _compare_protein_sequences(self, seq1, seq2):
        """
        Compare the translated amino acid sequences of two DNA sequences.
        
        Parameters:
            seq1 (str): The first DNA sequence.
            seq2 (str): The second DNA sequence.
        
        Returns:
            bool: True if both translated protein sequences are identical, False otherwise.
        """
        # Translate both DNA sequences to protein sequences
        protein_seq1 = Seq(seq1).translate()
        protein_seq2 = Seq(seq2).translate()
    
        # Return the result of comparing the two protein sequences
        return protein_seq1 == protein_seq2
    
    def compute_n_robustness(self, n_sim, method, **kwargs):
        # set defaults and extract eval params
        kwargs.setdefault('maxlen', kwargs.get('n', 1))
        phi_eval = kwargs.get('phi_eval', 'nonsym')
        by_score = kwargs.get('by_score', False)
        threshold = kwargs.get('threshold', None)

        # generate trajectories 
        trajs = self.generate_trajectories(
            n_sim,
            method=method,
            eval_method='max_length',
            **kwargs
        )

        # check whether to use multiprocessing
        use_mp = kwargs.pop('multiprocessing', False)
        m = 0

        if use_mp:
            # decide number of processes
            default_n_cpu = max(1, (os.cpu_count() or 1) - 2)
            n_cpu = kwargs.pop('n_cpu', default_n_cpu)

            # prepare arguments for each trajectory
            arg_list = [
                (traj, phi_eval, by_score, threshold)
                for traj in trajs
            ]
            with Pool(n_cpu) as pool:
                contributions = pool.starmap(_traj_contribution, arg_list)
            m = sum(contributions)

        else:
            # serial 
            for traj in trajs:
                m += _traj_contribution((traj, phi_eval, by_score, threshold))
        
        if by_score:           
            robustness = 1/(m/n_sim)
        else:
            robustness = m / n_sim
        #self.trajs = trajs

        self.trajs[f"nrho-{kwargs.get('n', 1)}"] = {
            "robustness": robustness,
            "trajs": trajs,
            "parameters": {
                "n_sim": n_sim,
                "method": method,
                **kwargs
            }
        }

        return robustness

    # def compute_n_robustness(self, n_sim, method, **kwargs):    
    #     kwargs.setdefault('maxlen', kwargs.get('n', 1))
    #     trajs = self.generate_trajectories(n_sim, 
    #                                        method=method, 
    #                                        eval_method='max_length', 
    #                                        **kwargs)

    #     phi_eval = kwargs.get('phi_eval', 'nonsym')
    #     by_score = kwargs.get('by_score', False)

    #     m = 0
    #     for traj in trajs:
    #         records = traj.get_records()
    #         orig_seq = records[0].sequence[1:-1]
    #         last_seq = records[-1].sequence[1:-1] 
            
    #         if phi_eval == 'nonsym':
    #             if self._compare_protein_sequences(orig_seq, last_seq):
    #                 m += 1
    #         elif phi_eval == 'blosum':
    #             if by_score:
    #                 alignment = compute_alignment(orig_seq, last_seq)
    #                 m += alignment.score
    #             else:
    #                 threshold = kwargs.get('threshold', None)
    #                 if threshold is None:
    #                     raise ValueError("threshold must be provided")

    #                 dist = compute_align_distance(orig_seq, last_seq)
                
    #                 if dist <= threshold:
    #                     m += 1

    #     robustness = m/n_sim

    #     self.trajs = trajs
                
    #     return robustness

    