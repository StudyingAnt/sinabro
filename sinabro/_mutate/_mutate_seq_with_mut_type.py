import time
import numpy as np

from Bio.Seq import Seq, MutableSeq
from typing import Union

from .._helper._helper import (
    MutInfo,
    _unpack_complex_mut_type,
    _reverse_complement_mut_type,
    _get_target_of_mut_type,
    _get_result_of_mut_type,
    _get_motif_indices,
    _get_idx_offset_of_mut_type,
    _get_idx_target_from_idx_motif,
    _get_idx_motif_from_idx_target
)


def mutate_seq_with_mut_type(seq, mut_type, strand = "both", start=None, end=None):
    if isinstance(seq, str):
        seq = MutableSeq(seq)
        data_type = "str"
    elif isinstance(seq, Seq):
        seq = MutableSeq(seq)
        data_type = "Seq"
    elif isinstance(seq, MutableSeq):
        data_type = "MutableSeq"
    else:
        raise TypeError(
            "seq should be a string, Seq, or MutableSeq object"
        )

    if start is not None:
        if isinstance(start, int):
            pass                
        else:
            raise TypeError(
                "start should be an integer object"
            )
    else:
        start = 0
        
    if end is not None:
        if isinstance(end, int):
            pass                
        else:
            raise TypeError(
                "end should be an integer object"
            )
    else:
        end = len(seq)-1

    if strand != "both":
        if isinstance(strand, str):
            if strand == "single":
                pass
            else:
                raise ValueError(
                    "strand should be either single or both, default is both"
                )
        else:
            raise TypeError(
                "strand should be either single or both in string, default is both"
            )

    # check if mut_type has other symbols then A, C, G, or T    
    mut_type_nts = [nt for nt in mut_type if nt not in ["[", ">", "]"]]
    non_single_iupac_nts = [
        "R", "Y", "S", "W", "K", "M", 
        "B", "D", "H", "V",
        "N"
        ]
    
    if set(mut_type_nts) & set(non_single_iupac_nts):
        simple_mut_type = False
    else:
        simple_mut_type = True
        

    if simple_mut_type:
        # strand == single
        if strand == "single":           
            idx_motifs = _get_motif_indices(seq, mut_type)
            idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

            if not idx_motifs:
                return MutInfo(MutableSeq(""), -1, "", "", 1)
            
            # select specific motif to mutate
            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            idx_motif = np.random.choice(idx_motifs, 1)[0] 
            idx_target = _get_idx_target_from_idx_motif(idx_motif, mut_type)

            nt_before = seq[idx_target:idx_target+1]
            len_motif = len(_get_result_of_mut_type(mut_type))
            seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(mut_type)
            nt_after = seq[idx_target:idx_target+1]
            
            hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

            if data_type == "str":
                return MutInfo(str(seq), idx_target, hgvs_mrna, mut_type, 0)
            elif data_type == "Seq":
                return MutInfo(Seq(seq), idx_target, hgvs_mrna, mut_type, 0)
            elif data_type == "MutableSeq":
                return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)        
        else: # strand == both
            revcomp_mut_type = _reverse_complement_mut_type(mut_type)

            idx_motifs_pos = _get_motif_indices(seq, mut_type)
            idx_motifs_neg = _get_motif_indices(seq, revcomp_mut_type)

            idx_motifs = idx_motifs_pos + idx_motifs_neg
            idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

            if not idx_motifs:
                return MutInfo(MutableSeq(""), -1, "", "", 1)
            
            
            # select specific motif to mutate
            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            idx_motif = np.random.choice(idx_motifs, 1)[0] 

            # run different mut_type based on the strand
            if idx_motif in idx_motifs_pos:
                idx_target = _get_idx_target_from_idx_motif(idx_motif, mut_type)

                nt_before = seq[idx_target:idx_target+1]
                len_motif = len(_get_result_of_mut_type(mut_type))
                seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(mut_type)
                nt_after = seq[idx_target:idx_target+1]
                
                hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

                if data_type == "str":
                    return MutInfo(str(seq), idx_target, hgvs_mrna, mut_type, 0)
                elif data_type == "Seq":
                    return MutInfo(Seq(seq), idx_target, hgvs_mrna, mut_type, 0)
                elif data_type == "MutableSeq":
                    return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)
            else:
                idx_target = _get_idx_target_from_idx_motif(idx_motif, revcomp_mut_type)

                nt_before = seq[idx_target:idx_target+1]
                len_motif = len(_get_result_of_mut_type(revcomp_mut_type))
                seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(revcomp_mut_type)
                nt_after = seq[idx_target:idx_target+1]
                
                hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

                if data_type == "str":
                    return MutInfo(str(seq), idx_target, hgvs_mrna, revcomp_mut_type, 0)
                elif data_type == "Seq":
                    return MutInfo(Seq(seq), idx_target, hgvs_mrna, revcomp_mut_type, 0)
                elif data_type == "MutableSeq":
                    return MutInfo(seq, idx_target, hgvs_mrna, revcomp_mut_type, 0)
    else: # if mut type has non single nucleotide iupac symbols
        mut_types = _unpack_complex_mut_type(mut_type)
        # strand == single
        if strand == "single":
            # select one specific mut type
            cnt_mut_types = []
            sum = 0
            for unpacked_mut_type in mut_types:
                new_motifs = _get_motif_indices(seq, unpacked_mut_type)
                cnt_mut_types.append(len(new_motifs))
                sum = sum + len(new_motifs)

            if sum == 0:
                return MutInfo(MutableSeq(""), -1, "", "", 1)

            prob_mut_types = [cnt/sum for cnt in cnt_mut_types]

            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            mut_type_idx = np.random.choice(range(len(mut_types)), 1, p=prob_mut_types)[0]

            selected_mut_type = mut_types[mut_type_idx]

            # find idx of selected mut type
            idx_motifs = _get_motif_indices(seq, selected_mut_type)
            idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

            if not idx_motifs:
                return MutInfo(MutableSeq(""), -1, "", "", 1)
            
            # select specific motif to mutate
            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            idx_motif = np.random.choice(idx_motifs, 1)[0] 
            idx_target = _get_idx_target_from_idx_motif(idx_motif, selected_mut_type)

            nt_before = seq[idx_target:idx_target+1]
            len_motif = len(_get_result_of_mut_type(selected_mut_type))
            seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(selected_mut_type)
            nt_after = seq[idx_target:idx_target+1]
            
            hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

            if data_type == "str":
                return MutInfo(str(seq), idx_target, hgvs_mrna, selected_mut_type, 0)
            elif data_type == "Seq":
                return MutInfo(Seq(seq), idx_target, hgvs_mrna, selected_mut_type, 0)
            elif data_type == "MutableSeq":
                return MutInfo(seq, idx_target, hgvs_mrna, selected_mut_type, 0)
                    
        else: # strand == both
            # select one specific mut type
            cnt_mut_types = []
            sum = 0
            for unpacked_mut_type in mut_types:
                new_motifs = _get_motif_indices(seq, unpacked_mut_type)
                cnt_mut_types.append(len(new_motifs))
                sum = sum + len(new_motifs)
            
            if sum == 0:
                return MutInfo(MutableSeq(""), -1, "", "", 1)

            prob_mut_types = [cnt/sum for cnt in cnt_mut_types]

            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            mut_type_idx = np.random.choice(range(len(mut_types)), 1, p=prob_mut_types)[0]

            selected_mut_type = mut_types[mut_type_idx]

            # find idx of selected mut type
            revcomp_mut_type = _reverse_complement_mut_type(selected_mut_type)

            idx_motifs_pos = _get_motif_indices(seq, selected_mut_type)
            idx_motifs_neg = _get_motif_indices(seq, revcomp_mut_type)

            idx_motifs = idx_motifs_pos + idx_motifs_neg
            idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

            if not idx_motifs:
                return MutInfo(MutableSeq(""), -1, "", "", 1)
            
            # select specific motif to mutate
            seed1 = time.time()
            seed2 = (seed1 - int(seed1))*(10**7)
            np.random.seed(int(seed1+seed2))
            idx_motif = np.random.choice(idx_motifs, 1)[0] 

            # run different mut_type based on the strand
            if idx_motif in idx_motifs_pos:
                idx_target = _get_idx_target_from_idx_motif(idx_motif, selected_mut_type)

                nt_before = seq[idx_target:idx_target+1]
                len_motif = len(_get_result_of_mut_type(selected_mut_type))
                seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(selected_mut_type)
                nt_after = seq[idx_target:idx_target+1]
                
                hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

                if data_type == "str":
                    return MutInfo(str(seq), idx_target, hgvs_mrna, selected_mut_type, 0)
                elif data_type == "Seq":
                    return MutInfo(Seq(seq), idx_target, hgvs_mrna, selected_mut_type, 0)
                elif data_type == "MutableSeq":
                    return MutInfo(seq, idx_target, hgvs_mrna, selected_mut_type, 0)
            else:
                idx_target = _get_idx_target_from_idx_motif(idx_motif, revcomp_mut_type)

                nt_before = seq[idx_target:idx_target+1]
                len_motif = len(_get_result_of_mut_type(revcomp_mut_type))
                seq[idx_motif:idx_motif+len_motif] = _get_result_of_mut_type(revcomp_mut_type)
                nt_after = seq[idx_target:idx_target+1]
                
                hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

                if data_type == "str":
                    return MutInfo(str(seq), idx_target, hgvs_mrna, revcomp_mut_type, 0)
                elif data_type == "Seq":
                    return MutInfo(Seq(seq), idx_target, hgvs_mrna, revcomp_mut_type, 0)
                elif data_type == "MutableSeq":
                    return MutInfo(seq, idx_target, hgvs_mrna, revcomp_mut_type, 0)




        