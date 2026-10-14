from sinabro.evaluation.policies import protein as dna_to_aa, aligner, alignment_distance, BlosumDistance

def compute_alignment(seq1, seq2, mode='global', gap_open=-10, gap_extend=-0.5, matrix_name='BLOSUM80'):
    return aligner(matrix_name, mode, gap_open, gap_extend).align(dna_to_aa(seq1), dna_to_aa(seq2))[0]

def compute_align_distance(seq1, seq2, **kwargs):
    return alignment_distance(seq1, seq2, **kwargs)

def eval_blosum(records, mutinfo, **kwargs):
    if mutinfo.e:
        return (True, None)
    start = kwargs.get('cds_start', 1)
    end = kwargs.get('cds_end', len(records[0].sequence) - 1)
    policy = BlosumDistance(**{k: kwargs[k] for k in ('threshold', 'matrix_name', 'mode', 'gap_open', 'gap_extend', 'stop_policy') if k in kwargs})
    result = policy(str(records[0].sequence[start:end]), str(records[-1].sequence[start:end]), mutinfo)
    return (result.stop, result.score)
