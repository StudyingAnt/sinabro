"""Compatibility evaluator; external BLASTP is an explicit unsupported legacy mode."""
from .policies import BlosumDistance

def get_stop_condition(cur_seq, ori_score=None, **kwargs):
    if kwargs.get('measure') == 'blastp':
        raise ValueError('legacy BLASTP mode is not supported; use explicit BLOSUM evaluation')
    original = kwargs.get('original_seq')
    if original is None:
        raise ValueError('original_seq required for BLOSUM evaluation')
    result = BlosumDistance(kwargs['threshold'])(str(original[1:-1]), str(cur_seq[1:-1]), None)
    return (result.stop, result.score)
