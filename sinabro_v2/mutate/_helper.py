from functools import wraps
from Bio.Seq import Seq, MutableSeq
from ..types.types import MutInfo, MutationRecord

def preserve_seq_type(func):
    """
    Decorator that preserves the original sequence type.
    
    This decorator converts the input sequence to a MutableSeq for processing
    and converts it back to its original type before returning the result.
    """

    @wraps(func)
    def wrapper(seq, *args, **kwargs):
        if isinstance(seq, str):
            convert = str
            mutable_seq = MutableSeq(seq)
        elif isinstance(seq, Seq):
            convert = Seq
            mutable_seq = MutableSeq(seq)
        elif isinstance(seq, MutableSeq):
            convert = lambda x: x
            mutable_seq = MutableSeq(str(seq))
        else:
            raise TypeError('seq should be a string, Seq, or MutableSeq object')
        result = func(mutable_seq, *args, **kwargs)
        converted_seq = convert(result.new_seq)
        return MutInfo(converted_seq, result.idx_target, result.hgvs_mrna, result.mut_type, result.e, result.reason)
    return wrapper
