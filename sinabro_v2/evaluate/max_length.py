from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord

def eval_maxlen(records, mutinfo, **kwargs):
    maxlen = kwargs.get('maxlen', 10)

    if len(records)-1 == maxlen:
        return True, None
    else:
        return False, None