from typing import Union
from Bio.Seq import Seq, MutableSeq

from .types.types import MutInfo, MutationRecord

def eval_maxlen(records, **kwargs):
    