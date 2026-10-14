"""Optional pairwise sequence dN/dS; this is not tumor-selection inference."""
import numpy as np
from ..core import STOP_CODONS, coding_bounds, dna

def _has_inframe_stop(seq_without_terminal):
    text = dna(seq_without_terminal)
    if len(text) % 3:
        raise ValueError('CDS must have complete codons')
    return any((text[i:i + 3] in STOP_CODONS for i in range(0, len(text), 3)))

def compute_dnds(records, **kwargs):
    if not records:
        raise ValueError('at least one record is required')
    include = kwargs.get('threshold_hit', False)
    if not include and len(records) < 2:
        raise ValueError('excluding terminal record requires at least two records')
    start = kwargs.get('cds_start', 1)
    end = kwargs.get('cds_end', len(records[0].sequence) - 1)
    (ref, _, _) = coding_bounds(records[0].sequence, start, end)
    (cur, _, _) = coding_bounds(records[-1 if include else -2].sequence, start, end)
    (ref, cur) = (ref[start:end], cur[start:end])
    if len(ref) != len(cur) or len(ref) < 6 or ref[-3:] not in STOP_CODONS:
        raise ValueError('pairwise dN/dS requires equal-length CDS with a reference terminal stop and at least one sense codon')
    stoploss = cur[-3:] not in STOP_CODONS
    (ref, cur) = (ref[:-3], cur[:-3])
    premature = _has_inframe_stop(cur)
    if premature or _has_inframe_stop(ref):
        return (np.nan, np.nan, np.nan, premature, stoploss)
    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    (dn, ds) = cal_dn_ds(CodonSeq(ref), CodonSeq(cur), method=kwargs.get('method', 'NG86'))
    if not np.isfinite(dn) or not np.isfinite(ds) or dn < 0 or ds < 0:
        return np.nan, np.nan, np.nan, premature, stoploss
    omega = dn / ds if ds > 0 else np.inf if dn > 0 else np.nan
    return (dn, ds, omega, premature, stoploss)
