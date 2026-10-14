from sinabro.evaluation.policies import protein

def eval_nonsym(records, mutinfo, **kwargs):
    if mutinfo.e:
        return (True, None)
    start = kwargs.get('cds_start', 1)
    end = kwargs.get('cds_end', len(records[0].sequence) - 1)
    return (protein(records[0].sequence[start:end]) != protein(records[-1].sequence[start:end]), None)
