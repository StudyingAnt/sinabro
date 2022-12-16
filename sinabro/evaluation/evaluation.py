import subprocess

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

def get_stop_condition(cur_seq, ori_score, **kwargs):
    if "measure" in kwargs.keys():
        measure = kwargs["measure"]
    else:
        raise ValueError(
                    "measure must be given"
                )

    if measure == "blastp":
        threshold = kwargs["threshold"]
    
        cur_seq = cur_seq[1:len(cur_seq)-1]
        cur_pp_seq = cur_seq.translate()
        record = SeqRecord(cur_pp_seq, id="curseq", description="")
        SeqIO.write(record, "tmp_cur_seq.fa", "fasta")
    
        rlt_score = subprocess.run(
            [
                "blastp",
                "-subject", "tmp_std_seq.fa",
                "-query", "tmp_cur_seq.fa",
                "-outfmt", "6 bitscore"
            ],
            stdout=subprocess.PIPE
        ).stdout.decode("utf-8")
        if not rlt_score:
            return True, "0"
        else:
            rlt_score = float(rlt_score.replace("\n", ""))

        diff = ori_score-rlt_score
        if diff > threshold:
            return True, str(rlt_score)
        else:
            return False, str(rlt_score)

