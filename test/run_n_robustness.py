import warnings
from Bio import BiopythonExperimentalWarning

warnings.filterwarnings('ignore', category=BiopythonExperimentalWarning)

from pathlib import Path
import os
import sys

file_path = Path(__file__).resolve()
scr_path = file_path.parent.parent
sys.path.insert(0, str(scr_path))

import sinabro_v2 as snbr
from Bio import SeqIO
import signal
import h5py
from datetime import datetime

data_path = file_path.parent.parent.parent / "data"
lgenes_fasta_file = data_path / "gencode.v40.pc_transcripts.nopary.cdsplus.longest.fa"


seq_records = list(SeqIO.parse(lgenes_fasta_file, "fasta"))
partial_seq_records = seq_records[:10]

# lrho
with h5py.File('lrho_sbs2.h5', 
               mode='a',
               libver='latest') as f:
    f.swmr_mode = True
    if 'lrho' not in f:
        f.create_dataset(
            'lrho',             # dataset name
            shape=(0,),                  # initial length = 0
            maxshape=(None,),            # unlimited resize in first dim
            dtype='float64',             # float 결과
            chunks=True,                 # chunked storage (required for resize+compression)
            compression='gzip',          # gzip 압축
            compression_opts=9           # 압축 레벨 (1~9)
        )
    ds = f['lrho']
    try:
        for seq_record in partial_seq_records:
            transcript_name = seq_record.id.split("|")[4]
            seq = seq_record.seq

            # now
            now = datetime.now()
            formatted_time = now.strftime("%H:%M:%S.%f")[:-3]
            print(f"{formatted_time} RUNNING   {transcript_name}")
            rhocom = snbr.RobustnessComputer(transcript_name, seq)
        
            lrho = rhocom.compute_l_robustness(n_sim=1000, 
                                   method="signature", 
                                   eval_method="blosum", 
                                   mutational_signature='SBS2', 
                                   threshold = 20,
                                   multiprocessing=True)
            ds.resize((ds.shape[0] + 1,))
            # Append the new float value
            ds[-1] = lrho
    
            # Optional: 즉시 디스크에 flush
            f.flush()
            os.fsync(f.id.get_vfd_handle())

            now = datetime.now()
            formatted_time = now.strftime("%H:%M:%S.%f")[:-3]
            print(f"{formatted_time} COMPLETE  {transcript_name}")
    except Exception:
        # 예외 시에도 지금까지 쓴 데이터는 저장
        f.flush()
        os.fsync(f.id.get_vfd_handle())
        raise