"""Small portable package CLI. Manuscript-specific analysis lives separately."""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
from pathlib import Path
import platform
import sys
import numpy as np
import Bio
from Bio import SeqIO
from . import __version__
from .engine import RobustnessComputer
from .core import signature_path, integer

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='FASTA of CDS-plus sequences')
    parser.add_argument('--output', type=Path, required=True, help='New output CSV; refuses to overwrite')
    parser.add_argument('--signatures', nargs='+', default=['SBS2'])
    parser.add_argument('--n-trajectories', type=int, default=100)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--evaluation', choices=['protein_identity', 'blosum'], default='protein_identity')
    parser.add_argument('--threshold', type=float)
    parser.add_argument('--matrix', default='BLOSUM80')
    parser.add_argument('--stop-policy', choices=['score', 'terminate'], default='score')
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--cds-start', type=int, default=1)
    parser.add_argument('--cosmic-version', default='3.3.1')
    parser.add_argument('--genome-ref', default='GRCh37')
    args = parser.parse_args(argv)
    try:
        integer(args.n_trajectories, 'n-trajectories', 1)
        integer(args.max_steps, 'max-steps', 1)
        integer(args.workers, 'workers', 1)
        integer(args.seed, 'seed')
        integer(args.cds_start, 'cds-start')
        if args.evaluation == 'blosum' and args.threshold is None:
            raise ValueError('--threshold is required for BLOSUM evaluation')
        if args.evaluation == 'protein_identity' and args.threshold is not None:
            raise ValueError('--threshold applies only to BLOSUM evaluation')
        if args.output.exists() or args.output.with_suffix('.manifest.json').exists():
            raise ValueError('output already exists; choose a new output path')
        records = list(SeqIO.parse(args.input, 'fasta'))
        if not records:
            raise ValueError('input contains no FASTA records')
        if len(set((r.id for r in records))) != len(records):
            raise ValueError('FASTA IDs must be unique')
        if len(set(args.signatures)) != len(args.signatures):
            raise ValueError('signatures must be unique')
        rows = []
        for record in records:
            for signature in args.signatures:
                digest = hashlib.sha256(f'{args.seed}\x00{record.id}\x00{signature}'.encode()).digest()
                run_seed = int.from_bytes(digest[:8], 'big')
                end = len(record.seq) if args.cds_start == 0 else len(record.seq) - 1
                computer = RobustnessComputer(record.id, record.seq, seed=run_seed, cds_start=args.cds_start, cds_end=end)
                config = dict(mutational_signature=signature, cosmic_version=args.cosmic_version, genome_ref=args.genome_ref, max_steps=args.max_steps, multiprocessing=args.workers > 1, n_cpu=args.workers)
                if args.evaluation == 'blosum':
                    config.update(threshold=args.threshold, matrix_name=args.matrix, stop_policy=args.stop_policy)
                trajectories = computer.generate_trajectories(args.n_trajectories, 'signature', 'blosum' if args.evaluation == 'blosum' else 'nonsynonymous', **config)
                for t in trajectories:
                    rows.append(dict(sequence_id=record.id, signature=signature, trajectory_id=t._traj_id, evaluation=t.evaluation, terminal_reason=t.terminal_reason, n_steps=t._length, accepted_steps=t.accepted_steps, synonymous_steps=t.synonymous_steps, score=t.score, stream_seed=run_seed))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        data_path = signature_path(args.cosmic_version, args.genome_ref)
        manifest = dict(sinabro_version=__version__, python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__, input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(), signature_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(), output_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(), configuration={k: str(v) if isinstance(v, Path) else v for (k, v) in vars(args).items()}, rows=len(rows), complete=sum((r['accepted_steps'] is not None for r in rows)), incomplete=sum((r['accepted_steps'] is None for r in rows)))
        args.output.with_suffix('.manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        print(f'Wrote {len(rows)} trajectories to {args.output}')
        return 0
    except (ValueError, TypeError, OSError) as exc:
        parser.error(str(exc))
if __name__ == '__main__':
    sys.exit(main())
