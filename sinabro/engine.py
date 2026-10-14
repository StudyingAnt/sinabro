"""Trajectory API: mutate once, evaluate, then repeat until termination."""
from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import math
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .core import MutationRecord, MutInfo, EvaluationResult, coding_bounds, integer, dna, rng_from, get_codon, amino_acid, random_mutation, motif_mutation, signature_mutation
from .evaluation.policies import ProteinIdentity, BlosumDistance, evaluate_callable, protein

class Trajectory:

    def __init__(self, id, data, note='.', *, seed=None, rng=None, cds_start=1, cds_end=None):
        (text, self.cds_start, self.cds_end) = coding_bounds(data, cds_start, cds_end)
        self._id = self._traj_id = id
        self._original_seq = Seq(text)
        self._seq_length = len(text)
        self._records = [MutationRecord(self._original_seq, note=note)]
        self.rng = rng_from(rng, seed)
        self.seed = seed
        self.terminal_reason = None
        self.evaluation = None
        self.evaluation_parameters = {}
        self.score = None
        self._sync()

    def _sync(self):
        self._data = [r.sequence for r in self._records]
        self._hgvs_mrnas = [r.hgvs_mrna for r in self._records]
        self._hgvs_aa = [r.hgvs_aa for r in self._records]
        self._mut_types = [r.mut_type for r in self._records]
        self._notes = [r.note for r in self._records]
        self._length = len(self._records) - 1

    def get_records(self):
        """Return a list of mutation records, including the original sequence."""
        return [replace(record) for record in self._records]

    def get_length(self):
        return self._length

    def __str__(self):
        rows = [f'ID: {self._id}', f'Trajectory length: {self._length}', f'Terminal reason: {self.terminal_reason}', 'Sequence\tHGVS_mRNA\tHGVS_Protein\tMutation_Type\tNote']
        rows.extend(('\t'.join(map(str, (r.sequence, r.hgvs_mrna, r.hgvs_aa, r.mut_type, r.note))) for r in self._records))
        return '\n'.join(rows)

    def show(self, verbose=True, noseq=False, show_idx=False):
        if verbose:
            print(str(self) if not noseq else '\n'.join(('\t'.join(map(str, (r.hgvs_mrna, r.hgvs_aa, r.mut_type, r.note))) for r in self._records)))
        else:
            print(f'ID: {self._id}\nTrajectory length: {self._length}')

    def append(self, seq, hgvs_mrna='.', mut_type='.', note='.'):
        """Append exactly one in-CDS substitution and synchronize compatibility arrays."""
        text = dna(seq)
        previous = str(self._records[-1].sequence)
        if len(text) != self._seq_length:
            raise ValueError('trajectory sequences must have the same length')
        changed = [i for (i, (a, b)) in enumerate(zip(previous, text)) if a != b]
        if len(changed) != 1:
            raise ValueError('each record must represent exactly one substitution')
        idx = changed[0]
        if not self.cds_start <= idx < self.cds_end:
            raise ValueError('mutation target is outside the CDS')
        old_codon = get_codon(previous, idx, self.cds_start, self.cds_end)
        new_codon = get_codon(text, idx, self.cds_start, self.cds_end)
        (old, new) = (amino_acid(old_codon), amino_acid(new_codon))
        if hgvs_mrna == '.':
            hgvs_mrna = f'c.{idx - self.cds_start + 1}{previous[idx]}>{text[idx]}'
        pos = (idx - self.cds_start) // 3 + 1
        annotation = '.' if old == new else f'p.{pos}{amino_acid(old_codon, True)}>{amino_acid(new_codon, True)}'
        self._records.append(MutationRecord(Seq(text), hgvs_mrna, annotation, mut_type, note))
        self.terminal_reason = self.score = None
        self._sync()
        return self._data

    def add_record(self, record):
        if not isinstance(record, MutationRecord):
            raise TypeError('record must be MutationRecord')
        self.append(record.sequence, record.hgvs_mrna, record.mut_type, record.note)
        return self.get_records()

    def remove_last_record(self):
        if self._length == 0:
            raise IndexError('cannot remove the original sequence')
        result = self._records.pop()
        self.terminal_reason = self.score = None
        self._sync()
        return result

    def pop(self):
        self.remove_last_record()
        return self._data

    def _configure_bounds(self, kwargs):
        if 'start' in kwargs and 'cds_start' in kwargs:
            raise ValueError('provide start or cds_start, not both')
        if 'end' in kwargs and 'cds_end' in kwargs:
            raise ValueError('provide end or cds_end, not both')
        start = kwargs.pop('start', kwargs.pop('cds_start', self.cds_start))
        end = kwargs.pop('end', kwargs.pop('cds_end', self.cds_end))
        (_, start, end) = coding_bounds(self._original_seq, start, end)
        if self._length and (start, end) != (self.cds_start, self.cds_end):
            raise ValueError('cannot change CDS coordinates after mutations')
        (self.cds_start, self.cds_end) = (start, end)

    def add_mutated_sequence(self, method='random', **kwargs):
        """Sample one event with this trajectory RNG; return MutInfo."""
        self._configure_bounds(kwargs)
        note = kwargs.pop('note', '.')
        if 'seed' in kwargs or 'rng' in kwargs:
            raise ValueError('pass seed/rng to Trajectory, not per mutation')
        common = dict(start=self.cds_start, end=self.cds_end, cds_start=self.cds_start, rng=self.rng)
        seq = self._records[-1].sequence
        if method == 'random':
            if kwargs:
                raise TypeError(f'unknown random mutation arguments: {sorted(kwargs)}')
            info = random_mutation(seq, **common)
        elif method in ('mut_type', 'mut_types'):
            required = ('mut_type',) if method == 'mut_type' else ('mut_types', 'mut_type_probs')
            if any(name not in kwargs for name in required):
                raise ValueError('missing required motif mutation arguments: ' + ', '.join(required))
            motifs = [kwargs.pop('mut_type')] if method == 'mut_type' else kwargs.pop('mut_types')
            weights = None if method == 'mut_type' else kwargs.pop('mut_type_probs')
            info = motif_mutation(seq, motifs, weights, **common, **kwargs)
        elif method == 'signature':
            info = signature_mutation(seq, **common, **kwargs)
        else:
            raise ValueError('method must be random, mut_type, mut_types, or signature')
        if info.e:
            self.terminal_reason = info.reason or 'no_available_event'
        else:
            self.append(info.new_seq, info.hgvs_mrna, info.mut_type, note)
        return info

    def autofill(self, condition='max_length', max_length=10, *, max_steps=10000, evaluator=None, **kwargs):
        """Generate substitutions until a criterion, fixed length, no event, or safety cap."""
        if self.terminal_reason not in (None, 'max_length', 'max_steps', 'no_available_event'):
            raise ValueError('trajectory already terminated; create a new trajectory or remove the terminal record')
        self._configure_bounds(kwargs)
        method = kwargs.pop('method', 'random')
        if 'seed' in kwargs:
            if self._length or self.seed is not None:
                raise ValueError('trajectory seed already initialized')
            self.seed = kwargs.pop('seed')
            self.rng = rng_from(seed=self.seed)
        if 'rng' in kwargs:
            raise ValueError('pass rng to Trajectory constructor')
        max_steps = integer(max_steps, 'max_steps', 1)
        if condition == 'max_length':
            max_length = integer(max_length, 'max_length')
            policy = None
            self.evaluation = 'fixed_steps'
            self.evaluation_parameters = {'max_length': max_length}
        elif condition == 'nonsynonymous':
            if evaluator is not None:
                raise ValueError('use condition=evaluation for a custom evaluator')
            policy = ProteinIdentity()
            self.evaluation = policy.name
            self.evaluation_parameters = {}
        elif condition in ('evaluation', 'blosum'):
            if evaluator is None:
                measure = kwargs.pop('measure', 'blosum')
                if measure not in ('blosum', 'BLOSUM80'):
                    raise ValueError('supported built-in measure is blosum; BLASTP is not silently substituted')
                settings = {k: kwargs.pop(k) for k in ('threshold', 'matrix_name', 'gap_open', 'gap_extend', 'mode', 'stop_policy') if k in kwargs}
                if 'threshold' not in settings:
                    raise ValueError('BLOSUM evaluation requires threshold')
                policy = BlosumDistance(**settings)
                self.evaluation = policy.name
                self.evaluation_parameters = asdict(policy)
            else:
                if not callable(evaluator):
                    raise TypeError('evaluator must be callable')
                policy = evaluator
                self.evaluation = getattr(evaluator, 'name', 'custom')
                self.evaluation_parameters = {}
        else:
            raise ValueError('unknown stop condition')
        if self._length and isinstance(policy, (ProteinIdentity, BlosumDistance)):
            original = str(self._original_seq[self.cds_start:self.cds_end])
            current = str(self._records[-1].sequence[self.cds_start:self.cds_end])
            if policy(original, current, None).stop:
                raise ValueError('current sequence already violates the selected neutrality criterion')
        self.terminal_reason = self.score = None
        for _ in range(max_steps):
            if condition == 'max_length' and self._length >= max_length:
                self.terminal_reason = 'max_length'
                break
            info = self.add_mutated_sequence(method, **kwargs)
            if info.e:
                break
            if policy is not None:
                result = evaluate_callable(policy, str(self._original_seq[self.cds_start:self.cds_end]), str(self._records[-1].sequence[self.cds_start:self.cds_end]), info)
                self.score = result.score
                if result.score is not None:
                    self._records[-1].note = f'score: {result.score}'
                    self._sync()
                if result.stop:
                    self.terminal_reason = result.reason
                    break
            elif self._length >= max_length:
                self.terminal_reason = 'max_length'
                break
        else:
            self.terminal_reason = 'max_steps'
        return self._data

    @property
    def accepted_steps(self):
        if self.terminal_reason not in ('max_length', 'max_steps', 'no_available_event', None):
            return self._length - 1
        return None

    @property
    def synonymous_steps(self):
        return sum((r.hgvs_aa == '.' for r in self._records[1:]))

    def save(self, output_dir_path, prefix=None, suffix=None):
        directory = Path(output_dir_path)
        directory.mkdir(parents=True, exist_ok=True)
        name = (f'{prefix}_' if prefix else '') + 'trajectory' + (f'_{suffix}' if suffix else '') + '.fa'
        records = [SeqRecord(r.sequence, id=str(self._id), description=';'.join(map(str, (r.hgvs_mrna, r.hgvs_aa, r.mut_type, r.note)))) for r in self._records]
        path = directory / name
        SeqIO.write(records, path, 'fasta')
        return path

def _trajectory_job(job):
    (seq, idx, seed, method, condition, kwargs, coords) = job
    traj = Trajectory(idx, seq, seed=seed, cds_start=coords[0], cds_end=coords[1])
    traj.autofill(condition, method=method, **kwargs)
    return traj

class RobustnessComputer:

    def __init__(self, gene_name, gene_seq, *, seed=None, cds_start=1, cds_end=None):
        (text, self.cds_start, self.cds_end) = coding_bounds(gene_seq, cds_start, cds_end)
        (self.gene_name, self.gene_seq) = (gene_name, Seq(text))
        self.seed = seed
        self.trajs = {}

    @staticmethod
    def _condition(eval_method):
        if eval_method == 'blosum':
            return 'evaluation'
        if eval_method in ('nonsynonymous', 'max_length', 'evaluation'):
            return eval_method
        raise ValueError('unknown evaluation method')

    def generate_trajectory(self, traj_id, method, eval_method, **kwargs):
        seeds = np.random.SeedSequence(self.seed).spawn(integer(traj_id, 'traj_id') + 1)
        return _trajectory_job((self.gene_seq, traj_id, seeds[-1], method, self._condition(eval_method), kwargs, (self.cds_start, self.cds_end)))

    def generate_trajectories(self, n_traj, method, eval_method, **kwargs):
        """Generate indexed independent streams, invariant to worker scheduling."""
        n_traj = integer(n_traj, 'n_traj', 1)
        use_mp = kwargs.pop('multiprocessing', False)
        workers = integer(kwargs.pop('n_cpu', 2), 'n_cpu', 1)
        if 'maxlen' in kwargs:
            kwargs['max_length'] = kwargs.pop('maxlen')
        seed = kwargs.pop('seed', self.seed)
        seeds = np.random.SeedSequence(seed).spawn(n_traj)
        condition = self._condition(eval_method)
        jobs = [(self.gene_seq, i, seeds[i], method, condition, dict(kwargs), (self.cds_start, self.cds_end)) for i in range(n_traj)]
        if use_mp:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                return list(pool.map(_trajectory_job, jobs))
        return [_trajectory_job(job) for job in jobs]

    def compute_l_robustness(self, n_sim, method, eval_method, **kwargs):
        """Average complete pre-terminal neutral lengths; reject incomplete paths."""
        trajs = self.generate_trajectories(n_sim, method, eval_method, **kwargs)
        incomplete = [t for t in trajs if t.accepted_steps is None]
        result = {'trajs': trajs, 'parameters': dict(n_sim=n_sim, method=method, eval_method=eval_method, **kwargs), 'n_complete': len(trajs) - len(incomplete), 'n_incomplete': len(incomplete), 'evaluation': trajs[0].evaluation, 'robustness': None}
        self.trajs['lrho'] = result
        if incomplete:
            raise ValueError('stopping-time mean undefined for incomplete trajectories; inspect trajs[lrho]')
        values = np.array([t.accepted_steps for t in trajs], dtype=float)
        result['robustness'] = float(np.mean(values))
        result['median'] = float(np.median(values))
        result['sd'] = float(np.std(values, ddof=1)) if len(values) > 1 else None
        result['mc_se'] = result['sd'] / math.sqrt(len(values)) if result['sd'] is not None else None
        result['p_zero'] = float(np.mean(values == 0))
        result['metric'] = 'RUMS' if result['evaluation'] == 'protein_identity' else 'criterion_neutral_stopping_time'
        return result['robustness']

    def _compare_protein_sequences(self, seq1, seq2):
        return protein(seq1) == protein(seq2)

    def compute_n_robustness(self, n_sim, method, **kwargs):
        """Measure endpoint neutrality after exactly n substitutions."""
        n = integer(kwargs.pop('n', kwargs.pop('maxlen', 1)), 'n')
        phi = kwargs.pop('phi_eval', 'nonsym')
        by_score = kwargs.pop('by_score', False)
        threshold = kwargs.pop('threshold', None)
        settings = {k: kwargs.pop(k) for k in ('matrix_name', 'mode', 'gap_open', 'gap_extend') if k in kwargs}
        if phi not in ('nonsym', 'blosum'):
            raise ValueError('unknown phenotype evaluation')
        if phi == 'nonsym' and by_score:
            raise ValueError('by_score requires blosum')
        if phi == 'blosum' and (not by_score) and (threshold is None or not math.isfinite(float(threshold)) or threshold <= 0):
            raise ValueError('blosum phenotype evaluation requires positive finite threshold')
        trajs = self.generate_trajectories(n_sim, method, 'max_length', max_length=n, **kwargs)
        self.trajs[f'nrho-{n}'] = {'trajs': trajs, 'parameters': dict(n=n, phi_eval=phi, by_score=by_score, threshold=threshold, **settings, **kwargs), 'robustness': None}
        if any((t.terminal_reason != 'max_length' or t._length != n for t in trajs)):
            raise ValueError('fixed-n trajectories must complete exactly n substitutions')
        contributions = []
        from .evaluation.policies import alignment_distance
        for t in trajs:
            original = t._original_seq[t.cds_start:t.cds_end]
            current = t._records[-1].sequence[t.cds_start:t.cds_end]
            if phi == 'nonsym':
                contributions.append(float(protein(original) == protein(current)))
            else:
                distance = alignment_distance(original, current, **settings)
                contributions.append(distance if by_score else float(distance < threshold))
        mean = float(np.mean(contributions))
        result = (math.inf if mean == 0 else 1 / mean) if by_score else mean
        self.trajs[f'nrho-{n}']['robustness'] = result
        return result
