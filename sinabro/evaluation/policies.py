"""Neutrality criteria evaluated against the original coding sequence."""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import math
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
from ..core import EvaluationResult, dna

@lru_cache(maxsize=16)
def _matrix(name):
    return substitution_matrices.load(name)

def protein(seq):
    text = dna(seq)
    if len(text) % 3:
        raise ValueError('translation requires complete codons')
    return str(Seq(text).translate())

def aligner(matrix_name='BLOSUM80', mode='global', gap_open=-10, gap_extend=-0.5):
    if mode not in ('global', 'local'):
        raise ValueError('alignment mode must be global or local')
    if not all((math.isfinite(float(x)) and x <= 0 for x in (gap_open, gap_extend))):
        raise ValueError('gap penalties must be finite and nonpositive')
    result = PairwiseAligner()
    result.mode = mode
    result.open_gap_score = gap_open
    result.extend_gap_score = gap_extend
    result.substitution_matrix = _matrix(matrix_name)
    return result

@lru_cache(maxsize=128)
def _self_score(original_protein, matrix_name, mode, gap_open, gap_extend):
    return float(aligner(matrix_name, mode, gap_open, gap_extend).score(original_protein, original_protein))

def alignment_distance(original_cds, current_cds, *, matrix_name='BLOSUM80', mode='global', gap_open=-10, gap_extend=-0.5):
    (original, current) = (protein(original_cds), protein(current_cds))
    a = aligner(matrix_name, mode, gap_open, gap_extend)
    baseline = _self_score(original, matrix_name, mode, gap_open, gap_extend)
    return float(baseline - a.score(original, current))

@dataclass(frozen=True)
class ProteinIdentity:
    name: str = 'protein_identity'

    def __call__(self, original_cds, current_cds, mutinfo):
        return EvaluationResult(protein(original_cds) != protein(current_cds), reason='nonsynonymous')

@dataclass(frozen=True)
class BlosumDistance:
    threshold: float
    matrix_name: str = 'BLOSUM80'
    gap_open: float = -10
    gap_extend: float = -0.5
    mode: str = 'global'
    stop_policy: str = 'score'
    name: str = 'blosum_distance'

    def __post_init__(self):
        if isinstance(self.threshold, bool) or not math.isfinite(float(self.threshold)) or self.threshold <= 0:
            raise ValueError('threshold must be positive and finite')
        if self.stop_policy not in ('score', 'terminate'):
            raise ValueError("stop_policy must be 'score' or 'terminate'")
        aligner(self.matrix_name, self.mode, self.gap_open, self.gap_extend)

    def __call__(self, original_cds, current_cds, mutinfo):
        (original, current) = (protein(original_cds), protein(current_cds))
        a = aligner(self.matrix_name, self.mode, self.gap_open, self.gap_extend)
        baseline = _self_score(original, self.matrix_name, self.mode, self.gap_open, self.gap_extend)
        score = float(baseline - a.score(original, current))
        if self.stop_policy == 'terminate' and any(((x == '*') != (y == '*') for (x, y) in zip(original, current))):
            return EvaluationResult(True, score, 'stop_change')
        return EvaluationResult(score >= self.threshold, score, 'evaluation')

def evaluate_callable(evaluator, original, current, mutinfo):
    result = evaluator(original, current, mutinfo)
    if isinstance(result, EvaluationResult):
        if not isinstance(result.stop, bool) or not isinstance(result.reason, str) or (not result.reason):
            raise ValueError('invalid custom EvaluationResult')
        if result.score is not None and (not math.isfinite(float(result.score))):
            raise ValueError('evaluation score must be finite')
        if result.stop and result.reason in ('max_length', 'max_steps', 'no_available_event'):
            raise ValueError('custom terminal reason is reserved for incomplete/fixed-step trajectories')
        return result
    if isinstance(result, bool):
        return EvaluationResult(result)
    raise TypeError('custom evaluator must return bool or EvaluationResult')
