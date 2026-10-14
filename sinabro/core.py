"""Shared mutation primitives. Indices here are zero based, with exclusive ends."""
from __future__ import annotations
import csv
import itertools
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional
import numpy as np
from Bio.Seq import Seq, MutableSeq
BASES = 'ACGT'
TRIPLETS = tuple((''.join(t) for t in itertools.product(BASES, repeat=3)))
TRIPLET_INDEX = {t: i for (i, t) in enumerate(TRIPLETS)}
IUPAC = dict(zip('RYSWKMBDHVN', ('AG', 'CT', 'GC', 'AT', 'GT', 'AC', 'CGT', 'AGT', 'ACT', 'ACG', 'ACGT')))
PATTERN = re.compile('^([ACGTRYSWKMBDHVN]*)\\[([ACGT])>([ACGT])\\]([ACGTRYSWKMBDHVN]*)$')
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

@dataclass
class MutInfo:
    new_seq: object
    idx_target: int
    hgvs_mrna: str
    mut_type: str
    e: int = 0
    reason: Optional[str] = None

@dataclass
class MutationRecord:
    sequence: object
    hgvs_mrna: str = '.'
    hgvs_aa: str = '.'
    mut_type: str = '.'
    note: str = '.'

@dataclass(frozen=True)
class EvaluationResult:
    stop: bool
    score: Optional[float] = None
    reason: str = 'evaluation'

def dna(seq):
    """Validate ACGT DNA and normalize case without mutating the input."""
    if not isinstance(seq, (str, Seq, MutableSeq)):
        raise TypeError('sequence must be str, Seq, or MutableSeq')
    text = str(seq).upper()
    if not text or set(text) - set(BASES):
        raise ValueError('sequence must be nonempty DNA containing only A, C, G, T')
    return text

def preserve_type(text, original):
    if isinstance(original, MutableSeq):
        return MutableSeq(text)
    if isinstance(original, Seq):
        return Seq(text)
    return text

def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)

def bounds(length, start=1, end=None):
    """Validate zero-based eligible indices with an exclusive end."""
    start = integer(start, 'start')
    end = length - 1 if end is None else integer(end, 'end')
    if not 0 <= start < end <= length:
        raise ValueError('require 0 <= start < end <= sequence length (end exclusive)')
    return (start, end)

def coding_bounds(seq, start=1, end=None):
    """Validate the coding slice independently of biological annotation filters."""
    text = dna(seq)
    (start, end) = bounds(len(text), start, end)
    if (end - start) % 3:
        raise ValueError('coding region length must be a multiple of three; specify CDS coordinates')
    return (text, start, end)

def codon_start(idx, start=1):
    integer(idx, 'idx_target')
    integer(start, 'cds_start')
    if idx < start:
        raise ValueError('mutation is outside the coding region')
    return start + 3 * ((idx - start) // 3)

def get_codon(seq, idx, start=1, end=None):
    """Return the complete codon containing a raw mutation target."""
    text = dna(seq)
    end = len(text) - 1 if end is None else end
    first = codon_start(idx, start)
    if idx >= end or first + 3 > end:
        raise ValueError('mutation is outside a complete coding codon')
    return text[first:first + 3]

def amino_acid(codon, three=False):
    text = dna(codon)
    if len(text) != 3:
        raise ValueError('codon must contain exactly three bases')
    aa = str(Seq(text).translate())
    if not three:
        return aa
    from Bio.Data.IUPACData import protein_letters_1to3
    return '*' if aa == '*' else protein_letters_1to3[aa]

def synonymous(codon1, codon2):
    return amino_acid(codon1) == amino_acid(codon2)

def parse_motif(motif):
    if not isinstance(motif, str):
        raise TypeError('mutation type must be a string')
    match = PATTERN.fullmatch(motif.upper())
    if not match:
        raise ValueError('mutation type must look like T[C>T]A or NN[C>T]')
    (left, ref, alt, right) = match.groups()
    if ref == alt:
        raise ValueError('a substitution must change the base')
    return (left, ref, alt, right)

def expand_motif(motif):
    parse_motif(motif)
    return [''.join(t) for t in itertools.product(*(IUPAC.get(c, c) for c in motif.upper()))]

def reverse_motif(motif):
    (left, ref, alt, right) = parse_motif(motif)
    return f'{Seq(right).reverse_complement()}[{Seq(ref).complement()}>{Seq(alt).complement()}]{Seq(left).reverse_complement()}'

def motif_target(motif):
    (left, ref, _, right) = parse_motif(motif)
    return left + ref + right

def motif_result(motif):
    (left, _, alt, right) = parse_motif(motif)
    return left + alt + right

def motif_indices(seq, motif):
    text = dna(seq)
    target = motif_target(motif)
    pattern = ''.join(('[' + IUPAC.get(c, c) + ']' for c in target))
    return [m.start() for m in re.finditer('(?=' + pattern + ')', text)]

def rng_from(rng=None, seed=None):
    if rng is not None and seed is not None:
        raise ValueError('provide rng or seed, not both')
    if rng is not None and (not isinstance(rng, np.random.Generator)):
        raise TypeError('rng must be numpy.random.Generator')
    return np.random.default_rng(seed) if rng is None else rng

def event(seq, idx, alt, motif, cds_start):
    """Construct one substitution with CDS-relative HGVS coordinates."""
    text = dna(seq)
    before = text[idx]
    if before == alt:
        raise ValueError('event must be a substitution')
    new = text[:idx] + alt + text[idx + 1:]
    return MutInfo(preserve_type(new, seq), idx, f'c.{idx - cds_start + 1}{before}>{alt}', motif)

def no_event(seq):
    return MutInfo(preserve_type(dna(seq), seq), -1, '', '', 1, 'no_available_event')

def random_mutation(seq, *, start=1, end=None, rng=None, seed=None, cds_start=None):
    text = dna(seq)
    (start, end) = bounds(len(text), start, end)
    generator = rng_from(rng, seed)
    idx = int(generator.integers(start, end))
    alt = str(generator.choice([b for b in BASES if b != text[idx]]))
    return event(seq, idx, alt, f'[{text[idx]}>{alt}]', start if cds_start is None else cds_start)

def motif_mutation(seq, motifs, weights=None, *, strand='both', strand_bias=0.5, start=1, end=None, rng=None, seed=None, cds_start=None):
    """Sample from all weighted, overlapping motif opportunities."""
    text = dna(seq)
    (start, end) = bounds(len(text), start, end)
    if strand not in ('both', 'single'):
        raise ValueError("strand must be 'single' or 'both'")
    bias = float(strand_bias)
    if not np.isfinite(bias) or not 0 <= bias <= 1:
        raise ValueError('strand_bias must be between zero and one')
    motifs = list(motifs)
    if not motifs:
        raise ValueError('at least one mutation type is required')
    weights = np.ones(len(motifs)) if weights is None else np.asarray(list(weights), dtype=float)
    if weights.shape != (len(motifs),) or not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError('motif weights must be finite, nonnegative, and match the motif count')
    if weights.sum() <= 0:
        raise ValueError('motif weights must have positive total mass')
    candidates = []
    for (motif, weight) in zip(motifs, weights):
        for simple in expand_motif(motif):
            specs = [(simple, 1.0)] if strand == 'single' else [(simple, bias), (reverse_motif(simple), 1 - bias)]
            for (oriented, split) in specs:
                (left, _, alt, _) = parse_motif(oriented)
                for first in motif_indices(text, oriented):
                    idx = first + len(left)
                    if start <= idx < end and weight * split > 0:
                        candidates.append((idx, alt, oriented, weight * split))
    if not candidates:
        return no_event(seq)
    p = np.array([c[3] for c in candidates], dtype=float)
    if not np.isfinite(p.sum()):
        raise ValueError('motif weights overflowed')
    selected = candidates[int(rng_from(rng, seed).choice(len(candidates), p=p / p.sum()))]
    return event(seq, selected[0], selected[1], selected[2], start if cds_start is None else cds_start)

def signature_path(version=3.3, genome_ref='GRCh37'):
    version = str(version)
    if version == '3.3':
        version = '3.3.1'
    if not re.fullmatch('[0-9]+(?:\\.[0-9]+)*', version) or not re.fullmatch('[A-Za-z0-9_]+', str(genome_ref)):
        raise ValueError('invalid signature version or reference')
    return Path(__file__).parent / 'data' / f'COSMIC_v{version}_SBS_{genome_ref}.txt'

@lru_cache(maxsize=128)
def _load_signature(path, column, bias, mtime, size):
    matrix = np.zeros((4, 64))
    seen = set()
    with open(path, newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError('signature TSV must have channel and profile columns')
        key = reader.fieldnames[0]
        column = reader.fieldnames[1] if column is None else column
        if column not in reader.fieldnames:
            raise ValueError(f'signature column not found: {column}')
        for row in reader:
            label = row[key]
            (left, ref, alt, right) = parse_motif(label)
            if len(left) != 1 or len(right) != 1 or ref not in 'CT' or any((b not in BASES for b in left + right)):
                raise ValueError('signature channels must be canonical SBS96 labels')
            if label in seen:
                raise ValueError('duplicate signature channel')
            seen.add(label)
            weight = float(row[column])
            if not np.isfinite(weight) or weight < 0:
                raise ValueError('signature weights must be finite and nonnegative')
            tri = left + ref + right
            matrix[BASES.index(alt), TRIPLET_INDEX[tri]] = bias * weight
            matrix[BASES.index(str(Seq(alt).complement())), TRIPLET_INDEX[str(Seq(tri).reverse_complement())]] = (1 - bias) * weight
    if len(seen) != 96 or matrix.sum() <= 0 or (not np.isfinite(matrix.sum())):
        raise ValueError('signature requires all 96 channels with positive finite total mass')
    matrix.setflags(write=False)
    return matrix

def full_signature(path, column=None, strand_bias=0.5):
    """Load a validated SBS96 profile; return a writable copy of cached weights."""
    bias = float(strand_bias)
    if not np.isfinite(bias) or not 0 <= bias <= 1:
        raise ValueError('strand_bias must be between zero and one')
    path = Path(path).resolve()
    stat = path.stat()
    return _load_signature(str(path), column, bias, stat.st_mtime_ns, stat.st_size).copy()

def probability_matrix(seq, signature, start=1, end=None):
    """Normalize event weights; return finite zeros when no event is possible."""
    text = dna(seq)
    (start, end) = bounds(len(text), start, end)
    matrix = np.asarray(signature, dtype=float)
    if matrix.shape != (4, 64) or not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError('signature matrix must be finite, nonnegative, and have shape (4,64)')
    for (j, tri) in enumerate(TRIPLETS):
        if matrix[BASES.index(tri[1]), j] != 0:
            raise ValueError('signature matrix must not assign weight to unchanged bases')
    weights = np.zeros((4, len(text)))
    for idx in range(max(1, start), min(len(text) - 1, end)):
        weights[:, idx] = matrix[:, TRIPLET_INDEX[text[idx - 1:idx + 2]]]
    total = weights.sum()
    if not np.isfinite(total):
        raise ValueError('signature weights overflowed')
    return weights / total if total > 0 else weights

def signature_mutation(seq, *, mutational_signature='SBS1', cosmic_version=3.3, genome_ref='GRCh37', custom_signature_path=None, column=None, strand_bias=0.5, start=1, end=None, rng=None, seed=None, cds_start=None):
    """Apply SBS weights to current sequence opportunities and sample one event."""
    text = dna(seq)
    (start, end) = bounds(len(text), start, end)
    if mutational_signature == 'custom':
        if custom_signature_path is None:
            raise ValueError('custom_signature_path is required')
        path = custom_signature_path
    else:
        path = signature_path(cosmic_version, genome_ref)
        column = mutational_signature
    matrix = full_signature(path, column, strand_bias)
    p = probability_matrix(text, matrix, start, end)
    if p.sum() == 0:
        return no_event(seq)
    flat = int(rng_from(rng, seed).choice(p.size, p=p.flatten('F')))
    (idx, alt) = (flat // 4, BASES[flat % 4])
    motif = f'{text[idx - 1]}[{text[idx]}>{alt}]{text[idx + 1]}'
    if text[idx] not in 'CT':
        motif = reverse_motif(motif)
    return event(seq, idx, alt, motif, start if cds_start is None else cds_start)
