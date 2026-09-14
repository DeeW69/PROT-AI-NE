"""Composite fitness function tying Couches 3 and 4 together."""

from __future__ import annotations

from .codon_tables import load_codon_table

DEFAULT_WEIGHTS = {"structure": 0.5, "codon_usage": 0.5}


def codon_usage_score(seq: str, organism: str = "human") -> float:
    """Cheap, always-on proxy for translational efficiency.

    Average of each codon's frequency relative to the most frequent
    synonymous codon for its amino acid (a simplified Codon Adaptation
    Index). 1.0 means every codon used is the organism's most frequent one.
    """
    table = load_codon_table(organism)
    codons = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    if not codons:
        return 0.0
    relative_freqs = [table.relative_frequency(codon) for codon in codons]
    return sum(relative_freqs) / len(relative_freqs)


def compute_fitness(
    seq: str,
    structure_score: float | None,
    codon_usage_score: float,
    weights: dict[str, float] | None = None,
) -> float:
    """Weighted composite fitness score.

    If ``structure_score`` is ``None`` (Couche 4 was skipped for this
    candidate to save compute, see ``genetic_engine``), its weight is
    redistributed onto ``codon_usage_score`` rather than assuming a value.
    """
    weights = weights or DEFAULT_WEIGHTS
    if structure_score is None:
        return codon_usage_score
    return weights["structure"] * structure_score + weights["codon_usage"] * codon_usage_score
