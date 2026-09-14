"""Weighted reverse-translation: protein -> mRNA (Couche 1).

Deterministic, near-instantaneous. This module never searches: for a
protein of length n it produces exactly one (or, in "weighted_sample" mode,
one randomly drawn) sequence out of the ~6^n theoretically possible mRNAs.
It is the starting point handed to Couche 2 (constraints) and Couche 3
(genetic search), never an attempt to explore that space itself.
"""

from __future__ import annotations

import random

from .codon_tables import STOP_SYMBOL, load_codon_table

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _clean_protein(protein_seq: str) -> str:
    seq = protein_seq.strip().upper()
    seq = seq.rstrip(STOP_SYMBOL)  # trailing '*' is a stop marker, not a residue
    invalid = set(seq) - VALID_AMINO_ACIDS
    if invalid:
        raise ValueError(
            f"Protein sequence contains non-standard amino acid codes: {sorted(invalid)}"
        )
    if not seq:
        raise ValueError("Protein sequence is empty")
    return seq


def reverse_translate(
    protein_seq: str,
    organism: str = "human",
    strategy: str = "most_frequent",
    include_stop: bool = True,
    rng: random.Random | None = None,
) -> str:
    """Reverse-translate a protein into an mRNA (5'->3', RNA alphabet).

    Args:
        protein_seq: protein sequence (one-letter codes, stop '*' optional).
        organism: codon usage table to use (see ``codon_tables``).
        strategy: "most_frequent" picks the single most common codon per
            amino acid (fully deterministic). "weighted_sample" draws a
            codon per amino acid with probability proportional to its usage
            frequency in the organism's table.
        include_stop: append the organism's most frequent stop codon.
        rng: optional ``random.Random`` instance (for reproducible sampling).

    Returns:
        The mRNA sequence as a string using RNA letters (A, C, G, U).
    """
    protein = _clean_protein(protein_seq)
    table = load_codon_table(organism)
    rng = rng or random

    codons: list[str] = []
    for amino_acid in protein:
        options = table.synonymous_codons(amino_acid)
        if strategy == "most_frequent":
            codons.append(options[0].codon)
        elif strategy == "weighted_sample":
            weights = [o.per_thousand for o in options]
            codons.append(rng.choices(options, weights=weights, k=1)[0].codon)
        else:
            raise ValueError(f"Unknown strategy: {strategy!r}")

    if include_stop:
        stop_options = table.synonymous_codons(STOP_SYMBOL)
        codons.append(stop_options[0].codon)

    return "".join(codons)
