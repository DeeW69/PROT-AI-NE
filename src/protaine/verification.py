"""Post-selection verification (Couche 5).

Applied only to the final top-N candidates retained after the genetic search
(Couche 3), never to the whole population - these checks are either exact-but-
cheap (back-translation) or slow/network-bound (BLAST) and would be wasteful
or infeasible to run on every candidate of every generation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .codon_tables import load_codon_table


def verify_translation(rna_seq: str, target_protein: str) -> bool:
    """Translate ``rna_seq`` and confirm it exactly matches ``target_protein``.

    Blocking check: a candidate that fails this (frameshift, premature stop,
    unintended amino acid substitution) must never be exported or reported as
    a valid result, regardless of its fitness score.
    """
    from Bio.Seq import Seq

    dna = rna_seq.upper().replace("U", "T")
    usable_len = len(dna) - len(dna) % 3
    translated = str(Seq(dna[:usable_len]).translate(to_stop=True))
    target = target_protein.strip().upper().rstrip("*")
    return translated == target


def compute_cai(rna_seq: str, organism: str = "human") -> float:
    """Codon Adaptation Index (Sharp & Li, 1987) relative to ``organism`` usage.

    Geometric mean, over all codons in ``rna_seq``, of each codon's relative
    adaptedness (its usage frequency divided by that of the most frequent
    synonymous codon for the same amino acid). 1.0 means every codon used is
    already the organism's most frequent choice for its amino acid; unknown
    codons (e.g. leftover partial codon) and the stop codon are excluded, per
    the standard CAI convention.
    """
    table = load_codon_table(organism)
    dna = rna_seq.upper().replace("T", "U")
    codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]

    weights = []
    for codon in codons:
        try:
            info = table.codon_info(codon)
        except ValueError:
            continue
        if info.amino_acid == "*":
            continue
        weights.append(table.relative_frequency(codon))

    if not weights:
        return 0.0
    log_mean = sum(math.log(w) for w in weights if w > 0) / len(weights)
    return math.exp(log_mean)


@dataclass
class OffTargetHit:
    accession: str
    description: str
    identity_pct: float
    e_value: float


def check_off_target_homology(
    rna_seq: str,
    organism: str = "human",
    e_value_threshold: float = 1e-6,
    hitlist_size: int = 5,
) -> list[OffTargetHit]:
    """Query remote NCBI BLAST (blastn/nt) for unintended matches in ``organism``.

    Informative only, never blocking: network-bound and can take from tens of
    seconds to several minutes per candidate (NCBI's shared queue), so this is
    meant to be run only on the final top-N, and callers should treat any
    exception here (network failure, service timeout) as "check unavailable"
    rather than a reason to reject the candidate.
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    dna = rna_seq.upper().replace("U", "T")
    entrez_query = f"{organism}[Organism]" if organism.lower() != "human" else "Homo sapiens[Organism]"

    result_handle = NCBIWWW.qblast(
        "blastn",
        "nt",
        dna,
        entrez_query=entrez_query,
        expect=e_value_threshold,
        hitlist_size=hitlist_size,
    )
    try:
        record = NCBIXML.read(result_handle)
    finally:
        result_handle.close()

    hits = []
    for alignment in record.alignments[:hitlist_size]:
        hsp = alignment.hsps[0]
        identity_pct = 100.0 * hsp.identities / hsp.align_length if hsp.align_length else 0.0
        hits.append(
            OffTargetHit(
                accession=alignment.accession,
                description=alignment.hit_def,
                identity_pct=identity_pct,
                e_value=hsp.expect,
            )
        )
    return hits
