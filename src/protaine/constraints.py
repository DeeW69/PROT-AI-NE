"""Constraint-based optimization on top of the reverse-translated mRNA (Couche 2).

Wraps DNAChisel, which already implements a mature constraint-satisfaction
engine for exactly this class of problem (CDS-preserving sequence design).
We do not reimplement that engine here, only the domain-specific defaults
(GC window, forbidden patterns, homopolymers, Kozak context) relevant to
mRNA therapeutic design.
"""

from __future__ import annotations

from dataclasses import dataclass, field

DEFAULT_RESTRICTION_SITES = [
    "EcoRI_site",
    "BamHI_site",
    "HindIII_site",
    "NotI_site",
    "XhoI_site",
]
KOZAK_CONTEXT = "GCCACC"  # canonical strong Kozak context placed upstream of the start codon


def _rna_to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


@dataclass
class ConstraintReport:
    sequence: str  # optimized mRNA (RNA alphabet), including Kozak context if requested
    gc_content: float
    constraints_ok: bool
    objectives_log: str = ""
    warnings: list[str] = field(default_factory=list)


def _ensure_kozak_plus4_g(dna_cds: str) -> str:
    """Best-effort: if a synonymous codon for residue 2 starts with G, use it.

    A strong Kozak context also wants a G right after the start codon
    (position +4). We only apply this when it doesn't change the encoded
    amino acid, and we don't fail if no such synonym exists.
    """
    from .codon_tables import load_codon_table

    if len(dna_cds) < 6:
        return dna_cds
    second_codon = dna_cds[3:6]
    table = load_codon_table("human")
    # Find the amino acid encoded by the current second codon, then look for
    # a synonymous codon starting with G.
    from Bio.Seq import Seq

    amino_acid = str(Seq(second_codon).translate())
    try:
        options = table.synonymous_codons(amino_acid)
    except ValueError:
        return dna_cds
    g_options = [o for o in options if _rna_to_dna(o.codon).startswith("G")]
    if not g_options:
        return dna_cds
    best_g_codon = _rna_to_dna(g_options[0].codon)
    return dna_cds[:3] + best_g_codon + dna_cds[6:]


def optimize_constraints(
    seq: str,
    gc_range: tuple[float, float] = (0.4, 0.6),
    avoid_patterns: list[str] | None = None,
    max_homopolymer_length: int = 4,
    add_kozak_context: bool = True,
    gc_window: int = 100,
) -> ConstraintReport:
    """Optimize an mRNA CDS under a default set of hard constraints.

    The amino acid sequence encoded by ``seq`` is always preserved (only
    synonymous codon substitutions are made). Returns a small,
    constraint-satisfying sequence (a single representative of the "valid
    neighborhood" Couche 3 will explore further), not an enumeration.
    """
    from dnachisel import (
        AvoidHairpins,
        AvoidPattern,
        CodonOptimize,
        DnaOptimizationProblem,
        EnforceGCContent,
        EnforceTranslation,
        NoSolutionError,
    )

    avoid_patterns = avoid_patterns or DEFAULT_RESTRICTION_SITES
    dna_cds = _rna_to_dna(seq)

    constraints = [
        EnforceTranslation(),
        EnforceGCContent(mini=gc_range[0], maxi=gc_range[1], window=gc_window),
        AvoidHairpins(stem_size=20, hairpin_window=200),
    ]
    for nucleotide in "ATGC":
        # forbid homopolymers longer than max_homopolymer_length
        constraints.append(AvoidPattern(f"{max_homopolymer_length + 1}x{nucleotide}"))
    for site in avoid_patterns:
        constraints.append(AvoidPattern(site))

    problem = DnaOptimizationProblem(
        sequence=dna_cds,
        constraints=constraints,
        objectives=[CodonOptimize(species="h_sapiens")],
        logger=None,
    )

    warnings: list[str] = []
    try:
        problem.resolve_constraints()
    except NoSolutionError as exc:
        warnings.append(f"Some constraints could not be fully satisfied: {exc}")
    problem.optimize()

    optimized_dna = str(problem.sequence)

    if add_kozak_context:
        optimized_dna = _ensure_kozak_plus4_g(optimized_dna)
        optimized_dna = KOZAK_CONTEXT + optimized_dna
        warnings.append(
            "Kozak context prepended upstream of the start codon; this is 5'UTR "
            "sequence, not part of the encoded CDS."
        )

    gc_content = (optimized_dna.count("G") + optimized_dna.count("C")) / len(optimized_dna)

    return ConstraintReport(
        sequence=_dna_to_rna(optimized_dna),
        gc_content=gc_content,
        constraints_ok=problem.all_constraints_pass(),
        objectives_log=problem.objectives_text_summary(),
        warnings=warnings,
    )
