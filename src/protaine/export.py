"""Export of the final top-N candidates (Couche 6 - export half).

Writes two sibling files:
- a FASTA file with one record per candidate (sequences only);
- a JSON file with per-candidate structural metadata (dot-bracket, MFE) and
  detailed scores, indexed by rank - this is what the viewer (fornac for the
  2D mode, later 3Dmol.js/NGL.js for the 3D mode) will consume.

Only ever called on the final top-N finalists, never on the full GA
population: recomputing the exact secondary structure here (rather than
reusing the mid-search heuristic score from ``genetic_engine``) is cheap at
this scale (tens of candidates), unlike doing it for an entire population
every generation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .structure import predict_structure


@dataclass
class ExportResult:
    fasta_path: Path
    structures_path: Path


def _default_structures_path(fasta_path: Path) -> Path:
    return fasta_path.with_name(fasta_path.stem + "_structures.json")


def export_top_n(
    candidates: list,
    record_id: str,
    fasta_path: str | Path,
    structures_path: str | Path | None = None,
    n: int | None = None,
    cai_scores: list[float] | None = None,
) -> ExportResult:
    """Write the FASTA and structural-metadata JSON for the top candidates.

    ``candidates`` must already be ranked best-first (as returned by
    ``GeneticOptimizer.run``); only the first ``n`` are exported (default:
    all of them, i.e. exactly the finalists the caller already narrowed down
    to). ``cai_scores``, if given, must be aligned with ``candidates`` (same
    order) - see ``verification.compute_cai``.
    """
    fasta_path = Path(fasta_path)
    structures_path = Path(structures_path) if structures_path else _default_structures_path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    structures_path.parent.mkdir(parents=True, exist_ok=True)

    selected = candidates[:n] if n is not None else candidates
    cai_by_rank = list(cai_scores[: len(selected)]) if cai_scores else [None] * len(selected)

    structures = []
    with fasta_path.open("w") as fh:
        for rank, (candidate, cai) in enumerate(zip(selected, cai_by_rank), start=1):
            header = f">{record_id}_candidate_{rank} fitness={candidate.fitness:.4f}"
            if cai is not None:
                header += f" cai={cai:.4f}"
            fh.write(header + "\n")
            fh.write(candidate.sequence + "\n")

            fold = predict_structure(candidate.sequence)
            structures.append(
                {
                    "rank": rank,
                    "sequence": candidate.sequence,
                    "dot_bracket": fold["structure"],
                    "mfe": fold["mfe"],
                    "length": fold["length"],
                    "scores": {
                        "fitness": candidate.fitness,
                        "codon_usage": candidate.codon_usage,
                        "structure": candidate.structure,
                        "cai": cai,
                    },
                }
            )

    with structures_path.open("w") as fh:
        json.dump({"protein_id": record_id, "candidates": structures}, fh, indent=2)

    return ExportResult(fasta_path=fasta_path, structures_path=structures_path)
