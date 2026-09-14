"""Codon usage tables (Couche 1).

Loads per-organism codon usage frequencies (source: Kazusa Codon Usage
Database, https://www.kazusa.or.jp/codon/) from ``data/codon_usage/*.csv``
and exposes them as synonymous-codon groups keyed by amino acid.

Codons are stored and returned in RNA notation (``U`` rather than ``T``)
since the whole pipeline reasons in terms of mRNA sequences.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "codon_usage"

STOP_SYMBOL = "*"


@dataclass(frozen=True)
class CodonUsage:
    codon: str
    amino_acid: str
    per_thousand: float
    fraction: float  # frequency relative to other codons for the same amino acid


class CodonTable:
    """Synonymous-codon groups for one organism, sorted by decreasing usage."""

    def __init__(self, entries: list[CodonUsage]):
        by_aa: dict[str, list[CodonUsage]] = {}
        for entry in entries:
            by_aa.setdefault(entry.amino_acid, []).append(entry)
        for aa, codons in by_aa.items():
            codons.sort(key=lambda c: c.per_thousand, reverse=True)
        self._by_aa = by_aa
        self._by_codon = {entry.codon: entry for entry in entries}

    def amino_acids(self) -> list[str]:
        return list(self._by_aa.keys())

    def synonymous_codons(self, amino_acid: str) -> list[CodonUsage]:
        try:
            return self._by_aa[amino_acid.upper()]
        except KeyError as exc:
            raise ValueError(f"Unknown amino acid code: {amino_acid!r}") from exc

    def best_codon(self, amino_acid: str) -> CodonUsage:
        return self.synonymous_codons(amino_acid)[0]

    def codon_info(self, codon: str) -> CodonUsage:
        try:
            return self._by_codon[codon.upper()]
        except KeyError as exc:
            raise ValueError(f"Unknown codon: {codon!r}") from exc

    def relative_frequency(self, codon: str) -> float:
        """Frequency of ``codon`` relative to the best codon for its amino acid.

        Used as a simple Codon-Adaptation-Index-like weight: 1.0 for the most
        frequent codon of an amino acid, lower for rarer synonyms.
        """
        info = self.codon_info(codon)
        best = self.best_codon(info.amino_acid)
        return info.per_thousand / best.per_thousand


def _load_csv(path: Path) -> list[CodonUsage]:
    entries: list[CodonUsage] = []
    totals: dict[str, float] = {}
    rows = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            codon = row["codon"].strip().upper()
            aa = row["amino_acid"].strip().upper()
            per_thousand = float(row["per_thousand"])
            rows.append((codon, aa, per_thousand))
            totals[aa] = totals.get(aa, 0.0) + per_thousand

    for codon, aa, per_thousand in rows:
        fraction = per_thousand / totals[aa] if totals[aa] else 0.0
        entries.append(CodonUsage(codon=codon, amino_acid=aa, per_thousand=per_thousand, fraction=fraction))
    return entries


@lru_cache(maxsize=None)
def load_codon_table(organism: str = "human") -> CodonTable:
    """Load (and cache) the codon usage table for ``organism``.

    Only ``"human"`` ships by default. Additional organisms can be added by
    dropping a ``<organism>.csv`` file (columns: codon, amino_acid,
    per_thousand) into ``data/codon_usage/``.
    """
    path = DATA_DIR / f"{organism.lower()}.csv"
    if not path.exists():
        available = sorted(p.stem for p in DATA_DIR.glob("*.csv"))
        raise FileNotFoundError(
            f"No codon usage table for organism {organism!r} (looked for {path}). "
            f"Available: {available}"
        )
    return CodonTable(_load_csv(path))
