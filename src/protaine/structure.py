"""Secondary-structure prediction and scoring (Couche 4).

Uses ViennaRNA's Python bindings (``RNA.fold``), not the pure-Python
Nussinov implementation from Gene Explorer: Nussinov's O(n^3) blows up well
before the ~900 nt typical of a 300-residue CDS, whereas ViennaRNA's
partition-function/MFE algorithms handle sequences of that length routinely.
"""

from __future__ import annotations

try:
    import RNA
except ImportError as exc:  # pragma: no cover - depends on local install
    RNA = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def _require_vienna() -> None:
    if RNA is None:
        raise ImportError(
            "The ViennaRNA Python bindings are required for structure prediction. "
            "Install with `pip install ViennaRNA` (Linux/macOS wheels on PyPI) or, "
            "on Windows, via conda-forge (`conda install -c bioconda viennarna`) "
            "or WSL."
        ) from _IMPORT_ERROR


def predict_structure(seq: str) -> dict:
    """Predict the MFE secondary structure of an RNA sequence.

    Returns a dict with the dot-bracket ``structure`` and the minimum free
    energy ``mfe`` (kcal/mol).
    """
    _require_vienna()
    dot_bracket, mfe = RNA.fold(seq.upper().replace("T", "U"))
    return {"structure": dot_bracket, "mfe": mfe, "length": len(seq)}


def structure_score(structure_result: dict, five_prime_window: int = 50) -> float:
    """Heuristic score in [0, 1]; higher means "less likely to hinder translation".

    This is a simplified proxy, not a validated predictor of translational
    efficiency: it favors sequences whose 5' region (ribosome
    loading/scanning site) is largely unpaired, and penalizes very stable
    (very negative MFE per nucleotide) overall folds.
    """
    structure = structure_result["structure"]
    mfe = structure_result["mfe"]
    length = len(structure) or 1

    window = structure[: min(five_prime_window, length)]
    unpaired_5p_fraction = window.count(".") / len(window) if window else 1.0

    mfe_per_nt = mfe / length
    # Empirically, mRNA MFE per nucleotide rarely goes far below ~ -0.35 kcal/mol/nt;
    # clamp so the penalty saturates instead of dominating the score.
    stability_penalty = min(1.0, max(0.0, -mfe_per_nt / 0.35))

    return 0.7 * unpaired_5p_fraction + 0.3 * (1.0 - stability_penalty)
