"""RhoFold+ local 3D structure inference (Couche 6, mode 3D).

RhoFold+ (Shen et al., "Accurate RNA 3D structure prediction using a
language model-based deep learning approach", Nature Methods 2024) is not a
pip package: it must be cloned from https://github.com/ml4bio/RhoFold into
its own conda environment (Linux only per its own README - Windows users
need WSL, same constraint already documented for ViennaRNA in this
project), with its pretrained checkpoint downloaded separately. This module
therefore shells out to its documented ``inference.py`` CLI - the only
interface RhoFold+ publishes with any stability guarantee - rather than
importing its internal (research-grade, unversioned) Python classes.

Because RhoFold+ lives in a separate environment, invocation deliberately
does NOT reuse this project's own Python interpreter (``sys.executable``):
see ``RHOFOLD_HOME``/``RHOFOLD_PYTHON`` below.

Only ever called on the final top-N finalists (see ``export.export_top_n``
and ``cli.py``'s ``--predict-3d``), never inside the GA loop (Couche 3) -
inference takes seconds to minutes per sequence even on GPU, and only makes
sense once a sequence has already been selected. Every failure mode
(RhoFold+ not installed/configured, no GPU, a crashed or stalled
subprocess, a single candidate erroring out) degrades to returning ``None``
for that candidate rather than raising - the viewer already treats a
missing PDB file as "structure 3D non disponible" (see
``viewer/viewer.js``'s ``render3d``), never as a crash.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

PDB_OUTPUT_NAME = "unrelaxed_model.pdb"


class RhoFoldNotFound(RuntimeError):
    """Raised internally when no usable RhoFold+ installation can be located."""


def _rhofold_home(rhofold_home: str | Path | None = None) -> Path:
    """Resolve the local RhoFold+ checkout: explicit arg, then $RHOFOLD_HOME, then ./RhoFold."""
    candidates = []
    if rhofold_home:
        candidates.append(Path(rhofold_home))
    env_path = os.environ.get("RHOFOLD_HOME")
    if env_path:
        candidates.append(Path(env_path))
    candidates.append(Path.cwd() / "RhoFold")

    for candidate in candidates:
        if (candidate / "inference.py").is_file():
            return candidate

    raise RhoFoldNotFound(
        "RhoFold+ introuvable. Cloner https://github.com/ml4bio/RhoFold "
        "(Linux uniquement - WSL sous Windows), suivre son installation "
        "(environnement conda dedie + `python setup.py install` + "
        "telechargement du checkpoint pretrained dans pretrained/), puis "
        "exporter RHOFOLD_HOME=/chemin/vers/RhoFold."
    )


def _rhofold_python() -> str:
    """Python interpreter to run RhoFold+ with.

    RhoFold+ needs its own conda environment (specific torch/CUDA versions,
    its own dependencies) - almost certainly NOT this project's own
    virtualenv. ``$RHOFOLD_PYTHON`` should point at that environment's
    interpreter (e.g. the output of `conda run -n RhoFold which python`).
    Falling back to this process's own interpreter is a last resort that
    will typically fail with an import error inside the subprocess (caught
    and reported like any other RhoFold+ failure), not a silent success.
    """
    return os.environ.get("RHOFOLD_PYTHON", sys.executable)


def predict_3d_structure(
    rna_seq: str,
    output_dir: str | Path,
    sequence_id: str = "candidate",
    timeout: int = 300,
    device: str | None = None,
    rhofold_home: str | Path | None = None,
) -> Path | None:
    """Run RhoFold+ on a single RNA sequence; return the PDB path, or ``None``.

    Uses single-sequence mode (``--single_seq_pred True``, no MSA): per
    RhoFold+'s own README this is explicitly framed as a lower-accuracy
    "quick reference structure" mode, but building a real MSA requires
    ~900GB of sequence databases, wildly out of scope for scoring
    computationally-generated mRNA candidates one at a time.

    This function must never raise for a "prediction unavailable" reason -
    only ever print a diagnostic to stderr and return ``None`` - so callers
    can treat 3D prediction purely as a best-effort enhancement.
    """
    output_dir = Path(output_dir)
    try:
        home = _rhofold_home(rhofold_home)
    except RhoFoldNotFound as exc:
        print(f"[prediction_3d] {exc}", file=sys.stderr)
        return None

    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / f"{sequence_id}.fasta"
    fasta_path.write_text(f">{sequence_id}\n{rna_seq.upper().replace('T', 'U')}\n")

    ckpt = home / "pretrained" / "RhoFold_pretrained.pt"
    cmd = [
        _rhofold_python(),
        str(home / "inference.py"),
        "--input_fas", str(fasta_path),
        "--output_dir", str(output_dir),
        "--single_seq_pred", "True",
        "--device", device or "cpu",
    ]
    if ckpt.is_file():
        cmd += ["--ckpt", str(ckpt)]

    try:
        result = subprocess.run(cmd, cwd=home, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"[prediction_3d] RhoFold+ a depasse le timeout de {timeout}s pour '{sequence_id}'.", file=sys.stderr)
        return None
    except OSError as exc:
        print(f"[prediction_3d] Impossible de lancer RhoFold+ ({exc}).", file=sys.stderr)
        return None

    if result.returncode != 0:
        print(
            f"[prediction_3d] RhoFold+ a echoue pour '{sequence_id}' (code {result.returncode}) :\n"
            f"{result.stderr[-2000:]}",
            file=sys.stderr,
        )
        return None

    pdb_path = output_dir / PDB_OUTPUT_NAME
    if not pdb_path.is_file():
        print(f"[prediction_3d] RhoFold+ n'a pas produit '{PDB_OUTPUT_NAME}' pour '{sequence_id}'.", file=sys.stderr)
        return None

    return pdb_path


def predict_top_n_3d(
    candidates: list,
    output_dir: str | Path,
    n: int | None = None,
    timeout: int = 300,
    device: str | None = None,
    rhofold_home: str | Path | None = None,
) -> dict[int, Path | None]:
    """Predict 3D structures for the top-N finalists, one candidate at a time.

    ``candidates`` must already be ranked best-first, matching
    ``export.export_top_n``'s convention (rank = 1-based position, not a
    field on the candidate). Output files land directly in ``output_dir``
    as ``candidate_<rank>.pdb`` - the exact path
    ``viewer/viewer.js``'s ``pdbUrlFor`` already looks for next to the
    exported structures JSON (``<stem>_structures3d/candidate_<rank>.pdb``).

    Never raises: a failure on one candidate (crashed/timed-out sequence,
    or RhoFold+ missing entirely) is isolated via
    ``predict_3d_structure``'s own error handling and does not prevent the
    remaining candidates from being attempted - mirroring the
    per-candidate isolation already used for off-target BLAST checks in
    ``verification.check_off_target_homology``.
    """
    output_dir = Path(output_dir)
    selected = candidates[:n] if n is not None else candidates
    results: dict[int, Path | None] = {}

    for rank, candidate in enumerate(selected, start=1):
        work_dir = output_dir / f"_rhofold_work_{rank}"
        pdb_path = predict_3d_structure(
            candidate.sequence,
            work_dir,
            sequence_id=f"candidate_{rank}",
            timeout=timeout,
            device=device,
            rhofold_home=rhofold_home,
        )
        if pdb_path is not None:
            output_dir.mkdir(parents=True, exist_ok=True)
            final_path = output_dir / f"candidate_{rank}.pdb"
            shutil.copyfile(pdb_path, final_path)
            results[rank] = final_path
        else:
            results[rank] = None
        shutil.rmtree(work_dir, ignore_errors=True)

    return results
