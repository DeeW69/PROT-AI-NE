import subprocess

import pytest

from protaine.genetic_engine import Candidate
from protaine.prediction_3d import predict_3d_structure, predict_top_n_3d

SEQ = "AUGAAAGUGUAA"


def _make_fake_rhofold_home(tmp_path):
    home = tmp_path / "RhoFold"
    home.mkdir()
    (home / "inference.py").write_text("# fake inference.py\n")
    (home / "pretrained").mkdir()
    return home


def test_returns_none_when_rhofold_not_installed(tmp_path, monkeypatch):
    monkeypatch.delenv("RHOFOLD_HOME", raising=False)
    monkeypatch.chdir(tmp_path)  # no ./RhoFold here either
    result = predict_3d_structure(SEQ, tmp_path / "out")
    assert result is None


def test_returns_pdb_path_on_success(tmp_path, monkeypatch):
    home = _make_fake_rhofold_home(tmp_path)
    output_dir = tmp_path / "out"

    def fake_run(cmd, cwd, capture_output, text, timeout):
        # Simulate RhoFold+ writing its output file into --output_dir.
        (output_dir / "unrelaxed_model.pdb").write_text("HEADER\nEND\n")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr("protaine.prediction_3d.subprocess.run", fake_run)

    result = predict_3d_structure(SEQ, output_dir, rhofold_home=home)

    assert result == output_dir / "unrelaxed_model.pdb"
    assert result.read_text().startswith("HEADER")
    # The input FASTA that was written for RhoFold+ to consume.
    assert (output_dir / "candidate.fasta").read_text().splitlines()[1] == SEQ


def test_returns_none_on_nonzero_exit_code(tmp_path, monkeypatch):
    home = _make_fake_rhofold_home(tmp_path)
    output_dir = tmp_path / "out"

    def fake_run(cmd, cwd, capture_output, text, timeout):
        return subprocess.CompletedProcess(cmd, returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr("protaine.prediction_3d.subprocess.run", fake_run)

    assert predict_3d_structure(SEQ, output_dir, rhofold_home=home) is None


def test_returns_none_on_timeout(tmp_path, monkeypatch):
    home = _make_fake_rhofold_home(tmp_path)
    output_dir = tmp_path / "out"

    def fake_run(cmd, cwd, capture_output, text, timeout):
        raise subprocess.TimeoutExpired(cmd, timeout)

    monkeypatch.setattr("protaine.prediction_3d.subprocess.run", fake_run)

    assert predict_3d_structure(SEQ, output_dir, rhofold_home=home, timeout=5) is None


def test_returns_none_when_success_but_pdb_missing(tmp_path, monkeypatch):
    home = _make_fake_rhofold_home(tmp_path)
    output_dir = tmp_path / "out"

    def fake_run(cmd, cwd, capture_output, text, timeout):
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr("protaine.prediction_3d.subprocess.run", fake_run)

    assert predict_3d_structure(SEQ, output_dir, rhofold_home=home) is None


def test_returns_none_when_subprocess_cannot_launch(tmp_path, monkeypatch):
    home = _make_fake_rhofold_home(tmp_path)

    def fake_run(cmd, cwd, capture_output, text, timeout):
        raise OSError("no such interpreter")

    monkeypatch.setattr("protaine.prediction_3d.subprocess.run", fake_run)

    assert predict_3d_structure(SEQ, tmp_path / "out", rhofold_home=home) is None


CANDIDATES = [
    Candidate(sequence="AUGAAAGUGUAA", codon_usage=0.9, structure=0.7, fitness=0.8),
    Candidate(sequence="AUGAAGGUCUGA", codon_usage=0.6, structure=0.5, fitness=0.55),
    Candidate(sequence="AUGAAAAUUUGA", codon_usage=0.4, structure=0.3, fitness=0.35),
]


def test_predict_top_n_3d_isolates_per_candidate_failures(tmp_path, monkeypatch):
    def fake_predict(rna_seq, out_dir, sequence_id, timeout, device, rhofold_home):
        out_dir = __import__("pathlib").Path(out_dir)
        if sequence_id == "candidate_2":
            return None  # simulate a crash/timeout on just this one
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "unrelaxed_model.pdb"
        path.write_text("HEADER\n" + sequence_id + "\nEND\n")
        return path

    monkeypatch.setattr("protaine.prediction_3d.predict_3d_structure", fake_predict)

    output_dir = tmp_path / "top_structures3d"
    results = predict_top_n_3d(CANDIDATES, output_dir)

    assert results[1] == output_dir / "candidate_1.pdb"
    assert results[1].exists()
    assert results[2] is None
    assert not (output_dir / "candidate_2.pdb").exists()
    assert results[3] == output_dir / "candidate_3.pdb"
    assert results[3].exists()

    # Per-candidate scratch directories must not leak into the final output.
    assert not (output_dir / "_rhofold_work_1").exists()
    assert not (output_dir / "_rhofold_work_2").exists()
    assert not (output_dir / "_rhofold_work_3").exists()


def test_predict_top_n_3d_respects_n_slicing(tmp_path, monkeypatch):
    calls = []

    def fake_predict(rna_seq, out_dir, sequence_id, timeout, device, rhofold_home):
        calls.append(sequence_id)
        return None

    monkeypatch.setattr("protaine.prediction_3d.predict_3d_structure", fake_predict)

    results = predict_top_n_3d(CANDIDATES, tmp_path / "out", n=1)

    assert calls == ["candidate_1"]
    assert list(results.keys()) == [1]
