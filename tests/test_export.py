import json

import pytest

pytest.importorskip("RNA")

from protaine.export import export_top_n
from protaine.genetic_engine import Candidate

SEQ_A = "AUGAAAGUGCUGGCCACCGGCUGCUGGGAACGCACCAUGAAAGUGCUGGCCACCGGCUGCUGGGAACGCACC"
SEQ_B = "AUGAAGGUCCUCGCCACAGGGUGUUGGGAGCGGACUAUGAAGGUCCUCGCCACAGGGUGUUGGGAGCGGACU"

CANDIDATES = [
    Candidate(sequence=SEQ_A, codon_usage=0.9, structure=0.7, fitness=0.8),
    Candidate(sequence=SEQ_B, codon_usage=0.6, structure=0.5, fitness=0.55),
]


def test_export_writes_fasta_with_one_record_per_candidate(tmp_path):
    result = export_top_n(CANDIDATES, "myprotein", tmp_path / "top.fasta")

    content = result.fasta_path.read_text()
    assert content.count(">") == 2
    assert "myprotein_candidate_1" in content
    assert "myprotein_candidate_2" in content
    assert SEQ_A in content
    assert SEQ_B in content


def test_export_default_structures_path_is_sibling_json(tmp_path):
    result = export_top_n(CANDIDATES, "myprotein", tmp_path / "top.fasta")
    assert result.structures_path == tmp_path / "top_structures.json"
    assert result.structures_path.exists()


def test_export_structures_json_has_one_entry_per_candidate_with_dot_bracket(tmp_path):
    result = export_top_n(CANDIDATES, "myprotein", tmp_path / "top.fasta")

    data = json.loads(result.structures_path.read_text())
    assert data["protein_id"] == "myprotein"
    assert len(data["candidates"]) == 2

    first = data["candidates"][0]
    assert first["rank"] == 1
    assert set(first["dot_bracket"]) <= set(".()")
    assert isinstance(first["mfe"], float)
    assert first["scores"]["fitness"] == 0.8


def test_export_respects_n_slicing(tmp_path):
    result = export_top_n(CANDIDATES, "myprotein", tmp_path / "top.fasta", n=1)

    fasta_content = result.fasta_path.read_text()
    assert fasta_content.count(">") == 1

    data = json.loads(result.structures_path.read_text())
    assert len(data["candidates"]) == 1


def test_export_includes_cai_scores_when_given(tmp_path):
    result = export_top_n(
        CANDIDATES, "myprotein", tmp_path / "top.fasta", cai_scores=[0.87, 0.62]
    )

    fasta_content = result.fasta_path.read_text()
    assert "cai=0.8700" in fasta_content
    assert "cai=0.6200" in fasta_content

    data = json.loads(result.structures_path.read_text())
    assert data["candidates"][0]["scores"]["cai"] == 0.87
    assert data["candidates"][1]["scores"]["cai"] == 0.62


def test_export_omits_cai_from_header_when_not_given(tmp_path):
    result = export_top_n(CANDIDATES, "myprotein", tmp_path / "top.fasta")

    fasta_content = result.fasta_path.read_text()
    assert "cai=" not in fasta_content

    data = json.loads(result.structures_path.read_text())
    assert data["candidates"][0]["scores"]["cai"] is None


def test_export_creates_missing_output_directory(tmp_path):
    nested = tmp_path / "nested" / "dir"
    result = export_top_n(CANDIDATES, "myprotein", nested / "top.fasta")

    assert result.fasta_path.exists()
    assert result.structures_path.exists()


def test_export_honors_explicit_structures_path(tmp_path):
    explicit = tmp_path / "custom_structures.json"
    result = export_top_n(
        CANDIDATES, "myprotein", tmp_path / "top.fasta", structures_path=explicit
    )
    assert result.structures_path == explicit
    assert explicit.exists()
