import pytest

pytest.importorskip("RNA")

from protaine.structure import predict_structure, structure_score

SEQ = "AUGAAAGUGCUGGCCACCGGCUGCUGGGAACGCACCAUGAAAGUGCUGGCCACCGGCUGCUGGGAACGCACC"


def test_predict_structure_returns_matching_length_dot_bracket():
    result = predict_structure(SEQ)
    assert len(result["structure"]) == len(SEQ)
    assert set(result["structure"]) <= set(".()")
    assert isinstance(result["mfe"], float)


def test_structure_score_is_bounded():
    result = predict_structure(SEQ)
    score = structure_score(result)
    assert 0.0 <= score <= 1.0


def test_fully_unpaired_structure_scores_higher_than_fully_paired():
    unpaired = {"structure": "." * 40, "mfe": 0.0, "length": 40}
    paired = {"structure": "(" * 20 + ")" * 20, "mfe": -25.0, "length": 40}
    assert structure_score(unpaired) > structure_score(paired)
