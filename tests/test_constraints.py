import pytest
from Bio.Seq import Seq

pytest.importorskip("dnachisel")

from protaine.constraints import optimize_constraints
from protaine.reverse_translate import reverse_translate

PROTEIN = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPVPSQAMDDLMLSPDDIEQWFTEDPGPD"


def test_optimized_sequence_preserves_protein():
    seed = reverse_translate(PROTEIN, strategy="most_frequent", include_stop=False)
    report = optimize_constraints(seed, add_kozak_context=False)
    dna = report.sequence.replace("U", "T")
    assert str(Seq(dna).translate()) == PROTEIN


def test_gc_content_within_requested_range():
    seed = reverse_translate(PROTEIN, strategy="most_frequent", include_stop=False)
    report = optimize_constraints(seed, gc_range=(0.45, 0.55), add_kozak_context=False)
    assert 0.40 <= report.gc_content <= 0.60


def test_kozak_context_prepended_when_requested():
    seed = reverse_translate("MAKV", strategy="most_frequent", include_stop=False)
    report = optimize_constraints(seed, add_kozak_context=True)
    assert report.sequence.startswith("GCCACC")
    cds = report.sequence[len("GCCACC") :]
    assert cds.startswith("AUG")
