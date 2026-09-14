import random

import pytest
from Bio.Seq import Seq

from protaine.reverse_translate import reverse_translate


def test_length_includes_stop_codon():
    protein = "MKV"
    mrna = reverse_translate(protein, strategy="most_frequent")
    assert len(mrna) == 3 * (len(protein) + 1)


def test_length_without_stop_codon():
    protein = "MKV"
    mrna = reverse_translate(protein, strategy="most_frequent", include_stop=False)
    assert len(mrna) == 3 * len(protein)


def test_most_frequent_translates_back_to_original_protein():
    protein = "MKVLAT"
    mrna = reverse_translate(protein, strategy="most_frequent")
    dna = mrna.replace("U", "T")
    translated = str(Seq(dna).translate(to_stop=True))
    assert translated == protein


def test_most_frequent_is_deterministic():
    protein = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPVPSQAMDDLMLSPDDIEQWFTEDPGPD"
    first = reverse_translate(protein, strategy="most_frequent")
    second = reverse_translate(protein, strategy="most_frequent")
    assert first == second


def test_weighted_sample_translates_back_to_original_protein():
    protein = "MKVLATGCWERT"
    rng = random.Random(42)
    mrna = reverse_translate(protein, strategy="weighted_sample", rng=rng)
    dna = mrna.replace("U", "T")
    translated = str(Seq(dna).translate(to_stop=True))
    assert translated == protein


def test_weighted_sample_produces_variation():
    protein = "MKVLATGCWERTMKVLATGCWERT"
    rng = random.Random(1)
    samples = {reverse_translate(protein, strategy="weighted_sample", rng=rng) for _ in range(20)}
    assert len(samples) > 1


def test_rejects_invalid_amino_acid():
    with pytest.raises(ValueError):
        reverse_translate("MKVXZ")


def test_rejects_empty_sequence():
    with pytest.raises(ValueError):
        reverse_translate("")


def test_strips_trailing_stop_symbol():
    mrna_with_star = reverse_translate("MKV*", strategy="most_frequent")
    mrna_without_star = reverse_translate("MKV", strategy="most_frequent")
    assert mrna_with_star == mrna_without_star
