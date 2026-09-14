import random

import pytest
from Bio.Seq import Seq

from protaine.genetic_engine import GeneticOptimizer, _split_codons
from protaine.reverse_translate import reverse_translate

PROTEIN = "MKVLATGCWERTMKVLATGCWERT"


def _translate(mrna: str) -> str:
    return str(Seq(mrna.replace("U", "T")).translate(to_stop=True))


def _make_optimizer(**overrides) -> GeneticOptimizer:
    kwargs = dict(protein_seq=PROTEIN, population_size=8, rng_seed=7)
    kwargs.update(overrides)
    return GeneticOptimizer(**kwargs)


def test_mutation_never_changes_encoded_protein():
    optimizer = _make_optimizer(mutation_rate=1.0)
    seq = reverse_translate(PROTEIN, strategy="most_frequent")
    mutated = optimizer.mutate(seq)
    assert _translate(mutated) == PROTEIN
    assert mutated != seq  # mutation_rate=1.0 with several synonymous codons available


def test_crossover_preserves_protein_and_uses_only_parent_codons():
    optimizer = _make_optimizer()
    rng = random.Random(0)
    parents = [
        reverse_translate(PROTEIN, strategy="weighted_sample", rng=rng) for _ in range(3)
    ]
    child = optimizer.crossover(parents)
    assert _translate(child) == PROTEIN

    parent_codon_sets = [set(_split_codons(p)[i] for p in parents) for i in range(len(_split_codons(child)))]
    child_codons = _split_codons(child)
    for i, codon in enumerate(child_codons):
        assert codon in parent_codon_sets[i]


def test_crossover_requires_at_least_two_parents():
    optimizer = _make_optimizer()
    with pytest.raises(ValueError):
        optimizer.crossover([reverse_translate(PROTEIN)])


def test_initial_population_all_encode_target_protein():
    optimizer = _make_optimizer()
    for seq in optimizer.population:
        assert _translate(seq) == PROTEIN


def test_run_returns_valid_ranked_candidates():
    pytest.importorskip("RNA")
    optimizer = _make_optimizer(population_size=6, structure_eval_fraction=0.5)
    finalists = optimizer.run(generations=2, final_top_n=3)

    assert len(finalists) == 3
    assert all(_translate(c.sequence) == PROTEIN for c in finalists)
    assert all(c.structure is not None for c in finalists)
    fitness_values = [c.fitness for c in finalists]
    assert fitness_values == sorted(fitness_values, reverse=True)
    assert len(optimizer.history) == 2
