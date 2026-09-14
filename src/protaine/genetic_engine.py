"""Genetic search over synonymous codon choices (Couche 3).

Every population member encodes exactly the same protein — mutation only
ever substitutes a codon for a synonymous one, crossover only ever mixes
codons position-by-position across parents that already encode that
protein — so the search never leaves the "silent" subspace consistent with
the target. It operates on the small neighborhood handed over by Couche 2,
never on the full ~6^n theoretical codon space.

Couche 4 (secondary-structure scoring) is comparatively expensive, so per
the project's design it is only recomputed for the top fraction of the
population each generation (``structure_eval_fraction``), not for everyone.
Candidates that haven't had it computed yet fall back to their codon-usage
score alone as an approximate fitness (see ``fitness.compute_fitness``).
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .codon_tables import load_codon_table
from .fitness import codon_usage_score as _codon_usage_score
from .fitness import compute_fitness
from .reverse_translate import reverse_translate
from .structure import predict_structure
from .structure import structure_score as _structure_score


def _split_codons(seq: str) -> list[str]:
    return [seq[i : i + 3] for i in range(0, len(seq), 3)]


@dataclass
class Candidate:
    sequence: str
    codon_usage: float
    structure: float | None
    fitness: float


@dataclass
class GenerationRecord:
    generation: int
    best_fitness: float
    mean_fitness: float
    best_sequence: str


class GeneticOptimizer:
    def __init__(
        self,
        protein_seq: str,
        organism: str = "human",
        population_size: int = 100,
        seed_sequences: list[str] | None = None,
        mutation_rate: float = 0.05,
        elite_fraction: float = 0.1,
        structure_eval_fraction: float = 0.2,
        parents_per_crossover: tuple[int, int] = (2, 4),
        weights: dict[str, float] | None = None,
        rng_seed: int | None = None,
    ):
        self.protein_seq = protein_seq.strip().upper().rstrip("*")
        self.organism = organism
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.elite_fraction = elite_fraction
        self.structure_eval_fraction = structure_eval_fraction
        self.parents_per_crossover = parents_per_crossover
        self.weights = weights
        self.rng = random.Random(rng_seed)
        self.table = load_codon_table(organism)

        self.history: list[GenerationRecord] = []
        self.population: list[str] = self._init_population(seed_sequences or [])

    def _init_population(self, seed_sequences: list[str]) -> list[str]:
        population = list(dict.fromkeys(seed_sequences))  # de-dup, keep order
        while len(population) < self.population_size:
            population.append(
                reverse_translate(
                    self.protein_seq,
                    organism=self.organism,
                    strategy="weighted_sample",
                    rng=self.rng,
                )
            )
        return population[: self.population_size]

    def mutate(self, seq: str) -> str:
        """Randomly substitute some codons for synonymous alternatives."""
        codons = _split_codons(seq)
        for i, codon in enumerate(codons):
            if self.rng.random() >= self.mutation_rate:
                continue
            amino_acid = self.table.codon_info(codon).amino_acid
            options = [o.codon for o in self.table.synonymous_codons(amino_acid) if o.codon != codon]
            if options:
                codons[i] = self.rng.choice(options)
        return "".join(codons)

    def crossover(self, parents: list[str]) -> str:
        """Uniform per-codon crossover across 2-4 parents.

        Every parent encodes the same protein, so picking any parent's
        codon at each position always yields a valid, same-protein child.
        """
        if len(parents) < 2:
            raise ValueError("crossover requires at least 2 parents")
        parent_codons = [_split_codons(p) for p in parents]
        n_codons = len(parent_codons[0])
        child = [self.rng.choice(parent_codons)[i] for i in range(n_codons)]
        return "".join(child)

    def _evaluate_cheap(self, seq: str) -> float:
        return _codon_usage_score(seq, organism=self.organism)

    def _evaluate_structure(self, seq: str) -> float:
        return _structure_score(predict_structure(seq))

    def _select_parents(self, ranked: list["Candidate"]) -> list[str]:
        n = self.rng.randint(*self.parents_per_crossover)
        pool = ranked[: max(n, len(ranked) // 2)] or ranked
        candidates = list(pool)
        weights = [max(c.fitness, 1e-6) for c in candidates]
        chosen: list[str] = []
        for _ in range(min(n, len(candidates))):
            picked = self.rng.choices(candidates, weights=weights, k=1)[0]
            idx = candidates.index(picked)
            chosen.append(candidates.pop(idx).sequence)
            weights.pop(idx)
        return chosen

    def _rank_population(self, population: list[str]) -> list[Candidate]:
        cheap_scores = [self._evaluate_cheap(seq) for seq in population]
        order = sorted(range(len(population)), key=lambda i: cheap_scores[i], reverse=True)

        top_k = max(1, int(len(population) * self.structure_eval_fraction))
        structure_scores: dict[int, float] = {}
        for i in order[:top_k]:
            structure_scores[i] = self._evaluate_structure(population[i])

        candidates = []
        for i, seq in enumerate(population):
            structure = structure_scores.get(i)
            fitness = compute_fitness(seq, structure, cheap_scores[i], weights=self.weights)
            candidates.append(
                Candidate(sequence=seq, codon_usage=cheap_scores[i], structure=structure, fitness=fitness)
            )
        candidates.sort(key=lambda c: c.fitness, reverse=True)
        return candidates

    def run(self, generations: int, final_top_n: int = 20) -> list[Candidate]:
        """Run the search and return the best candidates found.

        The returned candidates always have a fully computed structure
        score (a final structural pass is run on the finalists so the
        reported ranking isn't left relying on the codon-usage
        approximation used mid-search).
        """
        ranked = self._rank_population(self.population)
        for gen in range(generations):
            n_elite = max(1, int(len(ranked) * self.elite_fraction))
            elites = [c.sequence for c in ranked[:n_elite]]

            next_population = list(elites)
            while len(next_population) < self.population_size:
                parents = self._select_parents(ranked)
                child = self.crossover(parents) if len(parents) >= 2 else parents[0]
                child = self.mutate(child)
                next_population.append(child)

            self.population = next_population[: self.population_size]
            ranked = self._rank_population(self.population)

            mean_fitness = sum(c.fitness for c in ranked) / len(ranked)
            self.history.append(
                GenerationRecord(
                    generation=gen,
                    best_fitness=ranked[0].fitness,
                    mean_fitness=mean_fitness,
                    best_sequence=ranked[0].sequence,
                )
            )

        final_top_n = min(final_top_n, len(ranked))
        n_finalists = max(final_top_n, max(1, int(len(ranked) * self.elite_fraction)))
        finalists = ranked[:n_finalists]
        for c in finalists:
            if c.structure is None:
                c.structure = self._evaluate_structure(c.sequence)
                c.fitness = compute_fitness(c.sequence, c.structure, c.codon_usage, weights=self.weights)
        finalists.sort(key=lambda c: c.fitness, reverse=True)
        return finalists[:final_top_n]
