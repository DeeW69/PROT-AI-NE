"""DuckDB persistence for GA runs, generations and candidates."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import duckdb

DEFAULT_DB_PATH = "protaine.duckdb"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS generations (
    run_id VARCHAR,
    generation INTEGER,
    candidate_rank INTEGER,
    sequence VARCHAR,
    codon_usage_score DOUBLE,
    structure_score DOUBLE,
    fitness DOUBLE,
    created_at TIMESTAMP
);
"""


def get_connection(db_path: str | Path = DEFAULT_DB_PATH) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(db_path))
    con.execute(_SCHEMA)
    return con


def save_generation(
    con: duckdb.DuckDBPyConnection,
    run_id: str,
    generation: int,
    candidates: list[dict],
) -> None:
    """Persist one generation's candidates.

    Each item in ``candidates`` is expected to have ``sequence``,
    ``codon_usage``, ``structure`` (may be ``None``) and ``fitness`` keys
    (matching ``genetic_engine.Candidate``'s fields).
    """
    now = datetime.now(timezone.utc)
    rows = [
        (
            run_id,
            generation,
            rank,
            c["sequence"],
            c["codon_usage"],
            c["structure"],
            c["fitness"],
            now,
        )
        for rank, c in enumerate(candidates)
    ]
    con.executemany(
        """
        INSERT INTO generations
            (run_id, generation, candidate_rank, sequence, codon_usage_score, structure_score, fitness, created_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def get_best_candidates(con: duckdb.DuckDBPyConnection, run_id: str, n: int = 10) -> list[dict]:
    result = con.execute(
        """
        SELECT sequence, codon_usage_score, structure_score, fitness, generation
        FROM generations
        WHERE run_id = ?
        ORDER BY fitness DESC
        LIMIT ?
        """,
        [run_id, n],
    )
    columns = [c[0] for c in result.description]
    return [dict(zip(columns, row)) for row in result.fetchall()]


def get_history(con: duckdb.DuckDBPyConnection, run_id: str) -> list[dict]:
    result = con.execute(
        """
        SELECT generation, MAX(fitness) AS best_fitness, AVG(fitness) AS mean_fitness
        FROM generations
        WHERE run_id = ?
        GROUP BY generation
        ORDER BY generation
        """,
        [run_id],
    )
    columns = [c[0] for c in result.description]
    return [dict(zip(columns, row)) for row in result.fetchall()]
