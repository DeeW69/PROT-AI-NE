"""Command-line entry point tying all four layers together end to end."""

from __future__ import annotations

import argparse
import sys
import uuid

DISCLAIMER = (
    "=" * 78 + "\n"
    "PROT-AI-NE explore un sous-espace guide de sequences mRNA candidates ;\n"
    "il ne garantit PAS l'optimalite globale et ne fait AUCUNE recherche\n"
    "exhaustive (l'espace theorique des mRNA possibles pour une proteine\n"
    "de taille n est de l'ordre de 6^n, non enumerable).\n"
    "\n"
    "Les sequences produites sont des candidats computationnels NON VALIDES\n"
    "EXPERIMENTALEMENT : aucune garantie de securite, d'efficacite ou\n"
    "d'expression reelle en cellule. Toute application medicale requiert\n"
    "une validation en laboratoire puis des essais cliniques.\n" + "=" * 78
)


def _read_fasta_protein(path: str) -> tuple[str, str]:
    from Bio import SeqIO

    record = next(SeqIO.parse(path, "fasta"))
    return record.id, str(record.seq)


def run_command(args: argparse.Namespace) -> None:
    from .constraints import optimize_constraints
    from .export import export_top_n
    from .genetic_engine import GeneticOptimizer
    from .prediction_3d import predict_top_n_3d
    from .reverse_translate import reverse_translate
    from .storage import get_connection, save_generation
    from .verification import check_off_target_homology, compute_cai, verify_translation

    print(DISCLAIMER)

    record_id, protein_seq = _read_fasta_protein(args.protein_fasta)
    print(f"\n[1/5] Proteine '{record_id}' chargee ({len(protein_seq)} acides amines).")

    seed = reverse_translate(protein_seq, organism=args.organism, strategy="most_frequent")
    print(f"[1/5] Retro-traduction de reference generee ({len(seed)} nt).")

    print("[2/5] Optimisation sous contraintes (DNAChisel)...")
    try:
        constrained = optimize_constraints(seed, gc_range=(args.gc_min, args.gc_max))
        seed_sequences = [constrained.sequence.replace("GCCACC", "", 1) if constrained.sequence.startswith("GCCACC") else constrained.sequence]
        print(f"       GC% = {constrained.gc_content:.1%}, contraintes satisfaites = {constrained.constraints_ok}")
    except ImportError as exc:
        print(f"       DNAChisel indisponible ({exc}); poursuite sans optimisation de contraintes.")
        seed_sequences = [seed]

    print(f"[3/5] Recherche genetique : {args.generations} generations, population {args.population}...")
    weights = {"structure": args.structure_weight, "codon_usage": args.codon_weight}
    print(f"       Ponderation fitness : structure={weights['structure']:.2f}, codon_usage={weights['codon_usage']:.2f}")
    optimizer = GeneticOptimizer(
        protein_seq=protein_seq,
        organism=args.organism,
        population_size=args.population,
        seed_sequences=seed_sequences,
        weights=weights,
        rng_seed=args.seed,
    )
    try:
        finalists = optimizer.run(generations=args.generations, final_top_n=args.top_n)
    except ImportError as exc:
        print(f"\nErreur: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"[4/5] Recherche terminee. Meilleurs candidats (top {len(finalists)}) :\n")
    for rank, candidate in enumerate(finalists, start=1):
        print(
            f"  #{rank}  fitness={candidate.fitness:.4f}  "
            f"structure={candidate.structure:.4f}  codon_usage={candidate.codon_usage:.4f}"
        )

    print(f"\n[5/5] Verification post-selection (retro-traduction de controle + CAI)...")
    cai_scores: list[float] = []
    for rank, candidate in enumerate(finalists, start=1):
        if not verify_translation(candidate.sequence, protein_seq):
            print(
                f"\nErreur: le candidat #{rank} ne se retro-traduit pas en la proteine "
                "cible (frameshift, stop premature ou substitution) ; run avorte.",
                file=sys.stderr,
            )
            sys.exit(1)
        cai_scores.append(compute_cai(candidate.sequence, organism=args.organism))
    print(f"       {len(finalists)}/{len(finalists)} candidats verifies (retro-traduction exacte).")
    for rank, cai in enumerate(cai_scores, start=1):
        print(f"  #{rank}  CAI={cai:.4f}")

    if args.check_off_target:
        print("\n       Recherche d'homologie off-target (BLAST distant NCBI, peut prendre plusieurs minutes)...")
        for rank, candidate in enumerate(finalists, start=1):
            try:
                hits = check_off_target_homology(candidate.sequence, organism=args.organism)
            except Exception as exc:
                print(f"  #{rank}  verification off-target indisponible ({exc})")
                continue
            if hits:
                print(f"  #{rank}  {len(hits)} hit(s) potentiel(s) :")
                for hit in hits:
                    print(f"       {hit.accession}  identity={hit.identity_pct:.1f}%  E={hit.e_value:.2g}  {hit.description}")
            else:
                print(f"  #{rank}  aucun hit significatif.")

    run_id = str(uuid.uuid4())
    con = get_connection(args.db_path)
    save_generation(
        con,
        run_id=run_id,
        generation=args.generations,
        candidates=[
            {
                "sequence": c.sequence,
                "codon_usage": c.codon_usage,
                "structure": c.structure,
                "fitness": c.fitness,
            }
            for c in finalists
        ],
    )
    print(f"\nRun enregistre dans {args.db_path} (run_id={run_id}).")

    if args.output:
        result = export_top_n(finalists, record_id, args.output, cai_scores=cai_scores)
        print(f"Candidats exportes en FASTA : {result.fasta_path}")
        print(f"Metadonnees structurales exportees en JSON : {result.structures_path}")

        if args.predict_3d:
            structures3d_dir = result.fasta_path.with_name(result.fasta_path.stem + "_structures3d")
            print(
                f"\nPrediction de structure tertiaire (RhoFold+, top {len(finalists)})... "
                "indicatif, non valide experimentalement, cf. avertissement viewer."
            )
            pdb_paths = predict_top_n_3d(finalists, structures3d_dir, timeout=args.predict_3d_timeout)
            n_ok = sum(1 for p in pdb_paths.values() if p is not None)
            print(f"       {n_ok}/{len(finalists)} structures 3D generees dans {structures3d_dir}")
            for rank, path in pdb_paths.items():
                status = str(path) if path is not None else "non disponible"
                print(f"  #{rank}  {status}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="protaine", description="Generation et scoring de candidats mRNA a partir d'une proteine cible.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Lancer un pipeline complet de bout en bout.")
    run_parser.add_argument("--protein-fasta", required=True, help="Fichier FASTA contenant la proteine cible.")
    run_parser.add_argument("--organism", default="human", help="Table d'usage des codons a utiliser (defaut: human).")
    run_parser.add_argument("--generations", type=int, default=50)
    run_parser.add_argument("--population", type=int, default=200)
    run_parser.add_argument("--top-n", type=int, default=10, help="Nombre de candidats finaux a rapporter.")
    run_parser.add_argument("--gc-min", type=float, default=0.4)
    run_parser.add_argument("--gc-max", type=float, default=0.6)
    run_parser.add_argument(
        "--structure-weight",
        type=float,
        default=0.5,
        help="Poids du score de structure secondaire dans la fitness composite (defaut: 0.5).",
    )
    run_parser.add_argument(
        "--codon-weight",
        type=float,
        default=0.5,
        help="Poids du score d'usage des codons dans la fitness composite (defaut: 0.5).",
    )
    run_parser.add_argument("--seed", type=int, default=None, help="Graine aleatoire pour la reproductibilite.")
    run_parser.add_argument("--db-path", default="protaine.duckdb")
    run_parser.add_argument("--output", default=None, help="Chemin de sortie FASTA pour les meilleurs candidats.")
    run_parser.add_argument(
        "--check-off-target",
        action="store_true",
        help="Verifie l'homologie off-target via BLAST distant NCBI sur les finalistes (lent, reseau, informatif seulement).",
    )
    run_parser.add_argument(
        "--predict-3d",
        action="store_true",
        help=(
            "Predit la structure tertiaire des finalistes via RhoFold+ local "
            "(necessite une installation separee, cf. README ; Linux ou WSL "
            "uniquement ; ignore le candidat en cas d'echec/timeout)."
        ),
    )
    run_parser.add_argument(
        "--predict-3d-timeout",
        type=int,
        default=300,
        help="Timeout en secondes par candidat pour la prediction 3D (defaut: 300).",
    )
    run_parser.set_defaults(func=run_command)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
