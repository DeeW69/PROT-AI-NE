"""PROT-AI-NE: generation et scoring de candidats mRNA a partir d'une proteine cible.

Pipeline en 4 couches (voir docs/architecture.md) :
    1. reverse_translate  - retro-traduction ponderee par usage des codons
    2. constraints        - optimisation sous contraintes (DNAChisel)
    3. genetic_engine      - recherche genetique dans le voisinage contraint
    4. structure           - validation de structure secondaire (ViennaRNA)

Avertissement : les candidats produits sont computationnels et n'ont pas de
valeur experimentale ou clinique sans validation en laboratoire.
"""

__version__ = "0.1.0"
