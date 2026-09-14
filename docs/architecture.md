# Architecture

## Pourquoi pas une recherche exhaustive

Pour une proteine de longueur `n`, le nombre d'ARNm possibles (degenerescence
du code genetique) est de l'ordre de `6^n`. Pour une proteine de 300 acides
amines, cela depasse `10^233` combinaisons - un nombre sans rapport avec la
puissance de calcul disponible, meme cumulee sur l'age de l'univers. Il n'y a
par ailleurs aucune symetrie biologique (type "sequence inversee") permettant
de reduire cet espace : l'ARN a un sens de lecture chimique fixe (5'->3') et
est chiral, donc inverser une sequence donne un polymere different qui code
une proteine differente.

PROT-AI-NE ne tente jamais d'enumerer cet espace. Chaque couche du pipeline
reduit drastiquement l'espace avant que la suivante ne travaille dessus.

## Pipeline en 6 couches (viewer et prediction 3D restant a construire)

```text
Proteine (FASTA)
      |
      v
[1] reverse_translate.py   -- retro-traduction ponderee (deterministe, ms)
      |  une sequence de reference "raisonnable"
      v
[2] constraints.py          -- DNAChisel : GC%, sites de restriction,
      |                        homopolymeres, contexte Kozak
      |  un petit voisinage de sequences valides (dizaines-centaines)
      v
[3] genetic_engine.py       -- population, mutation (codons synonymes),
      |  <----------------+  crossover (2-4 parents), selection
      |                   |
      v                   |
[4] structure.py ----------+ ViennaRNA : structure secondaire / MFE,
      |    (uniquement sur le top X% de la population a chaque generation,
      |     pour ne pas exploser le cout de calcul)
      v
Top N finalistes seulement
      |
      v
[5] verification.py         -- retro-traduction de controle (bloquant),
      |                        score CAI, homologie off-target (BLAST,
      |                        informatif, --check-off-target)
      v
Top N candidats verifies (storage.py, DuckDB)
      |
      v
[6] export.py                -- top-N.fasta + top-N_structures.json
      |  (dot-bracket, MFE, scores detailles par candidat)
      v
      |
      v
[6b] prediction_3d.py (optionnelle, --predict-3d) -- RhoFold+ en sous-processus,
      |  un candidat a la fois, jamais dans la boucle GA ; degrade en None
      |  par candidat (RhoFold+ absent, timeout, crash) sans jamais lever
      v
viewer/ (page statique) -- selecteur top1..topN + toggle 2D/3D
```

### Couche 1 - Rétro-traduction pondérée (`reverse_translate.py`)

Un acide amine -> un codon (le plus frequent, ou un tirage pondere par
frequence) selon la table d'usage des codons de l'organisme cible
(`codon_tables.py`, donnees Kazusa). Aucune recherche : sortie unique et
quasi instantanee.

### Couche 2 - Optimisation sous contraintes (`constraints.py`)

Wrapper autour de DNAChisel (moteur de contraintes deja concu pour ce type
de probleme - non reimplemente ici). Contraintes par defaut : traduction
preservee, %GC dans une fourchette (fenetre glissante), absence
d'homopolymeres longs, absence de sites de restriction courants, objectif
d'optimisation de l'usage des codons humain. Un contexte Kozak canonique
(`GCCACC`) est prepend en amont de l'ATG (5'UTR, hors CDS).

### Couche 3 - Recherche génétique (`genetic_engine.py`)

Coeur du calcul. Population de sequences encodant toutes la meme proteine :
- **mutation** : substitution d'un codon par un synonyme aleatoire (jamais
  de changement d'acide amine, sauf mode explicite "variant d'antigene" -
  non implemente en V0) ;
- **crossover** : 2 a 4 parents, choix aleatoire du codon de l'un des
  parents a chaque position (invariant : le resultat encode toujours la
  meme proteine) - meme philosophie que le moteur genetique du bot de
  trading (population, elitisme, selection ponderee par fitness) ;
- **fitness composite** (`fitness.py`) : usage des codons (bon marche,
  calcule pour toute la population) + score de structure (couteux, calcule
  seulement sur le top `structure_eval_fraction` de la population a chaque
  generation - boucle de retroaction avec la couche 4). Ponderation
  configurable via `--structure-weight`/`--codon-weight` (0.5/0.5 par
  defaut) : utile quand le score d'usage des codons converge deja vers une
  valeur quasi uniforme sur toute la population (le seed de la couche 2
  ayant deja bien optimise), auquel cas seule la structure differencie
  encore reellement les candidats et merite plus de poids.

### Couche 4 - Validation structurale (`structure.py`)

Utilise les bindings Python de ViennaRNA (`RNA.fold`), pas l'implementation
Nussinov pure-Python de Gene Explorer (plafonnee a ~400 nt en O(n^3),
insuffisante pour une CDS de ~900 nt). Le score retourne est une heuristique
(pas un predicteur valide d'efficacite de traduction) qui favorise une
region 5' peu structuree et penalise les replis globaux tres stables.

### Couche 5 - Verification post-selection (`verification.py`)

Appliquee uniquement sur le top N final (jamais sur toute la population, a
la difference des couches 1-4) :

- **retro-traduction de controle** (`verify_translation`, bloquant) :
  retraduit chaque candidat retenu (Biopython `Seq.translate`) et compare a
  la proteine cible ; un run s'arrete en erreur plutot que d'exporter un
  resultat invalide (frameshift, stop premature, substitution) ;
- **score CAI** (`compute_cai`) : moyenne geometrique des frequences
  relatives de codons (Sharp & Li, 1987), ajoutee au rapport final et au
  FASTA exporte ;
- **homologie off-target** (`check_off_target_homology`) : BLAST distant
  NCBI (`blastn`/`nt`, via Biopython `Bio.Blast.NCBIWWW`) contre l'organisme
  cible ; informatif seulement (jamais bloquant), desactive par defaut
  (`--check-off-target`) car couteux en reseau et en temps (jusqu'a
  plusieurs minutes par candidat, file d'attente partagee NCBI).

### Couche 6 - Export (`export.py`)

Ecrit deux fichiers a partir du top N verifie : un FASTA (sequences + score
de fitness et CAI en en-tete) et un JSON compagnon (suffixe
`_structures.json` par defaut, deduit du chemin FASTA) contenant, pour
chaque candidat indexe par rang, sa structure secondaire (dot-bracket), son
energie libre (MFE) et le detail de ses scores. La structure secondaire est
recalculee ici (pas reutilisee depuis la couche 3) : peu couteux a l'echelle
du top N, et garantit une structure a jour meme pour les candidats dont
seul le score heuristique avait ete calcule pendant la recherche genetique.
Ce JSON est le contrat de donnees consomme par le futur viewer (mode 2D via
fornac) ; le mode 3D (RhoFold+, `prediction_3d.py`, non implemente) y
ajoutera un chemin vers un fichier PDB par candidat.

### Couche 6b - Prediction 3D optionnelle (`prediction_3d.py`)

RhoFold+ (Shen et al., Nature Methods 2024) n'est pas un paquet pip : c'est
un depot a cloner (<https://github.com/ml4bio/RhoFold>), avec son propre
environnement conda et un checkpoint pre-entraine a telecharger. Ce module
ne reimplemente donc pas ses classes Python internes (code de recherche,
sans garantie de stabilite) mais invoque son unique interface documentee,
le script `inference.py`, en sous-processus - un candidat a la fois, en
mode sequence seule (`--single_seq_pred`, RhoFold+ le qualifie lui-meme de
"structure de reference rapide", la version avec MSA necessitant ~900 Go
de bases de sequences).

Points de conception :

- **jamais dans la boucle genetique** (couche 3) - uniquement sur le top N
  final, via `--predict-3d` (desactive par defaut, cf. `cli.py`) ;
- **environnement separe** : `RHOFOLD_HOME` (checkout RhoFold+) et
  `RHOFOLD_PYTHON` (interpreteur de l'environnement conda dedie de
  RhoFold+, avec ses propres versions de torch/CUDA) sont resolus via
  variables d'environnement - l'environnement virtuel de PROT-AI-NE n'a
  jamais besoin de PyTorch installe ;
- **degradation gracieuse systematique** : RhoFold+ absent/mal configure,
  timeout, code de sortie non nul, fichier PDB manquant en sortie - chaque
  cas est intercepte et transforme en `None` pour ce candidat uniquement
  (jamais d'exception propagee), avec un diagnostic sur stderr. Un echec
  sur un candidat n'empeche pas les autres (meme isolation que
  `verification.check_off_target_homology`) ;
- **pas de support Windows natif** (RhoFold+ : Linux uniquement, WSL en
  repli, meme contrainte deja documentee pour ViennaRNA) ;
- fichiers ecrits en `<sortie>_structures3d/candidate_<rang>.pdb`, le
  chemin exact que `viewer/viewer.js` recherche deja pour le mode 3D.

**Etat des tests :** couvert par des tests unitaires qui mockent
`subprocess.run` (installation absente, succes, code de sortie non nul,
timeout, fichier PDB manquant, echec de lancement) et l'orchestration
top-N (isolation des echecs, nettoyage des dossiers de travail
temporaires, respect de `n`). Un run reel contre une installation RhoFold+
n'a pas ete effectue dans le cadre de ce developpement (installation
Linux/conda + telechargement de plusieurs centaines de Mo hors de portee
de l'environnement de dev Windows utilise ici) - a valider manuellement
avant un usage en conditions reelles.

### Viewer (`viewer/`)

Page statique (HTML/CSS/JS vanilla, aucune dependance CDN - `fornac` et
`3Dmol.js` sont vendorises dans `viewer/vendor/`) qui consomme directement
le JSON produit par `export.py`. Chargement via selecteur de fichier
(fonctionne en `file://`, sans serveur) ou via une URL relative si la page
est servie par un serveur HTTP.

Deux controles independants :

- **selecteur top1..topN**, construit dynamiquement a partir de
  `data.candidates` (pas de limite figee a 25) ;
- **toggle 2D/3D** : le mode 2D rend le dot-bracket via `fornac`
  (`FornaContainer`/`addRNA`, calcul local deja fait en couche 4/6, rendu
  instantane) ; le mode 3D tente de recuperer un PDB a
  `<stem>_structures3d/candidate_<rang>.pdb` (a cote du JSON charge par
  URL) via `fetch`, et affiche "structure 3D non disponible" en cas
  d'echec (fichier absent, JSON charge localement, reseau indisponible) -
  jamais d'erreur bloquante. Un chargement manuel de fichier PDB est
  propose en repli. Un bandeau permanent rappelle que la structure
  tertiaire de l'ARN reste une prediction non validee.

Les scores cles (fitness composite, %GC calcule cote client depuis la
sequence, energie libre, CAI) sont affiches a cote du rendu, communs aux
deux modes.

**Etat des tests :** le cablage DOM (chargement JSON, selecteur de rang,
toggle 2D/3D, calcul des scores, repli 3D gracieux, upload manuel de PDB)
est verifie par un harnais jsdom (hors repo, execute ponctuellement en
developpement) qui simule `fornac`/`3Dmol` pour isoler la logique du
projet ; les fichiers vendorises reels ont ete verifies separement
(chargement sans erreur, globaux `window.fornac`/`window.$3Dmol` corrects,
construction effective du graphe SVG par `fornac`). Le rendu visuel final
(mise en page du graphe apres simulation de force, rendu WebGL 3Dmol) n'a
pas ete verifie dans un navigateur reel - jsdom n'implemente pas certaines
API SVG (`baseVal`) ni WebGL necessaires a une verification bout en bout.
A confirmer manuellement en ouvrant `viewer/index.html` dans un navigateur.

## Limites connues / pistes V1+

- La couche 2 retourne une seule sequence optimisee (pas encore un
  veritable "voisinage" de dizaines de variantes contraint-valides) ; le
  voisinage exploite par la couche 3 vient surtout de l'echantillonnage
  pondere de la couche 1.
- Le contexte Kozak insere est fixe (`GCCACC`) ; il n'y a pas de veritable
  modelisation de 5'UTR.
- Filtrage immunologique (prediction d'epitopes), repliement 3D de la
  proteine, et interface graphique : explicitement hors perimetre V0/V1
  (voir le brief de demarrage du projet).
