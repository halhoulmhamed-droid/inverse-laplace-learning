# Étude 3 — Régularisation des covariances apprises pour une reconstruction robuste de données de Laplace

**Projet scientifique de Master — Halhoul Mhamed.** Université Abdelmalek Essaâdi, Faculté Polydisciplinaire de Larache. Master « Modélisation mathématique et applications et apprentissages », M2 en cours. Professeur encadrant : Pr. Chakir Tajani. Édition du 11 octobre 2026.

**Résultat principal : gain utile non établi.** La contraction vers H1 améliore ponctuellement le risque complet face à sa version γ=0 dans les cinq entraînements, mais aucune des cinq bornes supérieures corrigées n'atteint le seuil utile de 0,97. Les cinq bornes de protection en observation partielle restent sous 1,03. Le protocole exige les deux conditions : sa conclusion est donc `NOT_ESTABLISHED`. Des dégradations persistent sous rugosité et erreur de hauteur. Aucune supériorité générale ni invention n'est revendiquée.

## Lire les livrables

- [Rapport PDF, 13 pages](docs/report.pdf) et [source LaTeX](docs/report.tex).
- [Slides PDF, 14 diapositives](docs/slides.pdf) et [source LaTeX](docs/slides.tex).
- [Provenance des figures](docs/provenance.md) : liens immuables vers NPZ, script, configuration et figure au commit scientifique `867144db98eba3c0f5fe7cece4261918a46f1c31`.
- [Protocole figé](PROTOCOL.md), [mathématiques](MATHEMATICS.md), [addendum logique après campagne](MATHEMATICS_ADDENDUM.md), [sources scientifiques](SOURCES.md).
- [Audit ciblé de l'étude 2](AUDIT_STUDY2.md) et [revue indépendante de l'étude 3](INDEPENDENT_REVIEW.md).
- [Tests](results/tests.json), [pilote](results/pilot/report.json), [bilan de campagne](results/final/report.json), [analyse et décision](results/analysis/report.json), [1 200 lignes de résultats](results/analysis/summary.csv), [contrôles de livraison](results/document_checks.json), [ressources](results/resources.json).

## Ce qui est comparé

Le modèle S utilise `Cγ=(1−γ)Cemp+γ CH1` dans le direct nominal A. Cemp est la covariance du signal propre représenté ; CH1 est proportionnelle à `diag(1+(πk)²)⁻¹`, avec la même trace que Cemp. Le centre est Aμ et la covariance de bruit utilisée est `α σ²/64 I`. Positivité et solvabilité sont vérifiées, y compris à γ=0 lorsque Cemp est singulière.

M0 est la même méthode à γ=0. E est le LMMSE empirique original de l'étude 2, construit avec des covariances croisées signal–observation bruitée et le centre empirique μy ; E n'est pas confondu avec M0. Les autres comparateurs sont Tikhonov H1 historique, une covariance empirique diagonale **cosinus**, et H1m, l'extrémité γ=1 au paramètre exactement converti et au centre Aμ. Les diagnostics M0_at_S/H1m_at_S partagent l'alpha choisi pour S et isolent la modification de covariance. Ils ne remplacent pas les comparateurs réglés séparément.

Les coefficients cosinus du générateur sont indépendants. Des termes hors diagonale estimés et des corrélations dans d'autres bases ne démontrent aucune corrélation physique.

## Protocole et résultats

Gel daté avant implémentation, pilote et test au commit `8c23fccb9e00d053af71fa9349c1b761c30ca4b4`, empreintes dans [FREEZE.json](FREEZE.json). Nouvelle racine `2026101103`, partitions distinctes, cinq entraînements, tailles 32/64/256 imbriquées, validation nominale de 256 objets avec deux bruits, test de 400 objets avec quatre bruits. Même test entre entraînements ; les bruits sont appariés entre méthodes. Aucun réglage sur les vérités de test.

S examine 140 couples (7 γ × 20 α), contre 20 paramètres par comparateur. Ce budget supérieur est déclaré. Les seuils utiles de 3 % et de non-dégradation de 3 % sont conventionnels et préspécifiés. L'inférence principale rééchantillonne les objets physiques : 4 000 bootstraps par entraînement, quantile supérieur 99,5 %, dix inégalités sous Bonferroni. La couverture est approximative et conditionnelle aux entraînements observés ; la dispersion entre cinq entraînements ne fournit pas une confiance populationnelle.

| Jeu | Ratio principal S/M0 | Borne supérieure | Ratio partiel S/M0 | Borne supérieure |
|---:|---:|---:|---:|---:|
| 0 | 0,960694 | 0,984436 | 1,000000 | 1,000000 |
| 1 | 0,975915 | 1,000097 | 1,006009 | 1,029130 |
| 2 | 0,970663 | 0,990274 | 1,000431 | 1,015771 |
| 3 | 0,982300 | 1,009720 | 1,003675 | 1,007362 |
| 4 | 0,970879 | 0,985807 | 1,012636 | 1,026748 |

Les comparaisons secondaires couvrent bruit, masques, rugosité et hauteur +10 %. À N=64, les rapports des risques moyens S/M0 sont 0,9720 au nominal complet, 1,0137 en rugosité complète, 1,0266 sous hauteur +10 % complète, 1,0189 sous hauteur +10 % partielle. Diag bat S dans les cinq entraînements nominaux complets et partiels. Les résultats secondaires restent descriptifs ; ils n'ont pas servi à redéfinir le contraste principal.

## Vérifications et préservation

17 tests ciblés réussis ; pilote indépendant réussi ; 150 conditions, 90 ajustements, 1 920 000 reconstructions (dont 480 000 diagnostics), aucun candidat invalide ni perte manquante. La revue indépendante recompose les empreintes, 4 800 statistiques, 40 000 valeurs bootstrap, 540 choix de paramètres et 576 pertes indépendantes. L'audit antérieur utilise les archives de l'étude 2, sans relancer sa campagne.

La revue finale a corrigé une portée logique : trois identités proposées pour égaler E et M0 sont **suffisantes, pas nécessaires**. L'addendum est explicitement postérieur à la campagne. Le gel original, les méthodes, les données, les seuils et la décision sont inchangés.

Le rapport et les slides ont été compilés, leurs pages inspectées et leurs liens vers les données vérifiés. Logos locaux inchangés et proportions conservées. Les fichiers antérieurs au commit de base `1f11fea35d5367a355901c3866d21f5b3b21a5b0` sont contrôlés intégralement ; aucune modification hors `studies/study3/`. Les conclusions négatives antérieures, les sources LaTeX, les PDF, le portfolio et les autres dépôts sont préservés.

La campagne seule consomme 7,69 s CPU, un pic de 98,1 Mo et 63,5 Mo de sorties. Le bilan [resources.json](results/resources.json) additionne audit, tests, pilote, analyse, figures, compilation et contrôles, avec une réserve conservatrice. Limites respectées : 4 h CPU, 2 Gio RAM, 200 Mo de nouvelles sorties. Aucun GPU, service payant, FreeFem, nouvelle installation ni jeu expérimental réel. Aucun workflow CI n'est ajouté ; les contrôles GitHub présents sont vérifiés sur le commit final après push.

## Reproduire sans écraser les archives

Python 3.12 et les versions épinglées dans le `requirements.txt` racine ont été utilisés. Depuis la racine du dépôt, un audit des archives peut être exécuté vers un fichier neuf :

```bash
python studies/study3/src/bounded.py python -m studies.study3.src.independent_audit --output /tmp/study3-independent-review.json
```

Le script refuse d'écraser cette destination. Le lanceur limite les ressources et garde les journaux hors du dépôt. Pour recompiler les documents existants, avec pdflatex et les paquets LaTeX installés :

```bash
python studies/study3/src/bounded.py python -m studies.study3.src.build_documents
```

Une reproduction scientifique complète doit se faire dans une **copie séparée** où les archives sont d'abord conservées ailleurs : tests, pilote, campagne, analyse puis figures. Les commandes originales sont les modules `studies.study3.src.test_runner`, `run --pilot`, `run`, `analysis`, `figures`, lancés dans cet ordre via `bounded.py`. `run` refuse tout écrasement d'une campagne existante. Les seules dépendances de production réutilisées sont les composants hachés du direct et de l'acquisition de l'étude 2 ; leur comportement est inchangé. Les sorties peuvent différer au dernier bit selon BLAS/versions ; les versions de calcul et les empreintes de sources sont archivées.

## Attribution et contact

Code, analyse et rédaction avec assistance d'IA ; validation académique relevant de l'auteur et de l'encadrant. Les références et leur portée de lecture sont conservées. Aucun document privé tiers ou journal opérationnel n'est publié. Aucune nouvelle licence n'est attribuée ; les logos institutionnels ne valent pas approbation des résultats et conservent leur provenance dans `assets/logos/SOURCES.md`.

Halhoul Mhamed

halhoul.mhamed@etu.uae.ac.ma

halhoulmhamed@gmail.com
