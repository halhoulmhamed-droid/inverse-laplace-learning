# Étude 2 — Covariances apprises et observations partielles

**Projet scientifique de Master — Halhoul Mhamed.** Une étude distincte, sous branche `study2-covariances-20261011`. Les 209 fichiers de l'étude 1, dont sa conclusion négative, ses sources et PDF, sont préservés octet pour octet. La branche `main` n'est pas modifiée.

- [Rapport français — PDF, 13 pages](docs/report.pdf) et [source LaTeX](docs/report.tex)
- [Slides — PDF, 14 pages](docs/slides.pdf) et [source LaTeX](docs/slides.tex)
- [Protocole figé](PROTOCOL.md), [gel horodaté et empreintes](FREEZE.json), [audit mathématique](MATHEMATICS.md)
- [Provenance des cinq figures](docs/provenance.md), [références et limites de lecture](SOURCES.md)
- [Analyse primaire](results/analysis/report.json), [toutes les méthodes et scénarios](results/analysis/summary.csv), [vérification indépendante](results/verification.json)
- [Contrôles documentaires](results/document_checks.json) et [budget de livraison](results/resources.json)

## Conclusion et portée

Le contraste préspécifié compare F (blocs de covariances complets) à D (leurs diagonales) **dans les mêmes 12 directions d'entrée et de sortie**, à données, moyennes et budgets de réglage identiques. À N=256, avec 32 capteurs contigus et bruit absolu 0,01, les cinq rapports de risques sont entre 0,3192 et 0,3218. Les dix bornes conditionnelles satisfont les critères figés : réduction utile de risque quadratique de 5 % et marge de non-dégradation de 5 % sur le q95 des pertes objets. Bootstrap par objet apparié, 4 000 tirages, correction des deux critères sur cinq entraînements ; couverture nominale approximative, pas de confiance populationnelle sur tous les entraînements.

LMMSE complet L bénéficie aussi de la prédiction **hors du sous-espace retenu**, isolée grâce à O/F et L/LP à ridge partagé. F seul reste moins performant que H1 dans le scénario principal. L/H1 vaut 0,753 en ratio de risques moyens, mais ce contraste secondaire est descriptif. L se dégrade face à H1 sur tous les capteurs (+9,2 %) et sur les objets rugueux (+10,7 %). Avec 16 capteurs contigus, le risque moyen baisse mais le q95 moyen augmente de 4,4 %. Les résultats négatifs sont conservés.

Les coefficients cosinus du générateur restent indépendants. Les covariances non diagonales dans la base singulière peuvent être induites par le changement de base, les masques et l'échantillonnage ; elles ne prouvent aucune corrélation physique ajoutée. Une seule troncature r=12 est étudiée. Les méthodes sont établies ; aucune nouveauté ni supériorité générale n'est revendiquée.

## Jalons vérifiables

- Base préservée de l'étude 1 : `9bdee8f250e1157a99d6f7dece42ff48bc2a64d0`.
- Protocole et mathématiques avant code/pilote/test : commit `20a663d`.
- Code, 11 tests ciblés et pilote séparé réussi : commit `0337094`.
- Campagne, analyse, cinq figures et contrôle indépendant : `57ab54ebc994029be3c24a8a34da56371547f4ac`.
- 150 conditions, 1 920 000 évaluations des huit méthodes et 1 200 000 diagnostics, sans panne numérique. Les diagnostics ne sont pas des objets indépendants supplémentaires.
- 468 valeurs scalaires recomposées via intégrales de sinus et graines, erreur maximale 1,98e-15 ; les 150 NPZ de pertes sont hachés.

## Lire les résultats sans relancer la campagne

Tous les chemins ci-dessous sont relatifs à ce dossier. `results/final/configuration.json` décrit graines, scénarios, versions et axes.

- `*_loss.npz` : `loss` de forme `(400 objets, 4 bruits, 13 méthodes/diagnostics, 3 composantes)` ; dernière dimension = total, projeté, complémentaire (modes 64–127 inclus). `status=0` indique une entrée valide.
- Les huit premières méthodes : D, F, O, L, LP, H1, signé plafonné, positif borné. Les cinq suivantes : D/F/O au λ de F puis L/LP au λ de L.
- `*_models.npz` : matrices affines sélectionnées, moyennes, bases U/V, covariances et valeurs singulières. `*_selection.json` conserve les scores des 160 candidats et les paramètres sélectionnés.
- `truth.npz` conserve les objets finaux nominaux et rugueux. Les autres données se régénèrent depuis les graines sans téléchargement.
- `results/analysis/bootstrap.npz` : cinq entraînements × 4 000 réplications × deux ratios.
- `figures/*.npz` : valeurs exactes représentées ; `figures/manifest.json` et `docs/provenance.md` relient données, scripts et configuration au commit scientifique.

## Reproduction dans le cloud

Dépendances gratuites déjà utilisées : Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 ; pdfTeX/LaTeX pour les documents. Aucun GPU, FreeFem ni service payant. Les versions exactes de Python/NumPy et le SHA du solveur importé sont conservés dans la configuration.

Depuis la racine du dépôt, on peut recompiler les documents avec :

```bash
python studies/study2/src/bounded.py python -m studies.study2.src.build_documents
```

Le lanceur impose BLAS mono-thread, limites de mémoire et CPU, et écrit son journal sous `/workspace/inverse-laplace-private/study2/`, hors du dépôt. Le répertoire privé doit être accessible à l'exécuteur. Les PDF sont issus des sources et des figures conservées ; aucune expérience n'est exécutée par cette commande.

Pour reproduire toute la campagne, travailler dans **une copie isolée** du dépôt, conserver une copie des résultats publiés, puis libérer seulement les sorties générées de cette copie (`results/pilot`, `results/final`, `results/analysis`). Les commandes refusent d'écraser une campagne ou une analyse existante. Relancer les tests génère le rapport de tests et les empreintes des sources présentes, puis :

```bash
python studies/study2/src/bounded.py python -m studies.study2.tests.test_study2
python studies/study2/src/bounded.py python -m studies.study2.src.run --pilot
python studies/study2/src/bounded.py python -m studies.study2.src.run
python studies/study2/src/bounded.py python -m studies.study2.src.analysis
python studies/study2/src/bounded.py python -m studies.study2.src.verify
python studies/study2/src/bounded.py python -m studies.study2.src.figures
```

Les temps, versions d'environnement et métadonnées PDF peuvent différer ; les graines et le protocole sont fixes. Une reprise interrompue conserve les sorties partielles et exige une inspection avant de recommencer. Le protocole n'est pas modifié selon les résultats.

Code, analyse et rédaction avec assistance d'IA ; validation académique relevant de l'auteur et de l'encadrant. Logos institutionnels réutilisés avec proportions conservées, sans approbation institutionnelle implicite ; [sources et empreintes](../../assets/logos/SOURCES.md). Aucun document tiers privé ni journal opérationnel n'est publié.
