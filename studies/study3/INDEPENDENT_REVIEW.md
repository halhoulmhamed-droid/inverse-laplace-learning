# Revue scientifique indépendante — étude 3

11 octobre 2026. Revue après campagne, à partir des archives figées ; aucun réentraînement ni nouvelle campagne. Le [script autonome](src/independent_audit.py) n'importe ni modèle, ni générateur, ni analyse scientifiques de production. Les [contrôles et valeurs recomposées](results/independent_review.json) permettent d'examiner le verdict.

## Conclusion contrôlée

**L'amélioration utile de 3 % n'est pas établie selon la règle figée.** Sur le contraste principal S/M0, les ratios ponctuels de risque sont tous inférieurs à 1, mais les cinq bornes supérieures dépassent 0,97. La protection de 3 % en observations partielles satisfait sa règle dans chacun des cinq entraînements. La conjonction requise échoue donc, conformément à `NOT_ESTABLISHED` dans l'analyse originale. Ce résultat n'établit pas une absence absolue d'effet.

| Entraînement | S/M0 complet | Borne supérieure | S/M0 partiel | Borne supérieure |
|---:|---:|---:|---:|---:|
| 0 | 0,960694 | 0,984436 | 1,000000 | 1,000000 |
| 1 | 0,975915 | 1,000097 | 1,006009 | 1,029130 |
| 2 | 0,970663 | 0,990274 | 1,000431 | 1,015771 |
| 3 | 0,982300 | 1,009720 | 1,003675 | 1,007362 |
| 4 | 0,970879 | 0,985807 | 1,012636 | 1,026748 |

Les bornes sont les quantiles percentile 99,5 % du bootstrap **par objet**, avec indices communs aux deux masques et aux deux méthodes, quatre bruits conservés par objet. La correction couvre les dix inégalités selon la règle préétablie, sous l'approximation du bootstrap. L'inférence est conditionnelle aux entraînements et validations observés. Les cinq entraînements partagent les objets de test ; aucune précision populationnelle sur l'aléa d'entraînement n'est revendiquée.

## Contrôles numériques et mathématiques

- Quatre empreintes du gel, huit de sources et 150 d'archives brutes concordantes ; aucune panne enregistrée ni valeur de perte non finie.
- Les 1 200 lignes de résumé, soit 4 800 valeurs (risque, modes k≥8, résidu nominal, q95 par objet), sont recomposées exactement depuis les pertes.
- Les 40 000 valeurs bootstrap sont reproduites avec un écart maximal de `2,23e-16`, par un code séparé et une taille de blocs différente.
- Les 90 ajustements archivés conservent 140 candidats pour S et 20 pour chacune des cinq autres familles. Les 540 choix de méthodes respectent scores de validation et règle d'ex æquo. L'audit contrôle ces scores archivés ; il ne réexécute pas l'ensemble des validations.
- Les deux diagnostics utilisent exactement l'alpha sélectionné pour S. Les six cas où S choisit gamma=0 sont identiques bit à bit à M0_at_S, y compris le centre. Les formules des opérateurs plug-in ont un écart maximal normalisé de `1,28e-9` par rapport aux matrices sauvegardées.
- L'extrémité gamma=1 est vérifiée directement dans l'espace des coefficients : `R = (AᵀA + alpha*nu/s Q)⁻¹ Aᵀ`, au centre `A mu`. Écart absolu maximal de `2,79e-9`. Cela ne rend pas identique le H1 historique, qui conserve sa grille et le centre empirique des observations.
- Les deux opérateurs E de N=64, jeu 0, complet et partiel, sont recomposés depuis les graines d'apprentissage et les covariances croisées empiriques bruitées. Écart maximal de `3,07e-11` ; aucun remplacement par la formule plug-in nominale.
- Une intégration analytique indépendante par sinus aux extrémités des intervalles et les nouvelles graines reproduisent 576 pertes scalaires : deux objets, deux bruits, huit opérateurs, trois métriques, six scénarios nominal/rugueux/hauteur. Écart maximal de `1,25e-15`.

### Clarification logique après campagne

La relecture documentaire a identifié une formulation trop forte dans le texte mathématique figé : les égalités `mu_y=A mu`, `Ccy=Cemp Aᵀ` et `Cyy+lambda*t_y I=A Cemp Aᵀ+beta I` sont des **conditions suffisantes** d'égalité entre E et M0, et non des conditions nécessaires en général. Des produits matriciels peuvent coïncider malgré des blocs différents ; l'écart de centre peut aussi être annulé par l'opérateur de reconstruction. Le [complément mathématique daté](MATHEMATICS_ADDENDUM.md) consigne cette correction de portée sans réécrire le gel. Les formules implémentées, paramètres, archives, seuils et conclusions ne changent pas. Les autres distinctions de familles et la portée conditionnelle des conclusions dans le rapport et les slides sont cohérentes avec le présent audit.

## Résultats secondaires à conserver

Les valeurs ci-dessous sont les **médianes de cinq ratios de risques conditionnels** à N=64, sans nouvel intervalle confirmatoire. Un ratio supérieur à 1 correspond à une dégradation ponctuelle de S face au comparateur indiqué. Toutes les tailles et tous les scénarios restent disponibles dans [summary.csv](results/analysis/summary.csv).

| Scénario | S/M0 | S/E empirique | S/H1 historique | S/Diag cosinus |
|---|---:|---:|---:|---:|
| Complet nominal | 0,970879 | 0,767400 | 0,995176 | 1,013949 |
| Partiel nominal | 1,003675 | 0,885513 | 0,665912 | 1,023538 |
| Rugueux complet | 1,011215 | 0,794543 | 1,014162 | 1,027969 |
| Rugueux partiel | 0,962663 | 0,833230 | 1,085918 | 0,985941 |
| Hauteur +10 %, complet | 1,024008 | 0,901667 | 0,997154 | 1,039298 |
| Hauteur +10 %, partiel | 1,022924 | 0,905791 | 0,630165 | 1,029992 |

Le transfert vers une hauteur différente dégrade S face à M0 dans les cinq entraînements complets (ratios 1,01736 à 1,03841), et dans quatre des cinq partiels (le cinquième vaut 1). En rugosité complète, quatre des cinq ratios S/M0 et S/H1 dépassent 1 ; en rugosité partielle, quatre ratios S/H1 dépassent 1. Les gains de S face à E dans ces archives ne constituent donc pas une supériorité générale. Diag cosinus bat également S dans les cinq entraînements nominaux complets et partiels : une covariance complète régularisée n'est pas systématiquement préférable lorsque les coefficients générateurs sont indépendants.

Dans le scénario principal, S/M0_at_S a pour médiane 0,609818, alors que S/M0 avec réglages séparés vaut 0,970879. Cet écart montre le rôle de l'alpha choisi : le diagnostic à alpha commun isole un changement de covariance, tandis que le contraste principal compare deux procédures réglées. S bénéficie de 140 candidats contre 20 pour M0. Le diagnostic ne remplace pas le contraste figé et n'autorise pas à amplifier sa conclusion.

## Reproduction et limites

Depuis la racine du dépôt, avec les dépendances du projet disponibles, le contrôle peut être répété vers une nouvelle destination :

```bash
python studies/study3/src/bounded.py python -m studies.study3.src.independent_audit --output /tmp/study3-independent-review.json
```

Le script refuse d'écraser son fichier de sortie. Le lanceur conserve mesures et journaux opérationnels hors du dépôt. Les contrôles portent sur ces archives et ces formules finies, sans nouveau théorème de robustesse, sans garantie de couverture finie et sans extrapolation à toute loi de signaux, cible H1, géométrie ou quantité de données. Les modes cosinus restent indépendants dans le générateur ; aucune corrélation physique n'est déduite des matrices apprises.
