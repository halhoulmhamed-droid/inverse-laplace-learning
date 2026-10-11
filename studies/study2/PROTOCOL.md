# Protocole figé — étude 2

Date : 11 octobre 2026. Version 1, fixée avant code scientifique nouveau, pilote et test final. L'horodatage UTC exact et les SHA-256 figurent dans FREEZE.json, puis dans un commit de protocole antérieur à la campagne. Toute correction ultérieure doit être déclarée, sans réécrire ce protocole.

Titre : **Complétion de données de Laplace : apport des covariances apprises et robustesse aux observations partielles**.

## Question et contraste principal unique
L'étude 1 motive une hypothèse, pas une confirmation. Sur le masque contigu de 32 capteurs sur 64, σ=0,01, H=0,25, p=2,5, N=256, r=12 : F (covariances complètes dans le sous-espace) améliore-t-il D (diagonales dans exactement les mêmes coordonnées) ? Les méthodes sont définies et auditées dans MATHEMATICS.md.

Perte inverse : erreur **quadratique L² absolue** sur la trace supérieure complète à 128 modes, par objet, moyennée sur quatre bruits. Risque R : moyenne de ces pertes objets. Le signal a E||q||²=1. Seuil utile : R_F/R_D ≤ 0,95, soit au moins 5 % de risque quadratique, environ 2,53 % en racine du risque. Ce seuil de décision conventionnel, fixé avant les nouvelles données, exige une amélioration modeste mais non négligeable en contrepartie d'un bloc 12×12 au lieu de 12 gains ; ce n'est ni un seuil clinique ni une valeur dérivée du test ou de l'ancien seuil de 10 %.

Protection contre les queues : q95 des pertes objets de F / q95 de D ≤ 1,05 (marge maximale 5 %). Ce second critère porte sur la même perte inverse, sans multiplier les contrastes principaux.

Cinq entraînements indépendants. Pour chaque entraînement, bootstrap apparié par **objet physique**, conservant les quatre répétitions et les deux méthodes, 4 000 tirages. Bornes supérieures percentile au niveau 99,5 % pour les deux ratios : Bonferroni 0,05/(5×2). Une conclusion confirmatoire positive exige les deux bornes sous leurs seuils dans chacun des cinq entraînements. Sinon : « non établi », même si les points sont favorables. Couverture nominale approximative du bootstrap, sans garantie finie ; pas de mélange des objets et entraînements en pseudo-réplications. Pas d'intervalle populationnel sur l'aléa d'entraînement avec seulement cinq tirages. Résumer séparément leur dispersion. Indiquer chaque ratio ponctuel >1 comme dégradation et les bornes dépassant les seuils comme incertitude/inadéquation au critère. En cas de dénominateur nul ou bootstrap non fini : résultat indécidable, sans suppression sélective.

## Données nouvelles, partitions et dimensions
Graine racine 2026101102, différente de l'étude 1. PCG64/SeedSequence [racine, partition, entraînement, objet, canal, répétition]. Partitions 10 apprentissage, 20 validation, 30 test nominal/hauteur, 31 test rugueux ; ajouter 100 à chaque partition pour le pilote. Canal 0 coefficients, 1 bruit. Le bootstrap utilise partition 40. Objets 0..N−1, sans chevauchement entre partitions. Entraînements g=0..4, tailles 64/256/1024 imbriquées dans chaque g ; les g sont indépendants. Validation : 256 nouveaux objets par g, deux bruits chacun, partagés entre tailles et méthodes. Test : 400 objets communs aux cinq entraînements (g=0 dans la graine de test), quatre bruits chacun. Nominal et hauteur utilisent les mêmes objets, rugosité une partition indépendante. Vecteurs de bruit de 64 capteurs identiques entre méthodes, masques et niveaux pour un objet/répétition, puis masqués et multipliés par σ/8.

Générateur analytique 128 modes, coefficients cosinus indépendants, p=2,5 nominal et p=1,5 rugueux. Inverse FD : 64 coefficients, nx=64, ny=16, H=0,25, 64 moyennes de capteurs, puis masque et division par 8. Réutiliser direct.py inchangé. Aucune vraie corrélation physique ajoutée. Entraînement et validation toujours nominal p=2,5/H=0,25 ; chacun est bruité une fois en apprentissage et deux fois en validation.

Dix scénarios, pas de produit cartésien supplémentaire :
1. principal : 32 contigus, σ=0,01, nominal ;
2. 64 complets ; 3. 32 régulièrement espacés ; 4. 16 contigus ; 5. 16 régulièrement espacés ;
6. 32 contigus, σ=0,003 ; 7. 32 contigus, σ=0,03 ;
8. 32 contigus, σ=0,01, rugueux ;
9. 32 contigus, σ=0,01, H vrai=0,2625 (+5 %) ;
10. idem H vrai=0,275 (+10 %).
Entraîner au masque et niveau de bruit connus, sans réentraînement pour rugosité/hauteur. Ces tests sont secondaires descriptifs, sans revendication confirmatoire multiple.

## Méthodes et réglage équitable
Huit méthodes : D, F, O, L, LP, Tikhonov H1, filtre signé plafonné et filtre positif borné/plafonné, ces deux derniers sur les mêmes 12 directions. Tous partagent moyennes empiriques et observations d'apprentissage.
20 λ=logspace(−8,1,20) pour D/F/O/L/LP avec échelle t adaptée à l'entrée ; 20 α=logspace(−10,1,20) pour H1, pénalité diag(1+(kπ)²) normalisée par son maximum ; 20 plafonds logspace(0,4,20) pour les filtres spectraux, gains cov(b_j,z_j)/var(z_j), puis clip signé ou [0,min(plafond,1/s_j)].
Sélection **supervisée sur validation uniquement** : minimiser la perte quadratique moyenne sur 256 objets × deux bruits. Les vérités de validation sont disponibles en simulation, celles de test ne sont jamais données au sélecteur. Ceci diffère volontairement de la validation observationnelle de l'étude 1 : aucune comparaison chiffrée directe entre études. Ex æquo à 1e−12 max(1,minimum) : norme opérateur minimale, puis indice de grille. Aucun agrandissement de grille après résultats. Signaler sélections aux extrémités.

Diagnostics à λ partagé : D/F/O au λ sélectionné pour F et L/LP au λ sélectionné pour L ; enregistrer pertes totales/projetées/complémentaires et vérifier les identités. Ces diagnostics isolent les restrictions d'espace sans attribuer à l'ablation un gain dû à un λ différent. D/F restent le contraste principal avec réglage indépendant égalitaire. Une unique base r=12 est étudiée, sans généralisation à tous les rangs.

## Tests, pilote, gel et ressources
Tests ciblés : optimum quadratique contre moindres carrés augmentés, D=F pour covariances diagonales, identités de projection à λ commun, décomposition L² avec modes omis, invariance de signe SVD, déterminisme et disjonction des partitions, sélecteur sans vérité de test, bootstrap par objet et échecs non éliminés.
Pilote séparé : 1 entraînement N=64, 24 objets validation, 24 test, deux répétitions, scénarios principal/complet/16 contigus/hauteur+10 %. Ne sert qu'à tester coût, finitude, stockage et identités. Aucun choix de contraste ou grille selon ses performances. Avant campagne, exiger tests réussis, aucune panne numérique du pilote, CPU extrapolée ×10 < 4 h, RAM <2 Gio, sorties projetées avec réserve documentaire <200 Mo. Sinon arrêter la campagne et documenter le blocage, sans examiner le test final.
Plafonds globaux de nouveaux calculs : 14 400 s CPU, 2 Gio RAM, 200 000 000 octets de sorties. BLAS mono-thread, CPU/RAM limités par le lanceur, ledger privé des commandes scientifiques et documentaires ; sorties numériques publiques, journaux privés. Préserver résultats négatifs, pannes et valeurs manquantes explicites. Aucune sélection de sous-ensemble après résultats.

## Livraison et préservation
Branche dédiée study2-covariances-20261011. Tous les nouveaux fichiers sous studies/study2/. Aucun fichier de l'étude 1 modifié, y compris ses PDF et sa conclusion. Rapport/slides distincts, sources LaTeX et PDF, figures liées aux NPZ/scripts/configuration/commit scientifique ; identité académique et logos déjà vérifiés. Pas de nouveauté ni supériorité générale revendiquées. Git normal, aucune modification de portfolio ou d'autre dépôt.
