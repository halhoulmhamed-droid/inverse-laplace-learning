# Protocole indépendant figé — étude 3

11 octobre 2026. Version 1, avant nouveau code scientifique, pilote et test. Horodatage UTC, SHA-256 et base Git consignés dans FREEZE.json. Aucune positivité du résultat présumée.

Titre : **Régularisation des covariances apprises pour une reconstruction robuste de données de Laplace**.

## Question et décision principale
L'objectif est de déterminer si le mélange S vers CH1 réduit une dégradation liée à l'estimation de covariance tout en préservant un scénario partiel favorable. Le LMMSE empirique E et le plug-in M0 sont des familles distinctes ; l'étude 2 motive l'hypothèse, aucune de ses données finales n'entre dans cette étude.

Contraste principal unique : **S contre M0**, N=64, 64 capteurs complets, σ=0,01, H=0,25, p=2,5. Métrique : erreur quadratique L² absolue sur la trace complète à 128 modes, moyenne sur quatre bruits puis sur 400 objets. Gain utile : ratio R_S/R_M0≤0,97 (3 % de risque, environ 1,51 % de racine du risque). Ce seuil conventionnel modeste vise une amélioration reproductible d'une méthode ayant un paramètre de plus ; ni optimum déterminé sur test ni reprise automatique des seuils antérieurs.

Marge de non-dégradation : sur le scénario **32 capteurs contigus nominal**, mêmes N=64, σ=0,01 et population de test, R_S/R_M0≤1,03. Les modèles sont ajustés séparément pour chaque masque connu. Il s'agit de préserver le risque de M0, pas d'affirmer d'avance que M0 bat H1. Les comparaisons effectives à E et H1 seront rapportées.

Cinq entraînements indépendants. Pour chacun, 4 000 bootstraps appariés par objet physique, conservant ensemble les quatre bruits, les deux méthodes et les deux masques. Bornes supérieures percentile à 99,5 %, correction Bonferroni de dix inégalités (5 jeux × gain/non-dégradation). Conclusion positive seulement si, dans chacun des cinq jeux, la borne principale≤0,97 et la borne de protection≤1,03. Sinon « non établi », avec distinction descriptive entre ratios >1, amélioration ponctuelle sans seuil atteint, et échec numérique/indécidable. Une sélection γ=0 est un résultat admissible, sans le transformer en succès utile.

Couverture bootstrap nominale approximative ; intervalles conditionnels aux apprentissages/validations observés. Aucune inférence populationnelle sur l'aléa d'entraînement avec cinq jeux seulement ; sa dispersion est présentée séparément. Dénominateur nul, panne ou tirage non fini : indécidable, aucune élimination des entrées défavorables.

## Modèle et méthodes
Formules, centres, positivité et extrémités dans MATHEMATICS.md. Aucune troncature singulière imposée aux méthodes nouvelles : toutes prédisent 64 coefficients. Générateur analytique 128 modes indépendants, E||q||²=1 ; inverse FD nominal 64 modes, nx=64, ny=16, H=0,25, mêmes composants directs inchangés. Observations et bruit sont divisés par 8, poids 1/64 après masque. Le bruit connu vaut ν=σ²/64 I.

S : Cγ=(1−γ)Cemp+γsQ⁻¹, Qkk=1+(πk)², s=tr(Cemp)/tr(Q⁻¹). γ∈{0;0,01;0,03;0,1;0,3;0,6;1}, α∈logspace(−3,3,20), β=αν, soit 140 candidats. M0 (γ=0), Diag cosinus et H1m (γ=1) ont chacun les mêmes 20 α. E empirique original : 20 λ=logspace(−8,1,20). H1 historique : 20 a=logspace(−10,1,20), pénalité Q/Qmax. Les deux derniers sont les familles exactes de l'étude 2, centrées sur μy empirique, les autres sur Aμ. Tous partagent μ signal et couples d'apprentissage.

Sélection sur validation synthétique uniquement : minimiser l'erreur quadratique moyenne sur les deux répétitions, vérités de validation autorisées ; aucun accès du sélecteur aux vérités de test. Ex æquo ≤1e−12 max(1,min score) : norme opérateur puis indice de grille. Pas d'élargissement des grilles après test ; signaler bornes choisies. Le supplément de recherche de S est explicite. Diagnostics **M0_at_S** et **H1m_at_S** sans réglage supplémentaire, au α de S, isolent l'effet du mélange. Ils ne sont pas un second contraste confirmatoire.

## Partitions nouvelles et budgets de données
Racine 2026101103, PCG64/SeedSequence [racine, partition, jeu, objet, canal, répétition]. Partitions 10 apprentissage, 20 validation, 30 test nominal/hauteur, 31 test rugueux ; ajouter 100 pour le pilote. Canal 0 coefficients, 1 bruit ; partition 40 pour bootstrap. Les premières composantes diffèrent des études 1/2. Identifiants objet 0..n−1 dans chaque namespace. Jeux d'apprentissage g=0..4 indépendants ; tailles **32,64,256** imbriquées dans chaque g, un bruit par objet. Validation : 256 objets par g, deux bruits, nominale seulement p=2,5/H=0,25. Mêmes objets validation entre tailles et méthodes.

Test : 400 objets, quatre bruits, jeu=0 dans la graine, communs aux entraînements. Nominal et erreur de hauteur partagent les objets ; rugosité p=1,5 utilise la partition31. Bruit64 apparié avant masquage entre méthodes, masques et niveaux. Le changement de hauteur et la rugosité ne sont pas utilisés pour régler les paramètres.

Dix scénarios, sans produit cartésien complet :
1. **primary** : 64 capteurs, σ=0,01, nominal ;
2. **partial** : 32 contigus, σ=0,01, nominal (protection) ;
3. **even32** : 32 espacés, nominal ;
4. **quarter** : 16 contigus, nominal ;
5. **noise_low** : complets, σ=0,003 ;
6. **noise_high** : complets, σ=0,03 ;
7. **rough_full** : complets, p=1,5 ;
8. **rough_partial** : 32 contigus, p=1,5 ;
9. **height_full** : complets, H vrai=0,275 (+10 %) ;
10. **height_partial** : 32 contigus, H vrai=0,275.
Bruit0,01 et autres paramètres nominaux sauf indication. Apprentissage/validation au masque et bruit connus ; modèles conservés sans réentraînement pour hauteur/rugosité. Nouvelles sorties brutes : perte L², erreur des modes k≥8 (y compris omis) et résidu nominal quadratique, statuts pour six méthodes + deux diagnostics. Les q95 sont descriptifs.

## Audit préalable, tests et pilote
Avant la campagne, audit des archives de l'étude 2 : gel, partitions, bootstrap, calculs scalaires et projections à paramètres partagés. Aucun recalcul exhaustif de sa campagne. Anomalie substantielle : résolution et addendum avant toute utilisation scientifique.

Tests : PSD/traces et limite γ0 déficiente ; SPD observation et résidu solveur ; identité γ1/H1 au paramètre exact et centre correct ; contre-exemple séparant LMMSE empirique et plug-in ; diagonalisation dans la base déclarée ; pertes Parseval ; sélection/grilles diagnostics ; graines disjointes, bruits appariés ; bootstrap par objet et pannes.

Pilote : un jeu N=32, 24 validation et 24 test, deux bruits, scénarios primary/partial/rough_full/height_partial. Partitions dédiées, aucun réglage du protocole selon ses performances. Gate : tests validés, aucune panne numérique du pilote, extrapolation CPU×10<4h, RSS<2Gio, stockage final projeté+réserve<200Mo. Sinon blocage scientifique à résoudre avant le test, sans diminuer opportunément le contraste.

## Ressources, gel, analyse indépendante et livraison
Maximum global :14 400 sCPU,2Gio RAM,200 000 000octets de nouvelles sorties ; BLASmono-thread, limitesprocessus, suivi privé incluant audits, tests, pilote, campagne, analyse et documents. Les réserves conservatrices couvrent les contrôles légers non mesurés individuellement. Aucun service payant/GPU/FreeFem.

Tous échecs, statuts et limites conservés. La campagne est figée avant évaluation ; une revue indépendante contrôle l'interprétation et les chiffres avant publication. Branch study3-covariance-h1-20261011, nouveaux fichiers uniquement studies/study3/. Préservation totale des études1/2 et de main. Rapport/slides propres, sourcesLaTeX/PDF, logos et identité existants, liens figures vers données/script/configuration/commit scientifique. Journaux opérationnels privés. Publication normale de la branche seulement, aucune revendication d'invention ou de résultat non établi.
