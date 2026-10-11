# Audit mathématique — étude 2

Date : 11 octobre 2026. Audit écrit avant implémentation et évaluation finale.

## Modèle et sens des covariances
Sur (0,1) × (0,H), trace supérieure q, dérivée normale nulle en bas et sur les côtés. La base orthonormale est φ₀=1, φₖ=√2 cos(kπx). Le transfert vers le bas multiplie par sech(kπH). Les capteurs sont les moyennes exactes sur 64 intervalles. Le générateur a 128 coefficients gaussiens **indépendants**, variance proportionnelle à (1+k²)^(-p), somme des variances égale à 1. L'inverse représente 64 coefficients, avec le solveur de différences finies de l'étude 1, inchangé. Toutes les observations sont divisées par 8 : les poids restent 1/64 après masquage.

Si D est la covariance diagonale dans la base cosinus et V une base singulière, Cov(Vᵀc)=VᵀDV, généralement non diagonale pour des capteurs partiels. Ce changement de coordonnées, le modèle direct et l'erreur d'échantillonnage ne prouvent aucune corrélation physique des coefficients cosinus. Le champ spatial peut naturellement avoir une covariance non locale. On ne confond pas ces deux affirmations.

## Ablation contrôlée
Soit A=UΣVᵀ l'opérateur inverse numérique masqué, Uᵣ,Vᵣ ses 12 premières colonnes (r=12), P=VᵣVᵣᵀ. Fixer les signes des vecteurs n'est pas nécessaire : les opérateurs et pertes sont invariants au changement simultané de signe. Les éventuelles valeurs singulières égales rendent la diagonalité dépendante du choix de base ; la SVD NumPy utilisée est enregistrée, sans invariance par rotation revendiquée.

Même jeu de N couples (c,y), mêmes moyennes μc, μy, covariances empiriques centrées avec diviseur N−1. Poser z=Uᵣᵀ(y−μy), b=Vᵣᵀ(c−μc), t=tr(Czz)/r>0, λ>0.

- **D** : R_D=Vᵣ diag(diag(Cbz)/(diag(Czz)+λt)) Uᵣᵀ.
- **F** : R_F=Vᵣ Cbz(Czz+λtI)⁻¹ Uᵣᵀ.
- **O** : R_O=Ccz(Czz+λtI)⁻¹ Uᵣᵀ.
- **L** : R_L=Ccy(Cyy+λ t_y I)⁻¹, t_y=tr(Cyy)/m.
- **LP** : R_LP=P R_L.

Dans tous les cas ĉ=μc+R(y−μy), y compris pour les versions projetées : on projette l'incrément, pas la moyenne. D et F ont exactement les mêmes entrées Uᵣ et sorties Vᵣ, même régularisation et 20 candidats chacun. D remplace les blocs Cbz,Czz par leurs diagonales. Cette opération correspond à une covariance jointe des couples (b_j,z_j) rendus indépendants entre indices : chaque bloc 2×2 reste semi-défini positif. Elle n'est PAS la seule diagonalisation de Ccc dans la base cosinus.

La minimisation de E||b−Gz||²+λt||G||²_F donne G(Czz+λtI)=Cbz, d'où F. La contrainte G diagonale donne D. Avec λ>0, la Hessienne est définie positive et la solution unique. Ces optima concernent le risque empirique pénalisé, pas une garantie hors échantillon. « LMMSE » désigne ici une estimation empirique régularisée de l'opérateur affine optimal, pas le vrai estimateur oracle ni une moyenne conditionnelle exacte hors modèle gaussien connu.

**Identités d'audit :** P R_O=R_F, P R_L=R_LP, R_F UᵣUᵣᵀ=R_F. La première vaut à λ identique : une sélection indépendante des λ détruirait la comparaison exactement emboîtée. Les diagnostics fixent donc λ à celui sélectionné pour F (D/F/O) et pour L (L/LP). Les courbes principales de méthodes utilisent leur propre sélection, avec budgets identiques. Le contraste O/F à λ commun isole la sortie hors P ; L/LP isole aussi cette sortie en conservant toutes les observations. L/O change l'espace d'entrée et ne mesure pas à lui seul les corrélations.

## Décomposition exacte de l'erreur
Pour q vrai à 128 coefficients et ĉ à 64, e=ĉ−c₀:₆₄,
||ĉ−q||²_L² = ||Pe||² + ||(I−P)e||² + Σ(k=64..127)c_k².
La dernière somme est incluse dans la composante complémentaire enregistrée. P est orthogonal par SVD. À λ commun, F et O ont la même erreur projetée ; leur différence de perte totale provient uniquement du complément. Même identité L/LP. Un gain dans le complément n'est pas une récupération déterministe garantie d'un noyau : le complément spectral retenu n'est pas nécessairement le noyau exact de A.

## Références et limites
Les formules ci-dessus sont dérivées directement par différentiation quadratique et orthogonalité. Les filtres spectraux et l'apprentissage de régularisateurs quadratiques sont des méthodes établies ; voir SOURCES.md. Aucun théorème de convergence nouveau n'est proposé. L'audit sera complété par tests de gradients, projections, diagonalité, bruit apparié et décomposition indépendante avant le pilote.
