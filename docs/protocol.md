# Projet scientifique de Master — Protocole scientifique

Halhoul Mhamed — 11 octobre 2026.

Édition scientifique du plan fixé avant les expériences. Les contrôles exécutés et leurs résultats figurent dans les rapports de calcul ; les formulations prospectives ci-dessous conservent les choix préalables. Les sources sont décrites dans `sources.md`.

### 2. Modèle continu et preuve du transfert

Soit Ω_H=(0,1)×(0,H), H=1/4. On cherche u harmonique avec u(x,H)=q(x), ∂_y u(x,0)=0 et ∂_x u(0,y)=∂_x u(1,y)=0. Le flux normal inférieur est donc −∂_y u=0. On prend φ₀=1 et φₖ(x)=√2 cos(kπx), k≥1, base orthonormée de L2(0,1). Pour q de développement cₖ, la solution formelle est
\[
u(x,y)=c_0+\sum_{k\ge1} c_k\phi_k(x)
\frac{\cosh(k\pi y)}{\cosh(k\pi H)}.
\]
Pour une somme finie, différentiation directe : les dérivées secondes en x et y s'annulent, les dérivées normales imposées sont nulles et la trace haute vaut q. L'unicité du problème direct suit de l'identité d'énergie pour la différence de deux solutions : intégrale de |∇w|² nulle, puis w constant et sa trace supérieure nulle. Pour q∈H^{1/2}, la solution faible et ses traces sont comprises au sens usuel ; pour q∈L2, le transfert inférieur se définit par prolongement borné de la formule modale, sans prétendre que toute donnée L2 engendre une solution H1 jusqu'au bord supérieur.

Le transfert continu T_H:L2→L2 vérifie T_Hφₖ=tₖφₖ, t₀=1, tₖ=sech(kπH)>0. Il est injectif et compact, avec tₖ→0. Son inverse sur son image n'est pas borné : ||φₖ||=1, ||T_Hφₖ||=tₖ→0. La condition d'appartenance à l'image est la sommabilité de |dₖ/tₖ|². On n'assimile pas toute donnée bruitée à une trace compatible.

Le mode constant est conservé, estimé et régularisé. Aucun retrait de moyenne spatiale implicite n'est autorisé. Le centrage statistique décrit plus bas est différent d'une suppression du mode k=0.

### 3. Capteurs, normes et opérateurs distincts

Pour j=0,…,63, I_j=[j/64,(j+1)/64], ℓ=1/64. Le capteur retourne S_jv=ℓ^{-1}∫_{I_j}v(x)dx. Les trois masques conservent tous les j, les j pairs, ou j=0,…,31. Ils sont connus et fixes. On note M leur matrice de sélection et W=ℓI_m, m=64 ou 32. La norme des observations physiques est ||y||_W=||W^{1/2}y||₂. Les poids restent ℓ pour un demi-bord ; ils ne sont pas renormalisés à une somme unité. Cette norme mesure l'énergie de valeurs moyennes par cellules, pas la norme L2 exacte d'une trace arbitraire.

Pour k≥1,
\[
S_j\phi_k=\sqrt2\frac{\sin(k\pi b_j)-\sin(k\pi a_j)}{k\pi\ell},
\quad S_j\phi_0=1.
\]
Le générateur principal G_{H,K}=MST_HJ_K emploie K=128 coefficients (k=0,…,127), J_K leur synthèse. Il calcule ces intégrales analytiquement. La projection P_nq=J_nc utilise n=64 coefficients (k=0,…,63). Les coefficients sont euclidiens : ||J_nc||_{L2}=||c||₂. La matrice inverse A_raw provient d'un solveur par différences finies sur 65×17 points et des moyennes de son interpolant linéaire inférieur. Cette intégration numérique est explicitement distincte de la formule exacte S_jφₖ.

On pose a=||W^{1/2}A_raw||₂>0, Q=W^{1/2}/a, A=QA_raw. On définit y_raw=ST_{H*}q+ε_raw∈R⁶⁴, y_mask=M y_raw∈Rᵐ et d=Qy_mask=QM y_raw ; A_raw∈R^{m×n}, Q∈R^{m×m}. On calcule A=UΣVᵀ pour chaque masque. Les singularités σ_j sont celles de A, non les tₖ continus. Le rang algébrique r₀ et le rang numérique r ne sont pas confondus : r retient σ_j>10^{-12}σ₁ ; tous les modes rejetés sont explicitement traités comme un sous-espace tronqué, qui peut inclure des modes physiquement non nuls.

P=diag(1+(kπ)²), P̄=P/max diag(P). Tikhonov H1 pénalise vᵀP̄v dans les coefficients physiques ; la normalisation de A ne transforme pas cette base en base H1. Dans les unités brutes, α correspond à a²α/max diag(P) devant P lorsque le résidu est pondéré par W. Les plafonds L ci-dessous s'appliquent à d. Leur amplification par rapport à ||·||_W vaut L/a ; par rapport à la norme brute euclidienne, L√ℓ/a.

### 4. Décomposition d'erreur avec offset, noyau et projection

Notons c les n premiers coefficients vrais, q_⊥=q−J_nc. Pour une hauteur vraie H*, écrivons les observations normalisées
\[
d=(A+E)c+b+\varepsilon,
\]
où E=Q(MST_{H*}J_n−A_raw), b=QMST_{H*}q_⊥, ε=QMε_raw∈Rᵐ, avec ε_raw∈R⁶⁴. E peut être scindé en changement de hauteur et erreur du direct inverse ; b est un résidu de troncature transmis. Pour une vérité finie K, q_⊥ contient n,…,K−1. Si on interprète cette vérité comme approximation d'un champ infini, la queue au-delà de K doit être ajoutée ; aucune garantie pour ce champ infini n'est revendiquée par les simulations finies.

Avec moyennes entraînement μ_c, μ_d, l'estimateur affine est
\[
\widehat c=\mu_c+R(d-\mu_d),\quad
R=V_r\operatorname{diag}(g)U_r^\top.
\]
En posant f=Aμ_c−μ_d, on obtient exactement
\[
\widehat c-c=(RA-I)(c-\mu_c)+R(Ec+b+\varepsilon+f).
\]
La présence de f est indispensable quand les moyennes empiriques de données ne coïncident pas avec Aμ_c. Puisque ||R||₂=L_g=max|g_j|,
\[
\|\widehat c-c\|_2\le
\|(RA-I)(c-\mu_c)\|_2+
L_g(\|E\|_2\|c\|_2+\|b\|_2+\|\varepsilon\|_2+\|f\|_2).
\]
Pour β=Vᵀ(c−μ_c), le carré du biais vaut
Σ_{j≤r}(g_jσ_j−1)²β_j² + Σ_{j>r}β_j²,
en complétant V en base orthonormée. Le second terme contient le noyau véritable et les directions rejetées numériquement. Il n'est jamais supprimé du bilan.

L'erreur fonctionnelle satisfait l'identité
\[
\|\widehat q-q\|_{L2}^2=\|\widehat c-c\|_2^2+\|q_\perp\|_{L2}^2,
\quad \widehat q=J_n\widehat c.
\]
On peut majorer sa racine par la borne précédente plus ||q_⊥||. L’identité affine et sa borne restent valides pour toute matrice R∈R^{n×m}, notamment R_H1=(AᵀA+αP̄)^{-1}Aᵀ et R_LMMSE=Ĉ_cd(Ĉ_dd+λs_CI)^{-1}. Dans la borne, remplacer L_g par ||R||₂, plus grande valeur singulière calculée de la matrice effective. Ni max|g| ni la formule diagonale du biais ne s’appliquent à ces deux familles en général ; conserver alors ||(RA−I)(c−μ_c)||₂. Leur image n’est pas supposée incluse dans V_r.

Cette borne est conditionnelle aux valeurs effectivement bornées de E,b,ε,f. Les normes exactes du modèle analytique discret peuvent être évaluées ; une estimation de discrétisation ou un arrondi non majoré ne devient pas un certificat. Les diagnostics numériques sont présentés comme tels.

### 5. Bruit, calibration et probabilités

Pour chaque jeu d'entraînement, s_train² est la moyenne des carrés des 64 observations propres nominales sur tous ses objets. Cette unique échelle est conservée pour les trois masques et tous les tests associés. Pour ν∈{0,001;0,01;0,05}, τ=νs_train. On ajoute ε_raw∼N(0,τ²I₆₄) **après** les moyennes, puis applique le masque et Q. Il s'agit de bruit de capteur sur les moyennes, non d'un bruit blanc continu moyenné. La covariance normalisée est C_ε=(τ²ℓ/a²)I_m. Il n'y a pas de recalibration sur un profil ou une géométrie test.

Les méthodes partagent exactement les mêmes bruits. On réutilise les 64 composantes pour comparer les masques et les mêmes normales standard pour les niveaux ν : les comparaisons entre conditions sont appariées, et ne sont pas traitées comme des réplications indépendantes. Les répétitions 1 et 2 d'un objet de validation sont indépendantes conditionnellement à celui-ci ; elles ne réutilisent aucun bruit d'entraînement.

Si η∼N(0,I_m), ||ε||=(τ√ℓ/a)||η||. Une borne probabiliste élémentaire, sans concentration fine, vient de Markov : avec probabilité ≥1−ζ,
||ε||≤(τ√ℓ/a)√(m/ζ).
Elle est conservatrice et ne vaut pas simultanément pour toute une campagne sans ajustement. L'analyse principale utilise la norme réalisée et les distributions empiriques. Aucun δ fini ne borne uniformément toutes les réalisations gaussiennes presque sûrement.

### 6. Apprentissage scalaire et Wiener

Sur l'entraînement seulement, calculer μ_c et μ_d, puis β_i=V_rᵀ(c_i−μ_c) et z_i=U_rᵀ(d_i−μ_d). Pour une composante centrée,
\[
\mathcal R(g)=E[(\beta-gz)^2]
=E[\beta^2]-2gE[\beta z]+g^2E[z^2].
\]
Si Var(z)>0, g*=Cov(β,z)/Var(z). Si Var(z)=0, z=0 presque sûrement et le risque ne dépend pas de g ; on choisit g=0. Pour les variables non centrées, la classe optimale affine a intercept Eβ−gEz ; imposer un intercept nul donnerait un ratio de moments non centrés, un autre problème.

Dans le modèle nominal exact z=σβ+ε, Cov(β,ε)=0, Eε=0, p=Varβ, v=Varε, le gain vaut σp/(σ²p+v). C'est le Wiener scalaire. Avec t=Cov(β,ε), le ratio devient (σp+t)/(σ²p+2σt+v). Sous E≠0 ou mélange de modes, on conserve le ratio de covariance mais on ne lui attribue pas l'expression nominale. Ces gains peuvent être négatifs : les interdire est une ablation, pas une conséquence générale de l'optimalité.

L'estimateur empirique est \hat g_j=Σ_iβ_ij z_ij / Σ_i z_ij². Le seuil sur la variance empirique est
θ_z=100 ε_machine max(1,max_j \widehat Var(z_j)),
dans les unités normalisées ; si \widehat Var(z_j)≤θ_z, gain nul et compteur de rejets. Ce seuil fixe ne sera pas ajusté au test. On rapporte aussi les gains non finis comme échecs ; ils ne sont pas silencieusement remplacés. Tous les gains, bases, moyennes et paramètres restent fixes à l'inférence.

Les variances empiriques et covariances scalaires centrées utilisent toutes le diviseur N_t−1 (N_t≥2), y compris le seuil θ_z. Le ratio appris reste inchangé car ses deux diviseurs se simplifient. Wiener nominal emploie \hat p_j=Σ_iβ_ij²/(N_t−1), variance entraînement des β_j, et la variance connue v_j=(U_rᵀC_εU_r)_{jj}. Il a les mêmes moyennes et bases que les filtres appris. Son coût comprend l'estimation des variances ; sa connaissance de la covariance du bruit est une hypothèse de calibration commune.

### 7. Onze méthodes et sept familles réglables

Pour h=d−μ_d, les méthodes restituent μ_c+v(h).

| Méthode | Reconstruction du coefficient centré | Réglage |
|---|---|---|
| Tikhonov L2 | (AᵀA+αI)^{-1}Aᵀh | 20 α |
| Tikhonov H1 | (AᵀA+αP̄)^{-1}Aᵀh | 20 α |
| TSVD | Σ_{j≤k}σ_j^{-1}(u_jᵀh)v_j | ≤20 rangs |
| Wiener | g_j=σ_j\hat p_j/(σ_j²\hat p_j+v_j) | fixe |
| Wiener plafonné | min(g_j,L) | 20 L |
| Appris signé | \hat g_j | fixe |
| Appris signé plafonné | clip(\hat g_j,−L,L) | 20 L |
| Appris positif | max(0,\hat g_j) | fixe |
| Appris positif borné | clip(\hat g_j,0,1/σ_j) | fixe |
| Appris positif borné plafonné | clip(\hat g_j,0,min(1/σ_j,L)) | 20 L |
| LMMSE complet | \hat C_cd(\hat C_dd+λs_C I)^{-1}h | 20 λ |

Les quatre méthodes fixes et sept familles réglables totalisent au plus 144 candidats par condition. Les gains sont nuls hors du rang r pour les familles spectrales. Tikhonov complet emploie toutes les colonnes ; son inverse régularisé n'est pas remplacé par un TSVD. H1 ne commute généralement pas avec AᵀA : aucune diagonale artificielle dans la SVD ordinaire ne sera utilisée. Une résolution par factorisation, jamais une inverse explicite, est prescrite.

Pour LMMSE, \hat C_cd=Σ(c_i−μ_c)(d_i−μ_d)ᵀ/(N_t−1), \hat C_dd=Σ(d_i−μ_d)(d_i−μ_d)ᵀ/(N_t−1), et s_C=max(tr(\hat C_dd)/m,10^{-15}). L'entraînement utilise une observation bruitée par objet, la même pour toutes les méthodes. Les covariances sont estimées sur ces paires seulement. Résoudre le système régularisé symétrique par Cholesky ; un échec de factorisation est consigné, pas corrigé par une pénalité choisie sur le test. Pour λ>0 il est défini en arithmétique exacte. Coût : O(N_t(nm+m²)) pour les moments, O(m³) par facteur, O(nm) par inférence après pré-calcul de R. Les corrélations peuvent prédire une composante du noyau sans rendre le problème identifiable.

En population, l'optimalité affine matricielle suit de la dérivation de tr(C_cc)−2tr(RC_dc)+tr(RC_ddRᵀ), donnant RC_dd=C_cd. La régularisation empirique introduit un biais et ne conserve pas automatiquement cette optimalité exacte. Comparer au même budget de données et de validation.

Les grilles sont fixes : α=10^{linspace(−10,1,20)} pour L2 et H1 ; L=10^{linspace(0,4,20)} pour les trois plafonds ; λ=10^{linspace(−10,1,20)} pour LMMSE ; rangs distincts obtenus en arrondissant linspace(1,r,20), triés et dédupliqués. Si r=0, TSVD restitue la moyenne avec un candidat. Pas d'extension. Un optimum au bord est un diagnostic de grille, conservé comme limite.

### 8. Objets, graines et validation observationnelle

Un profil de K coefficients est généré par cₖ=Zₖ(1+k²)^{-p/2}/√(Σ_{l=0}^{K−1}(1+l²)^{-p}), Zₖ indépendants N(0,1). Ainsi E||q||²=1 dans chaque classe finie, sans renormalisation individuelle. p=2,5 nominal ; p=1,5 rugueux. Cette loi, surtout l'énergie du mode constant, est une hypothèse forte ; le bénéfice ML ne sera pas extrapolé à des traces arbitraires.

Trois entraînements disjoints de 1000 objets, un jeu validation de 200 et deux jeux test de 200 (nominal/rugueux), tous disjoints. La hauteur changée réutilise les objets test nominaux pour isoler la géométrie. Les objets pilote sont entièrement distincts. Pour éviter l'ordre d'exécution comme source d'aléa, prévoir des générateurs PCG64 initialisés par SeedSequence([20261011,partition,jeu,objet,canal,repetition]). Partition : pilote=0, train=1, validation=2, test nominal=3, test rugueux=4 ; canal coefficients=0, bruit=1, bootstrap=2. Aucun identifiant pilote ne migre au test. Les configurations gardent ces tuples. Les indices objet commencent à 0. Pour train, jeu=0,1,2 ; pour validation et les deux partitions test partagées, jeu=0, indépendamment de l’entraînement évalué. Pour pilote, jeu=0 à chaque reprise corrective afin de garder les mêmes objets. Pour coefficients, repetition=0 ; pour bruit train, repetition=0 ; pour validation, repetition=0,1 ; pour test, repetition=0,…,9 (0,…,4 au repli), plus repetition=10 réservée au résidu sur observation indépendante (5 au repli). Les masques, ν, méthodes et hauteurs ne figurent pas dans la graine : ils partagent la même normale standard de dimension 64, ensuite mise à l’échelle par τ de l’entraînement. Le bruit pilote possède deux réalisations 0,1 et une répétition 2 réservée au résidu indépendant.

Pour le bootstrap uniquement, partition=5, canal=2, jeu=indice entraînement 0,1,2, objet=0, repetition=b avec b=0,…,1999 ; PCG64 tire N_test indices uniformes entiers dans [0,N_test). Ce même vecteur est réutilisé pour P1, P2 et les médianes nominales associées, avec toutes les répétitions bruitées des objets sélectionnés. Partition=5 ne produit jamais un profil physique. Versions NumPy, BitGenerator, SeedSequence et configuration sont conservées ; ces conventions rendent les entrées régénérables, sans promesse de reproductibilité binaire entre bibliothèques différentes.

Pour chaque entraînement, masque et ν, les moyennes et filtres sont ajustés ; chaque hyperparamètre minimise
\[
Q(\lambda)=N_v^{-1}\sum_i\|A\widehat c_\lambda(d_i^{(1)})-d_i^{(2)}\|_2^2.
\]
Le direct nominal A est utilisé même si le générateur est exact/finer. Conditionnellement à l'objet, l'indépendance de la seconde répétition donne E||A\hat c−d_clean−ε₂||²=E||A\hat c−d_clean||²+tr Cε. Cette identité mesure la prédiction des observations vraies, pas ||\hat c−c||² ; le résidu de modèle reste présent. Les vérités de validation ne sont pas consultées par la règle, ni pour départager un ex aequo.

Égalité numérique : scores à distance ≤10^{-12}max(1,|Q_min|) sont ex aequo ; choisir le reconstructeur de plus petite norme d'opérateur calculée, puis l'indice de grille le plus petit. Cette règle est fixée avant test. Aucun oracle n'est prévu dans le classement principal.

Les entraînements ont des profils et bruits différents ; il ne s'agit pas de relancer trois fois un estimateur déterministe sur les mêmes données. Test/validation communs aux trois entraînements permettent une comparaison de sensibilité au jeu appris. Les scores ne sont pas artificiellement présentés comme provenant de 600 objets test indépendants.

### 9. Direct indépendant, tolérances et contrôles futurs

Le solveur différences finies utilisera un stencil centré d'ordre deux, des nœuds fantômes réfléchis pour Neumann homogène et une trace supérieure Dirichlet imposée. Les coins supérieurs reçoivent Dirichlet ; les coins inférieurs cumulent les réflexions latérales/inférieures. L'intégration des capteurs découpe exactement l'interpolant linéaire à leurs bornes, même si un capteur est plus fin qu'une maille.

Contrôles prévus, à exécuter seulement en phase B :
1. φ₀ : direct et moyenne constants, tolérance absolue 10^{-10}, puis modes k=1,2,4,8. Cette tolérance du constant couvre l'arrondi d'une résolution de faible taille, non une erreur de discrétisation non nulle.
2. Moyennes analytiques contre une quadrature indépendante raffinée ; cible absolue 10^{-10} pour les intégrales de modes k≤8. Vérifier masque et poids sur de petites données choisies à la main.
3. Solveur indépendant sur 65×17,129×33,257×65 : accord relatif en norme vectorielle avec la solution modale DISCRÈTE ≤10^{-10} pour k=1,…,8, constant selon contrôle 1. Séparément, erreur relative envers le transfert CONTINU ≤0,005 sur 257×65 ; au repli, ≤0,02 sur 129×33. Ces enveloppes sont justifiées ci-dessous avant implémentation. Sur les deux niveaux les plus fins, exiger diminution et ratio grossier/fin entre 2 et 6 pour les erreurs >10^{-10}. Aucune division cellule par cellule près des zéros ; employer la norme euclidienne du vecteur des 64 moyennes, puis contrôler les masques. Cette plage de ratios n’est pas une preuve universelle d’ordre.
4. Valeurs singulières, seuil numérique et identités de reconstruction SVD : résidu relatif ≤10^{-12} sur les petites matrices ; vérifier indépendamment noyau réel et directions tronquées.
5. Tikhonov sur un exemple 2×2 dont P ne commute pas avec AᵀA ; comparer aux équations normales, pas à une réimplémentation du même filtre erroné.
6. Cas scalaires de gain appris/Wiener : variances données, signal nul, bruit nul, moyenne non nulle, covariance négative et dénominateur faible ; contrôle algébrique sans ajustement au résultat.
7. LMMSE sur un cas corrélé 2×1 où le second coefficient appartient au noyau ; vérifier offset et coût. Norme et risque comparés aux formules indépendantes.
8. Identifiants de partitions disjoints et test protégé contre tout accès par apprentissage/réglage ; sérialisation des seeds et reprise atomique.

#### Justification a priori des tolérances et du rang d’acquisition

Avec N_x intervalles horizontaux, N_y verticaux, h_x=1/N_x, h_y=H/N_y et s=kπ/N_x, le stencil appliqué au cosinus donne cosh γ=1+2(h_y/h_x)²sin²(s/2), soit γ=2 asinh((h_y/h_x)sin(s/2)). La réflexion Neumann inférieure impose une dépendance cosh(jγ) ; le transfert nodal discret est sech(N_yγ). Lorsque les bornes des capteurs coïncident avec les nœuds (N_x multiple de 64), sommer exactement l’interpolant linéaire donne un facteur B(s)=(s/2)cot(s/2) relativement à la moyenne cosinus exacte. L’égalité vaut entre vecteurs même si une cellule a une moyenne nulle. Ainsi

ρ_k = cosh(kπH)/cosh(N_yγ) × (s/2)cot(s/2), et erreur relative = |ρ_k−1| ; ρ₀=1 par limite.

Pour H=1/4, N_y=N_x/4, l’évaluation scalaire de cette formule (pas une expérience du projet) donne au mode 8 :

| Points | Erreur relative prévue k=8 |
|---|---:|
| 65×17 | 0,0685089686443 |
| 129×33 | 0,0170127428900 |
| 257×65 | 0,00424583589265 |

L’évaluation des huit entiers k=1,…,8 donne son maximum en k=8 aux deux maillages fins. Les enveloppes 0,005 et 0,02 laissent respectivement environ 17,8 % et 17,6 % au-dessus de cette prédiction, sans modifier la tolérance d’accord à la solution discrète. Pour h_x=h_y, γ=s−s³/12+O(s⁵) et ρ_k−1=s²[(kπH)tanh(kπH)−1]/12+O(s⁴) : ordre deux à k fixé. Les ratios au mode 8 sont environ 4,027 et 4,007. Ces nombres ne démontrent pas que le futur solveur passe les tests.

Au repli inverse 33×9, les 32 capteurs contigus sur [0,1/2] ne dépendent que de 17 nœuds inférieurs. Les deux demi-cellules par maille ont coefficients (3/4,1/4) et (1/4,3/4) ; leur somme et leur différence permettent de retrouver les deux valeurs nodales. La matrice d’acquisition sur ces 17 nœuds a donc rang exact 17. Le rang du direct composé est au plus 17, jamais automatiquement min(m,n)=32 ; son rang numérique peut être inférieur. Cette grille grossière n’a pas de bornes toutes alignées aux capteurs : utiliser l’intégration par morceaux, sans lui appliquer le facteur B comme formule générale. Les grilles de contrôle 65×17 et129×33 restent alignées. Rapporter séparément rang acquisition, rang algébrique composé lorsqu’établi et rang numérique mesuré.

Les modes hauts, potentiellement plus mal résolus, sont examinés en erreur absolue normalisée et via ||Q(G_n−A_raw)||₂. Une faible erreur relative sur les modes bas ne certifie pas les hautes fréquences. Le générateur exact permet ce diagnostic sans confondre queue et discrétisation.

Pour calculer cosh ratios sans débordement, utiliser une expression exponentielle stable ; le pilote vérifiera l'accord avec la formule directe sur les modes sûrs. Une divergence, une mémoire excessive ou un défaut de partition bloque le passage. Les tolérances ne sont pas élargies pour faire passer une erreur.

### 10. Budget fermé, pilote et repli

Le plan complet comporte 3 entraînements × 3 masques × 3 bruits =27 conditions d'ajustement. Sept familles ×20 candidats +4 fixes : ≤3888 évaluations de configuration, chacune sur 200 objets de validation. Quatre régimes test ×200 objets ×10 répétitions ×11 méthodes ×27 conditions : 2 376 000 reconstructions. Les contrôles sans bruit utilisent seulement 16 profils pilotes préfixés et les modes analytiques, sans grille supplémentaire.

Le stockage comprend coefficients, observations ou seeds suffisants, matrices/filtres, erreurs par objet/répétition/méthode et agrégats. Ne pas conserver tous les champs intérieurs ou toutes les reconstructions de traces. Réaliser une estimation d'octets avant allocation ; écrire en lots atomiques identifiés. Les documents tiers ne sont pas inclus dans les sorties expérimentales.

#### Budget de stockage préfixé, en octets décimaux

Le plafond vaut 200 000 000 octets, sans supposer de compression. Les cinq métriques float64 par reconstruction sont erreur L2 absolue (queue incluse), erreur visible, erreur du complément orthogonal, résidu nominal et résidu indépendant. L’erreur relative est dérivée de l’absolue et de la norme vraie, avec d₀ enregistré ; amplification et paramètres sont par méthode/condition. Axes fixes : entraînement, masque, ν, régime, objet, répétition, méthode. Aucune chaîne ni dictionnaire dupliqué par ligne.

| Tableau ou groupe borné | Forme/type ou règle | Octets maximaux |
|---|---|---:|
| Chacun des 5 tableaux de métriques | 2 376 000 float64 | 19 008 000 chacun ; 95 040 000 ensemble |
| Statut par reconstruction | 2 376 000 uint8, bitmask motifs numérique/métrique/famille/manquant/corrompu | 2 376 000 |
| Coefficients vrais train/validation/tests | (3000+200+400)×128 float64 | 3 686 400 |
| Chaque matrice R sélectionnée | 64×64 float64 maximum, 297 matrices | 32 768 chacune ; 9 732 096 ensemble |
| Autres matrices, gains, moyennes, normes vraies | tableaux numériques inventoriés avant allocation | 2 267 904 |
| Scores/normes/statuts de validation | 3888 float64 par score ou norme ; autres résumés inclus | 1 000 000 |
| Bootstrap | 3×4×2000 float64 et statuts | 1 000 000 |
| Échantillons de sorties brutes en échec | au plus 1024, chacun ≤64 coefficients+64 observations float64 et identifiant uint64[4] | 1 081 344 |
| Pilotes et deux reprises, données et résultats | cumul, pas un plafond par essai | 4 000 000 |
| Configurations, seeds, provenance, entêtes et journaux | cumul ; motifs codés et compteurs, pas une trace texte par échec | 3 000 000 |
| Figures et données de courbes retenues | cumul | 10 000 000 |
| Écriture atomique temporaire | un lot à la fois | 1 000 000 |
| Réserve corrections documentées et autres sorties | aucune allocation implicite au-delà | 30 000 000 |
| TOTAL plafonds | inclut réserves et fichiers temporaires | 164 183 744 |

Les tableaux d’un repli sont plus petits ; ces plafonds restent valables. Ne pas stocker les champs intérieurs complets, tous les reconstructeurs candidats ou toutes les observations bruitées : les régénérer à partir des coefficients, des tuples de graines et des versions. Le résidu indépendant utilise la répétition réservée au §8, non une répétition ayant servi à reconstruire.

Pour les échecs, toutes les lignes de statut demeurent, même si elles sont toutes en échec ; les bits définissent une table finie de motifs. Les valeurs des métriques défaillantes valent +∞ et le statut manquant demeure distinct. L’échantillon brut retient les 1024 premiers identifiants en ordre lexicographique global (phase, condition, candidat ou méthode, objet, répétition), mises à jour bornées ; jamais les 1024 premiers selon l’ordre d’exécution. Les éléments non retenus gardent leur provenance et leurs graines, mais ne sont pas annoncés comme archivés en brut. Une exception de code conserve une trace par signature dans le budget de journaux, plus compteur. Toute limite de journal atteinte est signalée par compteur ; aucune ligne scientifique de statut n’est supprimée.

Avant chaque tableau, calculer produit(shape)×itemsize, ajouter entête réelle, fichiers existants, réserve et temporaire ; vérifier le quota de son groupe et le plafond total. Les sorties utiles invalidées restent comptabilisées et conservées : une nouvelle campagne complète ne peut pas être lancée si conserver l’ancienne et la nouvelle dépasse le plafond. Arrêter avant allocation et déclarer le blocage, sans effacement silencieux. La borne mémoire de travail indépendante reste 2 Gio.

Pilote : n=32,K=64, 24 objets entraînement, 8 validation, 8 test pilote, deux réalisations bruitées par objet test pilote (un seul niveau ν), masques complet/contigu et ν=0,01 ; méthodes et règles principales inchangées. Il contrôle aussi les factorisations du maillage fin complet et chronomètre un lot représentatif, sans consulter les futurs objets test. Limite de chaque essai : 600 s écoulées ; au plus deux reprises correctives documentées. Les grilles de 20 sont conservées afin de mesurer le coût de réglage.

Le budget CPU total de 14400 s inclut pilotes et reprises, génération, entraînement, réglage, inférence et analyse expérimentale. Répartition de conception : préparation/direct 1800 s, pilotes 1800 s au total, apprentissage/réglage 3600 s, test 3600 s, analyse 1800 s, réserve 1800 s. La limite CPU agrégée prévaut sur les 600 s écoulées par essai. Un thread BLAS par processus et l'absence de campagnes simultanées facilitent le contrôle ; comptabiliser les enfants éventuels. Arrêt avant allocation au-delà de 2 Gio et avant sortie au-delà de 200 Mo.

Passage campagne seulement si : tests valides ; sauvegarde/reprise éprouvée ; extrapolation CPU de la totalité des lots restants ×1,5 plus CPU déjà consommé ≤14400 s ; mémoire prévue <2 Gio ; sorties prévues <200 Mo. Les facteurs d'extrapolation incluent préparation et méthodes complètes, pas seulement les gains. Le facteur 1,5 est une marge d'ingénierie, pas une garantie.

Si le complet échoue à ce critère de faisabilité, adopter le repli de DECISION_V2 avant test : n=32,K=64, inverse 33×9, contrôles 65×17 et129×33, 300/60/60+60 objets, cinq bruits, trois régimes (nominal,rugueux,H×1,05), mêmes 27 conditions et onze méthodes. Total test 267300 reconstructions ; validation ≤3888×60 objets. Si le repli dépasse aussi le budget, ne pas réduire encore silencieusement ni ouvrir les tests : diagnostiquer le blocage. Le pilote peut constater une tolérance incompatible ; toute correction doit être justifiée mathématiquement et versionnée avant gel, sans chercher une amélioration ML.

### 11. Critères, incertitude et cas d'échec

Erreur principale : e=||\hat q−q||_{L2}/max(||q||_{L2},d₀), avec d₀=10^{-12} médiane des normes d'entraînement. La queue n,…,K−1 entre dans le numérateur. Rapporter également erreur absolue, erreur coefficient visible/projetée sur V_r, complément orthogonal, résidu nominal et résidu sur une répétition indépendante, norme d'amplification dans les unités déclarées, fréquence e>1, nombres d'échecs non finis et paramètres au bord de grille.

Deux comparaisons principales préfixées, ν=0,01 :
P1 : appris signé plafonné versus appris signé, masque contigu, H nominal ;
P2 : mêmes méthodes, capteurs complets, hauteur vraie 1,05H.
Pour chaque jeu d'entraînement, calculer ratio des quantiles 0,95 de e (plafonné/non plafonné) et ratio des médianes nominales correspondantes (même masque, H nominal). Un dénominateur nul donne comparaison indéfinie, pas victoire automatique.

Bootstrap apparié de 2000 tirages des objets physiques avec toutes leurs réalisations bruitées ; réutiliser les mêmes objets tirés pour les deux méthodes et les scénarios liés. Quantile empirique défini par la statistique ordonnée de rang ceil(αN), sans interpolation, dans les réels étendus [0,+∞]. Pour quatre ratios (deux contrastes ×deux critères), utiliser des bornes supérieures bootstrap à 98,75 %, correction de Bonferroni nominale à 5 % pour les quatre ratios CONDITIONNELLEMENT À UN entraînement. Aucune couverture simultanée à 95 % des trois entraînements (douze ratios) n’est revendiquée. Ces intervalles sont approximatifs, conditionnels au jeu appris, pas des garanties finies exactes.

Le critère est satisfait pour un entraînement si les bornes supérieures sont <0,90 pour le ratio du q95 et <1,05 pour la médiane nominale. Rapporter séparément les trois entraînements et leur étendue. Une conclusion robuste au jeu d'entraînement exige la satisfaction dans chacun ; trois jeux ne permettent pas de revendiquer une loi précise de généralisation sur tous les entraînements possibles. Ne pas traiter leurs résultats communs de test comme indépendants.

Analyses secondaires : comparaison à Wiener/Wiener plafonné et LMMSE ; effet de positivité puis borne 1/σ ; rugosité ; H×1,02 ; autres bruits ; erreur de modèle et de troncature ; coûts. Aucun gain secondaire ne remplace un échec principal. Pas d'énoncé de supériorité globale à partir d'un seul contraste.

Illustrations fixées : pour chaque contraste principal et le premier entraînement, choisir les objets de rang médian, au plus proche du q95 et pire selon la moyenne de e sur bruits pour le filtre appris non plafonné ; égalités départagées par identifiant. Montrer ensuite toutes les méthodes pour ces mêmes objets, avec répétition d'indice 0. Ajouter les cas non finis s'ils existent et un cas de résidu faible/erreur élevée selon les déciles préfixés (résidu ≤décile1, e≥décile9), sinon noter son absence. Aucun choix manuel de « belles » courbes.

#### Traitement préfixé des échecs, valeurs non finies et données manquantes

Une reconstruction comportant NaN ou une valeur infinie dans ses coefficients est un échec numérique. Un calcul de métrique qui déborde malgré des coefficients finis est identifié séparément comme échec du calcul de métrique. Dans les deux cas, conserver chaque ligne de statut et son motif, avec méthode, condition, objet, répétition, graines et version retrouvables par les axes et le manifeste ; seules les sorties brutes échantillonnées dans le plafond du §10 sont conservées, sans collection illimitée ; poser e=+∞ pour les distributions d’erreur. Ne jamais supprimer la ligne, imputer une erreur favorable ou renormaliser les poids sur les seuls succès. Rapporter taux d’échec numérique, taux d’échec de métrique et leur union, avec dénominateur égal au nombre de reconstructions prévues. La fréquence e>1 inclut ces échecs. Un résumé conditionnel aux seuls succès peut être montré en analyse secondaire explicitement étiquetée ; il ne sert jamais au contraste principal.

Pour N erreurs ordonnées dans [0,+∞], Q_α=e_(ceil(αN)). Cette définition couvre exactement médiane et q95, sans opération ambiguë telle que 0×∞. La moyenne des erreurs d’un objet vaut +∞ dès qu’une de ses répétitions échoue ; les objets ex aequo infinis sont départagés par identifiant pour les figures, et le nombre total d’échecs reste affiché. Un résidu non fini ne satisfait pas le critère secondaire « résidu faible ».

Un ratio principal n’est interprétable comme une réduction relative que si son dénominateur est strictement positif et fini. Si son numérateur est +∞ et le dénominateur fini positif, le ratio vaut +∞ et le critère de succès est non satisfait. Si le dénominateur vaut zéro ou +∞, le ratio est indéfini, y compris lorsque le numérateur est fini : ne pas annoncer une amélioration de 100 % sur une référence en échec. Le verdict correspondant est NON CONCLUANT POUR LE CONTRASTE RELATIF ; les valeurs absolues, taux d’échec et comparaisons finies secondaires demeurent rapportés. Si les deux quantiles sont infinis, aucun classement par leur ratio n’est produit.

Le bootstrap rééchantillonne toutes les lignes par groupes d’objets, y compris les échecs. Aucun tirage indéfini ou infini n’est retiré ni remplacé. Si un des 2000 tirages rend un ratio indéfini, rapporter son nombre et ne pas produire de borne percentile interprétée comme preuve du succès : verdict NON CONCLUANT pour ce contraste. Sinon, conserver les ratios +∞ dans leur ordre et calculer la borne percentile par rang ceil(0,9875×2000). Si cette borne vaut +∞, le critère de succès est non satisfait. Cette règle volontairement prudente évite de sélectionner les seuls tirages favorables ; elle ne transforme pas un bootstrap approximatif en garantie exacte.

Pendant la validation, un candidat produisant une reconstruction ou un score non fini reçoit Q=+∞ et est marqué invalide. Il ne peut pas gagner contre un candidat valide. Si tous les candidats d’une famille sont invalides, ne pas choisir arbitrairement le premier : conserver l’échec de cette famille pour la condition et compter ses reconstructions test prévues comme échecs de méthode si une campagne déjà engagée atteint ce cas. Au pilote, tout échec numérique inattendu bloque le passage jusqu’à diagnostic dans les limites de reprises autorisées. Une famille fixe aux gains non finis est traitée de la même manière, sans remplacement silencieux par zéro ou par une autre méthode.

Une sortie absente par interruption n’est pas assimilée à une reconstruction numérique ratée. Elle constitue un lot manquant, à reprendre seulement si le budget le permet. Tant que des lots prévus restent manquants, la campagne est INCOMPLÈTE et aucun verdict principal global n’est annoncé. Les échecs numériques observés, les lots non exécutés et les données corrompues sont trois états distincts. Toute correction scientifique invalidant des lots est documentée ; les sorties initiales utiles et leur motif d’invalidation sont conservés.

### 12. Sources, corrections et limites

C1–C5 et R01–R14 sont définis dans sources.md. Les formules ci-dessus sont redérivées ; leurs origines classiques sont explicites. Les références R04/R05/R09 d'accès incomplet ne fournissent aucune garantie importée.

| DECISION.md initial | Correction V2 |
|---|---|
| Sech près de la SVD discrète | Deux opérateurs, moyenne exacte, quadrature et masque distingués |
| Borne sans centrage explicite | Offset f, noyau numérique et projection ajoutés |
| Gain supervisé sans traitement des moyennes | Classe affine, seuil de variance, covariance signal–bruit |
| Contrainte positive confondue avec plafond | Trois ablations distinctes et plafond signé principal |
| Wiener seul comparateur statistique | LMMSE complet entraîné sur mêmes paires |
| H1 insuffisamment spécifiée | P complète dans les équations normales |
| Bruit « pour chaque masque » ambigu | Calibration 64 capteurs entraînement, bruit après moyenne avant masque |
| Générateur numérique fin assimilable à vérité | Générateur analytique exact et contrôle différences finies indépendant |
| Grille 20 puis extension | 20 au total ; aucune extension |
| Budget qualitatif | Dénombrement candidats/reconstructions et enveloppe CPU fermée |
| Trois graines sans objet défini | Trois entraînements différents, tests partagés et bootstrap par objet |

Aucune donnée clinique, supériorité mesurée, garantie d'erreur de modèle réelle, convergence infinidimensionnelle ou nouveauté algorithmique n'est établie. Le manuscrit final devra séparer les preuves analytiques, diagnostics numériques, estimations statistiques et hypothèses physiques.

### Précisions de sélection des diagnostics

Les diagnostics illustratifs P1 et P2 utilisent uniquement l'entraînement 0, le bruit de 1 % et le filtre appris signé non plafonné comme référence de sélection. Pour chaque objet, l'erreur relative et le résidu nominal sont moyennés sur les dix répétitions ; un échec est conservé comme +∞. Les objets médian, q95 et pire sont choisis par rang supérieur sur les 200 moyennes d'erreur, avec identifiant croissant en cas d'égalité. Toutes les méthodes sont ensuite représentées sur le même objet, réalisation 0.

Pour le diagnostic « faible résidu, erreur élevée », calculer séparément pour chaque contraste les déciles D1 du résidu et D9 de l'erreur sur tous les objets, sans interpolation ni retrait des infinis. Parmi les objets dont les deux moyennes sont finies, retenir ceux vérifiant résidu ≤ D1 et erreur ≥ D9 ; choisir l'erreur décroissante, puis le résidu croissant, puis l'identifiant croissant. Au plus un objet est retenu par contraste. Une absence est rapportée.

Une éventuelle figure d'échec utilise le premier couple (objet, répétition) lexicographique où au moins une méthode échoue, au plus un panneau par contraste. Les courbes non finies ne sont pas remplacées par des courbes inventées. Six illustrations de quantiles, deux diagnostics et deux panneaux d'échec au maximum donnent un plafond de dix figures de traces ; les figures agrégées sont distinctes. Toutes les figures et données de courbes partagent le plafond de 10 Mo. Ces règles ont été précisées avant l'ouverture du test principal.

Dans le contrôle d'acquisition du repli, deux moyennes successives vérifient a=(3u+v)/4 et b=(u+3v)/4. Les deux valeurs nodales sont u=(3a−b)/2 et v=(3b−a)/2 ; la différence seule ne suffit pas. Le rang 17 annoncé concerne l'acquisition contiguë du repli, non un rang minimal imposé à sa SVD numérique.

Le pilote utilise la partition 0 : objets 0–23 pour l’entraînement, 24–31 pour la validation, 32–39 pour le test pilote et 40–55 pour les contrôles sans bruit. Ils ne sont pas réutilisés dans les partitions principales.
