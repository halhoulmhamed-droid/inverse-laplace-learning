# Théorie : complétion harmonique et filtres statistiques

**Projet scientifique de Master — Halhoul Mhamed**<br>
Université Abdelmalek Essaâdi — Faculté Polydisciplinaire de Larache<br>
Master « Modélisation mathématique et applications et apprentissages » — M2 en cours<br>
**Professeur encadrant : Pr. Chakir Tajani**

Cette note établit les propriétés du modèle employé. Ses preuves analytiques ne sont pas des résultats expérimentaux. Le rapport de campagne distingue réussite des contrôles numériques et réussite des critères comparatifs : un statut d'exécution `PASS` ne démontre aucune supériorité statistique.

## 1. Inconnue, conditions aux limites et solution

Sur $\Omega_H=(0,1)\times(0,H)$, avec $H=1/4$, on cherche

$$
\Delta u=0,\qquad u(x,H)=q(x),\qquad
\partial_yu(x,0)=0,\qquad \partial_xu(0,y)=\partial_xu(1,y)=0.
$$

La trace supérieure $q$ est inconnue ; les observations sont des moyennes de la trace inférieure. Le flux inférieur nul est une hypothèse physique contrôlée, pas une seconde mesure supposée sans bruit. Pour une différence de solutions de même trace supérieure, l'identité de Green donne $\int_{\Omega_H}|\nabla w|^2=0$ ; la constante restante est nulle au bord supérieur. Le problème direct est donc unique dans son cadre énergétique.

Dans la base orthonormée de $L^2(0,1)$,

$$
\phi_0=1,\quad \phi_k(x)=\sqrt2\cos(k\pi x),\qquad
u_k(x,y)=\phi_k(x)\frac{\cosh(k\pi y)}{\cosh(k\pi H)}\quad(k\ge1),
$$

et $u_0=1$. Chaque mode vérifie directement l'équation et les conditions. Pour une somme finie $q=\sum c_k\phi_k$, la solution est $u=\sum c_ku_k$. Pour $q\in H^{1/2}$, on adopte la solution faible usuelle ; pour $q\in L^2$, le transfert inférieur est défini par prolongement de la série, sans affirmer une régularité $H^1$ jusqu'au bord supérieur pour toute donnée $L^2$.

Le transfert continu $T_H$ satisfait

$$
T_H\phi_k=t_k\phi_k,\qquad t_0=1,\quad t_k=\operatorname{sech}(k\pi H)>0.
$$

Il est injectif, compact et exponentiellement lissant. La suite $q_k=\phi_k$ a norme unité tandis que $\|T_Hq_k\|=t_k\to0$ : l'inverse n'est pas borné sur son image. Une donnée de coefficients $d_k$ appartient à cette image exactement lorsque $\sum_k|d_k/t_k|^2<\infty$. Des données bruitées arbitraires ne satisfont pas nécessairement cette condition. Le mode constant est conservé ; le centrage statistique ultérieur ne le supprime pas.

## 2. Moyennes, masque et normalisation

Les 64 capteurs portent sur $I_j=[j/64,(j+1)/64]$, de longueur $\ell=1/64$. Pour $I_j=[a_j,b_j]$,

$$
S_jv=\ell^{-1}\int_{I_j}v(x)\,dx,\qquad
S_j\phi_k=\frac{\sqrt2[\sin(k\pi b_j)-\sin(k\pi a_j)]}{k\pi\ell},\quad S_j\phi_0=1.
$$

Le masque $M$ conserve tous les capteurs, les indices pairs, ou les indices 0 à 31. Avec $m$ observations conservées, $W=\ell I_m$ et $\|y\|_W=\|W^{1/2}y\|_2$. Les poids ne sont pas renormalisés lorsque la moitié du bord manque. Cette norme porte sur des moyennes, pas sur la trace continue exacte.

Le générateur analytique utilise $K=128$ coefficients ; l'inverse utilise $n=64$ coefficients et une matrice numérique $A_{\rm raw}\in\mathbb R^{m\times n}$. Définissons

$$
a=\|W^{1/2}A_{\rm raw}\|_2,\quad Q=W^{1/2}/a,\quad
A=QA_{\rm raw},\quad d=QM y_{\rm raw}.
$$

La SVD $A=U\Sigma V^T$ concerne cette matrice masquée et discrétisée. Ses singularités ne sont pas automatiquement les $t_k$. Le rang numérique retient $\sigma_j>10^{-12}\sigma_1$ ; les directions rejetées peuvent être physiquement non nulles. Le noyau exact et ce sous-espace numérique ne doivent pas être confondus.

Le bruit physique $\varepsilon_{\rm raw}\sim N(0,\tau^2I_{64})$ est ajouté **après** les moyennes et **avant** le masque. L'échelle $\tau=\nu s_{\rm train}$ utilise le RMS des observations propres des seuls objets d'entraînement. Elle reste fixe aux tests. Ainsi

$$
\varepsilon=QM\varepsilon_{\rm raw},\qquad
C_\varepsilon=\frac{\tau^2\ell}{a^2}I_m.
$$

Une norme de reconstruction $L$ exprimée pour $d$ correspond à $L/a$ pour les observations physiques pondérées, ou $L\sqrt\ell/a$ pour leur norme euclidienne brute.

## 3. Direct numérique indépendant

Le direct inverse emploie des différences finies centrées, des nœuds fantômes réfléchis aux bords de Neumann et une trace supérieure imposée. Les coins supérieurs sont Dirichlet ; les coins inférieurs cumulent les deux réflexions. Les moyennes de la trace calculée intègrent exactement son interpolant linéaire, avec subdivision aux bornes des capteurs. Le générateur analytique ne réutilise pas cette matrice.

Pour $N_x$ et $N_y$ intervalles, $h_x=1/N_x$, $h_y=H/N_y$ et $s=k\pi/N_x$, le stencil modal donne

$$
\cosh\gamma=1+2(h_y/h_x)^2\sin^2(s/2),\qquad
\gamma=2\operatorname{asinh}((h_y/h_x)\sin(s/2)).
$$

La réflexion inférieure impose la dépendance $\cosh(j\gamma)$ ; le transfert nodal vaut $\operatorname{sech}(N_y\gamma)$. Si les capteurs sont alignés aux nœuds, l'intégration linéaire apporte le facteur $B(s)=(s/2)\cot(s/2)$ relativement à la moyenne analytique. Le ratio des vecteurs de moyennes est donc

$$
\rho_k=\frac{\cosh(k\pi H)}{\cosh(N_y\gamma)}B(s),\qquad \rho_0=1.
$$

Il ne faut pas diviser cellule par cellule près des zéros. Pour $h_x=h_y$,

$$
\rho_k-1=\frac{s^2}{12}\{(k\pi H)\tanh(k\pi H)-1\}+O(s^4).
$$

Au mode 8, les erreurs relatives prédites sur les grilles de points 65×17, 129×33 et 257×65 sont respectivement 0,0685089686443, 0,0170127428900 et 0,00424583589265. Les rapports successifs sont environ 4,027 et 4,007. Les seuils fixés avant calcul distinguent accord à la solution **discrète** ($10^{-10}$) et erreur envers la solution **continue** (0,005 au maillage 257×65 ; 0,02 au repli 129×33). Une bonne précision des modes bas ne certifie pas les hautes fréquences.

Au repli inverse 33×9, les 32 capteurs contigus dépendent de 17 nœuds bas. Deux moyennes par maille valent $a=(3v_0+v_1)/4$ et $b=(v_0+3v_1)/4$ : $v_0=(3a-b)/2$, $v_1=(3b-a)/2$. Leur somme et leur différence établissent le rang d'acquisition 17 ; le rang du direct composé est au plus 17 et peut être numériquement inférieur. Le facteur $B$ ne s'applique pas directement à cette grille non alignée.

## 4. Identité d'erreur affine : offset, noyau et queue

Écrivons $q=J_nc+q_\perp$, où $J_n$ synthétise les $n$ premiers cosinus. Pour une hauteur vraie $H_*$,

$$
d=(A+E)c+b+\varepsilon,\qquad
E=Q(MST_{H_*}J_n-A_{\rm raw}),\qquad b=QMST_{H_*}q_\perp.
$$

$E$ sépare changement géométrique et erreur du direct inverse ; $b$ est la queue transmise. Les moyennes d'entraînement $\mu_c,\mu_d$ définissent l'estimateur

$$
\widehat c=\mu_c+R(d-\mu_d),\qquad f=A\mu_c-\mu_d.
$$

Une substitution exacte donne

$$
\widehat c-c=(RA-I)(c-\mu_c)+R(Ec+b+\varepsilon+f).
$$

Par conséquent,

$$
\|\widehat c-c\|\le\|(RA-I)(c-\mu_c)\|
+\|R\|\bigl(\|E\|\|c\|+\|b\|+\|\varepsilon\|+\|f\|\bigr).
$$

Pour $R=V_r\operatorname{diag}(g)U_r^T$, $\|R\|=\max_j|g_j|$. Dans une base $V$ complétée et avec $\beta=V^T(c-\mu_c)$, le premier terme au carré vaut

$$
\sum_{j\le r}(g_j\sigma_j-1)^2\beta_j^2+\sum_{j>r}\beta_j^2.
$$

Le second terme conserve noyau et directions rejetées. H1 et LMMSE utilisent leur véritable norme matricielle et leur biais $(RA-I)(c-\mu_c)$ ; une formule diagonale ne leur est pas attribuée.

L'orthogonalité apporte

$$
\|\widehat q-q\|_{L^2}^2=\|\widehat c-c\|_2^2+\|q_\perp\|_{L^2}^2.
$$

Le diagnostic du complément de $V_r$ porte sur les **n coefficients** ; la queue est ajoutée dans l'erreur totale. Pour une vérité finie $K$, aucune affirmation n'est faite sur une queue physique au-delà de $K$.

Conditionnellement à l'objet et à un estimateur fixé, si le bruit test est centré de covariance $C_\varepsilon$, poser $h=(RA-I)(c-\mu_c)+R(Ec+b+f)$ donne

$$
\mathbb E_\varepsilon\|\widehat c-c\|^2=\|h\|^2+
\operatorname{tr}(RC_\varepsilon R^T)
\le\|h\|^2+\|R\|^2\operatorname{tr}(C_\varepsilon).
$$

Cette identité requiert le centrage conditionnel du bruit. Le plafond borne le terme de variance, mais peut augmenter le biais ; il ne garantit donc pas une diminution du risque total. Pour un bruit gaussien, aucune constante finie ne borne uniformément toutes les réalisations. Une borne de Markov, avec probabilité au moins $1-\zeta$, est $\|\varepsilon\|\le(\tau\sqrt\ell/a)\sqrt{m/\zeta}$. Elle est conservatrice et non simultanée sur toute une campagne. Une erreur de discrétisation simplement mesurée ne devient pas un certificat rigoureux.

## 5. Gains appris, Wiener et contraintes distinctes

Pour une paire scalaire centrée $(\beta,z)$,

$$
\mathbb E(\beta-gz)^2=\operatorname{Var}\beta-2g\operatorname{Cov}(\beta,z)
+g^2\operatorname{Var}z.
$$

Si $\operatorname{Var}z>0$, le minimiseur est $g_* =\operatorname{Cov}(\beta,z)/\operatorname{Var}z$. Sinon le choix conventionnel est zéro. Dans le modèle exact $z=\sigma\beta+\varepsilon$, avec $p=\operatorname{Var}\beta$, $v=\operatorname{Var}\varepsilon$ et $t=\operatorname{Cov}(\beta,\varepsilon)$,

$$
g_* =\frac{\sigma p+t}{\sigma^2p+2\sigma t+v}.
$$

Wiener nominal correspond à $t=0$. Un opérateur erroné peut mélanger les modes ; le rapport de covariance demeure défini, sans nécessairement être Wiener nominal. Les gains peuvent être négatifs. Pour des moyennes non nulles, l'intercept affine est indispensable.

Les moments empiriques centrés utilisent tous le diviseur $N_t-1$. Le gain appris est nul si la variance observée est au plus $100\epsilon_{\rm machine}\max(1,\max_j\widehat{\operatorname{Var}}z_j)$. Les gains non finis sont des échecs explicites ; un plafond ne les transforme pas en candidats valides. Les ablations distinguent gain signé, partie positive, restriction $[0,1/\sigma_j]$, plafond signé $[-L,L]$, puis plafond appliqué au gain positif borné. Tous restent fixes à l'inférence.

## 6. Comparateurs quadratiques et matriciels

Tikhonov L2 et H1 résolvent respectivement

$$
R_{L2}=(A^TA+\alpha I)^{-1}A^T,\qquad
R_{H1}=(A^TA+\alpha\overline P)^{-1}A^T,
$$

avec $P=\operatorname{diag}(1+(k\pi)^2)$ et $\overline P=P/\max\operatorname{diag}P$. Ces écritures désignent des résolutions de systèmes, pas une inversion explicite. H1 ne commute généralement pas avec $A^TA$ après masque ; diagonaliser artificiellement sa pénalité dans la SVD ordinaire change le problème. Les systèmes complets conservent les colonnes numériquement rejetées par les familles spectrales.

TSVD retient les $k$ premières directions du rang numérique. Wiener et Wiener plafonné sont distincts des cinq variantes apprises. Le onzième estimateur, LMMSE, emploie

$$
R=\widehat C_{cd}(\widehat C_{dd}+\lambda s_CI)^{-1},\qquad
s_C=\max(\operatorname{tr}\widehat C_{dd}/m,10^{-15}),
$$

résolu par Cholesky. En population, dériver le risque matriciel donne $RC_{dd}=C_{cd}$. Une corrélation peut prédire un coefficient du noyau : si $c_2=2c_1$ et $d=c_1$, l'estimateur prédit les deux coefficients. Cela n'établit aucune identifiabilité déterministe hors de cette loi. Le ridge empirique introduit son propre biais. Les coûts des moments, factorisations et applications sont respectivement $O(N_t(nm+m^2))$, $O(m^3)$ et $O(nm)$.

## 7. Validation indépendante et interprétation

Toutes les méthodes partagent les objets, observations bruitées et moyennes d'entraînement. Les partitions sont séparées par objet avant bruit. Les vérités d'entraînement servent aux moments supervisés ; les paramètres sont choisis uniquement par

$$
Q=N_v^{-1}\sum_i\|A\widehat c(d_i^{(1)})-d_i^{(2)}\|^2,
$$

avec deux répétitions indépendantes conditionnellement au même objet. L'espérance conditionnelle du carré est l'erreur prédictive envers l'observation propre plus $\operatorname{tr}C_\varepsilon$. Ce risque n'est pas le risque inverse : une erreur du noyau peut rester invisible. Aucune vérité de validation ou de test ne départage les paramètres.

Les grilles sont fixées : 20 valeurs logarithmiques entre $10^{-10}$ et $10^1$ pour $\alpha,\lambda$, entre $1$ et $10^4$ pour $L$ ; au plus 20 rangs distincts pour TSVD. Sept familles réglables et quatre fixes donnent au plus 144 candidats par condition, sans extension. Les ex æquo de score à $10^{-12}\max(1,|Q_{\min}|)$ sont départagés par norme d'opérateur puis indice. Une famille entièrement invalide demeure en échec.

Les contrastes principaux comparent appris signé plafonné et non plafonné à 1 % de bruit : masque contigu nominal, puis masque complet avec hauteur vraie $1,05H$. Les ratios q95 et médianes nominales associés sont évalués par bootstrap apparié de 2000 tirages d'objets, gardant toutes leurs répétitions. Les quatre bornes supérieures à 98,75 % sont conditionnelles à **un** entraînement ; elles ne couvrent pas simultanément les trois entraînements. Les seuils sont respectivement 0,90 et 1,05. Les quantiles utilisent le rang supérieur sans interpolation ; échecs infinis et ratios indéfinis ne sont ni retirés ni assimilés à des victoires.

## 8. Portée et références

Le rectangle, les profils gaussiens à spectre prescrit, la calibration et le flux contrôlé constituent un banc synthétique. Ce travail ne valide ni mesure clinique, ni géométrie réaliste, ni stabilité uniforme au raffinement. Une absence de gain face à Wiener ou LMMSE est recevable. La régularisation spectrale apprise et l'apprentissage quadratique ont des antériorités ; aucune nouveauté algorithmique n'est revendiquée.

Les identifiants ci-dessous conservent ceux de l'étude bibliographique. Ils situent les méthodes ; les dérivations élémentaires ci-dessus sont explicites et ne prétendent pas reproduire toutes les preuves des publications.

- **C2.** Chakir Tajani et Jaafar Abouchabaka, *On the data completion problem for Laplace’s equation*, Annals of the University of Craiova 45(1), 11–36, 2018. Lien thématique avec la complétion de Cauchy ; aucun accord académique nouveau ni reproduction de son schéma itératif n'est supposé.
- **R01.** Martin Burger et Samira Kabri, *Learned Regularization for Inverse Problems: Insights from a Spectral Model*, [arXiv:2312.09845v2](https://arxiv.org/abs/2312.09845v2), 2024. Antériorité spectrale et dépendance aux moments des distributions.
- **R02.** Giovanni S. Alberti, Ernesto De Vito, Matti Lassas, Luca Ratti et Matteo Santacesaria, *Learning the optimal Tikhonov regularizer for inverse problems*, [NeurIPS 2021](https://papers.nips.cc/paper_files/paper/2021/hash/d3e6cd9f66f2c1d3840ade4161cf7406-Abstract.html). Notice et résumé examinés ; les garanties complètes ne sont pas transférées ici.
- **R06.** Michael J. Byrne et Rosemary A. Renaut, *Learning spectral windowing parameters for regularization using unbiased predictive risk and generalized cross validation techniques for multiple data sets*, [Inverse Problems and Imaging 17(4), 841–869, 2023](https://doi.org/10.3934/ipi.2023006). Motivation du réglage observationnel ; ne garantit pas son optimalité pour le risque inverse présent.

Ces références proviennent de la bibliographie préparatoire ; cette note n'en constitue pas une nouvelle vérification externe. Aucun PDF tiers n'est redistribué.

## Contact

halhoul.mhamed@etu.uae.ac.ma<br>
halhoulmhamed@gmail.com<br>
https://github.com/halhoulmhamed-droid/inverse-laplace-learning
