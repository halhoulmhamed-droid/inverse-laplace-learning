# Note de préparation à la soutenance

**Halhoul Mhamed — Projet scientifique de Master**<br>
Université Abdelmalek Essaâdi — Faculté Polydisciplinaire de Larache<br>
Master « Modélisation mathématique et applications et apprentissages » — M2 en cours<br>
Professeur encadrant : Pr. Chakir Tajani

## Présenter la question scientifique

L'objectif est de reconstruire un potentiel inaccessible à partir de mesures prises sur un autre bord, lorsque ces mesures sont bruitées, incomplètes et produites par un modèle légèrement différent de celui utilisé pour l'inversion. Le projet compare des régularisations classiques et des filtres spectraux appris. Il cherche notamment à savoir si limiter l'amplification d'un filtre appris améliore ses erreurs extrêmes, et à quel prix sur son erreur typique.

Il s'agit d'une étude synthétique sur un rectangle. Le filtrage appris, Wiener et la régularisation quadratique sont des méthodes établies : aucune invention de ces principes ni validation clinique n'est revendiquée. L'apport pédagogique consiste à rendre les hypothèses, les comparaisons et leurs limites vérifiables.

## Expliquer le modèle et les données

Dans le rectangle \(\Omega=(0,1)\times(0,H)\), avec \(H=0{,}25\), le potentiel vérifie \(\Delta u=0\). La trace supérieure \(q(x)=u(x,H)\) est inconnue. Le flux normal est nul en bas et sur les côtés. Cette condition représente une isolation supposée connue ; ce n'est pas une deuxième mesure miraculeusement exempte de bruit.

Le bord inférieur porte 64 capteurs. Chacun mesure une **moyenne sur un intervalle** de longueur \(1/64\), et non une valeur ponctuelle. Trois configurations conservent respectivement les 64 capteurs, les 32 d'indices pairs ou les 32 premiers capteurs contigus. Le bruit est ajouté aux 64 moyennes physiques avant de sélectionner les capteurs. Son échelle est calibrée sur l'entraînement, puis conservée au test.

La base \(\phi_0=1\), \(\phi_k=\sqrt2\cos(k\pi x)\) est orthonormée dans \(L^2(0,1)\). Pour une trace supérieure de coefficients \(c_k\), la solution modale est

\[
u(x,y)=c_0+\sum_{k\ge1}c_k\phi_k(x)
\frac{\cosh(k\pi y)}{\cosh(k\pi H)}.
\]

La substitution permet de vérifier l'équation, le flux nul et la trace supérieure. Au bord inférieur, chaque coefficient est multiplié par \(t_k=\operatorname{sech}(k\pi H)\). Ces facteurs tendent exponentiellement vers zéro. Des traces supérieures de norme unité peuvent donc produire des observations arbitrairement petites : inverser amplifie fortement le bruit. Le transfert continu est injectif, mais son inverse sur son image n'est pas borné. Avec un nombre fini de capteurs, la reconstruction fonctionnelle comporte en plus une non-unicité.

Les singularités de la **matrice discrète** ne sont pas simplement les \(t_k\) : moyennes, masque, discrétisation et normalisation interviennent. Les observations utilisent les poids \(1/64\), conservés même pour un demi-bord. Les coefficients physiques gardent leur norme euclidienne, égale à la norme \(L^2\) de leur synthèse. Le mode constant n'est pas supprimé par le centrage statistique.

## Montrer ce qui est appris

Pour la matrice normalisée \(A=U\Sigma V^\top\), un filtre diagonal reconstruit

\[
\widehat c=\mu_c+V_r\operatorname{diag}(g_j)U_r^\top(d-\mu_d).
\]

Les moyennes, la base et les gains sont estimés ou définis avant le test, puis restent fixes. L'apprentissage est ici une régression scalaire par direction singulière, sans réseau neuronal. Pour deux variables centrées, \(\beta\) et \(z\), minimiser \(\mathbb E(\beta-gz)^2\) donne

\[
g^*=\frac{\operatorname{Cov}(\beta,z)}{\operatorname{Var}(z)}.
\]

Sous le modèle exact \(z=\sigma\beta+\varepsilon\), avec signal et bruit décorrélés, ce gain devient celui de Wiener : \(\sigma p/(\sigma^2p+v)\). Un apprentissage qui retrouve Wiener est donc cohérent ; il ne constitue pas un échec. Les gains empiriques dépendent de l'échantillon et peuvent être négatifs. Positivité, borne \(1/\sigma_j\) et plafond global \(L\) sont des contraintes différentes, évaluées séparément.

Le plafonnement signé impose \(|g_j|\le L\). Il borne la norme du reconstructeur diagonal dans les unités normalisées et limite l'amplification des perturbations. Il peut aussi accroître le biais en atténuant des composantes utiles : aucun gain universel n'en découle.

Le LMMSE complet utilise les covariances croisées entre toutes les composantes. Il peut prédire statistiquement une direction non observable si l'entraînement la corrèle à une direction observée. Cela ne crée pas une identifiabilité déterministe. Tikhonov H1 utilise, lui, une pénalité liée aux dérivées de la trace ; après masquage, elle n'est généralement pas diagonale dans la SVD ordinaire. Pour H1 et LMMSE, l'amplification se mesure par la vraie norme de leur matrice, et non par un maximum de gains diagonaux fictifs.

## Décomposer l'erreur avant d'interpréter les figures

Écrivons \(d=(A+E)c+b+\varepsilon\), où \(E\) représente l'erreur d'opérateur et \(b\) les modes de vérité absents du modèle inverse. Pour tout reconstructeur affine de matrice \(R\), en posant \(f=A\mu_c-\mu_d\),

\[
\widehat c-c=(RA-I)(c-\mu_c)+R(Ec+b+\varepsilon+f).
\]

Le premier terme est le biais de reconstruction, y compris dans les directions rejetées. Le second rassemble erreur de modèle, troncature transmise, bruit et décalage des moyennes. Plafonner \(\|R\|\) agit sur ce second terme sans annuler le premier. L'erreur de la trace entière inclut aussi l'énergie des coefficients de vérité au-delà des 64 coefficients reconstruits. Les retirer artificiellement favoriserait les méthodes.

Le générateur emploie 128 modes et des moyennes analytiques. L'inverse utilise 64 modes et un direct par différences finies sur 65 × 17 points. Cette séparation évite de générer toutes les données par la matrice même que l'on inverse. Des contrôles de raffinement et une solution modale discrète permettent de distinguer erreur de discrétisation et défaut d'implémentation. Les perturbations de hauteur, de +2 % et +5 %, sont testées sans réentraînement ; une classe de traces plus rugueuses est examinée séparément.

## Défendre la validation et la portée des résultats

Le réglage minimise la distance entre la prédiction issue d'une observation bruitée de validation et une **seconde observation indépendante du même objet**. En espérance, le bruit de cette seconde observation ajoute une constante au risque prédictif. La vérité de validation ne sert pas à choisir le paramètre. Cependant, bien prédire les observations ne garantit pas de reconstruire les composantes invisibles : ce critère n'est pas un oracle de l'erreur inverse.

La campagne complète comporte trois entraînements indépendants de 1 000 objets, 200 objets de validation, 200 objets test nominaux et 200 objets rugueux. Les tests nominaux sont réutilisés sous changement de hauteur. Dix réalisations du bruit par objet servent à étudier sa variabilité ; elles ne sont pas dix objets indépendants. Les méthodes partagent objets et bruits pour permettre des comparaisons appariées. Le bootstrap rééchantillonne les objets avec leurs répétitions.

Le [compte rendu d'exécution](../results/full/report.json) enregistre une campagne complète de 2 376 000 reconstructions et zéro échec numérique détecté. Son indicateur `PASS` porte sur l'exécution et ses contrôles ; il ne démontre aucune supériorité du filtre appris. Les conclusions comparatives doivent provenir des métriques, des intervalles et des contrastes préfixés, avec la [configuration effective](../results/full/configuration.json).

Les contrastes principaux comparent le filtre signé plafonné au signé non plafonné à 1 % de bruit : masque contigu nominal et masque complet avec hauteur +5 %. La cible préfixée combine au moins 10 % de réduction du 95e percentile et moins de 5 % d'augmentation de la médiane nominale associée. Les bornes corrigées pour quatre ratios s'interprètent conditionnellement à chaque entraînement ; elles n'offrent pas une couverture simultanée des trois entraînements. Un résultat non concluant ou défavorable reste un résultat à présenter.

## Questions à préparer

**Pourquoi ne pas inverser directement la matrice ?** Les petites singularités amplifient le bruit et l'erreur de modèle. Une petite erreur dans les données peut devenir une grande erreur sur la trace.

**Pourquoi apprendre si Wiener existe déjà ?** Pour examiner l'estimation empirique des relations signal–observation et son comportement lorsque les hypothèses nominales sont imparfaites. Wiener est un comparateur nécessaire, pas un adversaire que l'apprentissage doit obligatoirement battre.

**Une faible erreur de données suffit-elle ?** Non. Les directions peu visibles ou dans le noyau peuvent avoir une grande erreur sans beaucoup changer les données. Il faut regarder l'erreur de la trace entière et ses composantes.

**Le plafond garantit-il la robustesse ?** Il borne une amplification dans des unités définies. Le biais, la loi des profils et l'amplitude effective de l'erreur de modèle restent déterminants. Une borne conditionnelle n'est pas une garantie uniforme sur tous les bruits gaussiens.

**Que peut-on conclure pour des données réelles ?** Seulement que le protocole révèle certains mécanismes dans le modèle synthétique étudié. Géométrie, flux contrôlé, loi des profils et calibration du bruit sont des hypothèses fortes ; une application réelle demanderait une validation spécifique.

**Quel résultat négatif serait utile ?** Un plafonnement sans amélioration, une sélection au bord de la grille, un Wiener aussi bon que l'apprentissage, ou un LMMSE meilleur que les filtres diagonaux. Chacun précise les limites de la méthode et évite d'attribuer à l'apprentissage ce que les données ne démontrent pas.
