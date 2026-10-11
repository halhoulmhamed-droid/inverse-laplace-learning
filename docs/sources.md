# Sources et limites de lecture

Bibliographie historique documentée le 10 octobre 2026 ; les notices ci-dessous ne constituent pas une nouvelle vérification des liens ni une lecture intégrale de chaque référence. T : passages techniques ; R : résumé ; M : métadonnées. Aucun texte tiers n’est redistribué.

## Sources du cours et lectures complémentaires

- **C1.** Chapitre « Problèmes mal posés », extrait de cours de 16 pages. Attribution et édition non établies : ne pas attribuer arbitrairement une édition. Texte examiné, figures non inspectées.
- **C2.** Chakir Tajani et Jaafar Abouchabaka, *On the data completion problem for Laplace’s equation*, Annals of the University of Craiova, Mathematics and Computer Science Series 45(1), 11–36, 2018. Texte des 26 pages examiné ; figures non inspectées. Antériorité physique ; la méthode itérative KMF n’est pas reproduite ici.
- **C3.** Maha Farah, *Problèmes inverses de sources et lien avec l’Electro-Encéphalo-Graphie*, UTC, 18 juin 2007. [Archive institutionnelle](https://theses.hal.science/tel-00372423v1). Texte accessible de 169 pages sur 172 examiné ; pages PDF 3, 5 et 11 sans texte retourné, nature visuelle non établie. Pas d’audit intégral des preuves ni des figures.
- **C4.** Azariel Paul Eyimi Minto’o Ebang, *Sur un problème inverse de type Cauchy en théorie des plaques minces élastiques*, Poitiers, 20 janvier 2011. [Archive institutionnelle](https://theses.hal.science/tel-00590516v1). Texte des 177 pages examiné ; pas d’audit intégral des preuves ni des figures.
- **C5.** Anne-Claire Egloffe, *Étude de quelques problèmes inverses pour le système de Stokes. Application aux poumons*, 19 novembre 2012. Support de soutenance de 57 pages, texte examiné ; figures non inspectées. La stabilité sur Stokes avec coefficient Robin ne se transfère pas automatiquement au rectangle scalaire.

Les résultats mathématiques utilisés pour le rectangle sont redémontrés dans `theory.md` et `protocol.md`. Ces lectures contextualisent le travail ; elles ne constituent pas une reproduction des trois thèses ou de résultats biomédicaux.

## Bibliographie externe

#### R01 — Filtrage appris et dépendance aux distributions

**Martin Burger, Samira Kabri.** *Learned Regularization for Inverse Problems: Insights from a Spectral Model.* [arXiv:2312.09845](https://arxiv.org/abs/2312.09845), v1 15/12/2023, v2 04/06/2024. Chapitre publié sous le titre abrégé *Learned regularization for inverse problems*, DOI [10.1515/9783111251233-002](https://doi.org/10.1515/9783111251233-002). Les notices consultées divergent entre 2024 et 2025 pour le volume : citer la **v2 datée** pour les résultats examinés.

**T** : [texte](https://arxiv.org/html/2312.09845v2), §2.1–2.3, §3, formule (9), remarque 2, lemme 1 ; tableau comparatif et description des expériences repérés. Opérateur compact fixé, moments finis, bruit décorrélé du signal. Source centrale sur l’équivalence à une pénalisation spectrale ; aucune nouveauté revendiquée pour cette formule. [Code des auteurs](https://github.com/samirak98/learnedSpectral), README lu : données préenregistrées annoncées à 5,6 Go ; elles ne sont pas nécessaires au minimum A.

#### R02 — Apprendre un régularisateur quadratique est déjà établi

**Giovanni S. Alberti, Ernesto De Vito, Matti Lassas, Luca Ratti, Matteo Santacesaria.** *Learning the optimal Tikhonov regularizer for inverse problems.* **NeurIPS 2021**, vol. 34. [Proceedings officiels](https://papers.nips.cc/paper_files/paper/2021/hash/d3e6cd9f66f2c1d3840ade4161cf7406-Abstract.html) ; [arXiv:2106.06513](https://arxiv.org/abs/2106.06513).

**R** : résumé original et métadonnées. Régularisateur optimal lié aux moments du signal dans le cadre étudié, apprentissage à partir d’un échantillon fini et bornes de généralisation. Sert à classer le filtrage appris comme antériorité, pas à importer une garantie sur notre erreur de modèle. Preuves non lues.

#### R03 — Antériorité spectrale très proche, prépublication 2026

**Hang-Cheng Dong, Pengcheng Cheng, Shuhuan Li.** *Interpretable Operator Learning for Inverse Problems via Adaptive Spectral Filtering: Convergence and Discretization Invariance.* **Prépublication**, 21/03/2026, [arXiv:2603.20602v1](https://arxiv.org/abs/2603.20602v1), DOI [10.48550/arXiv.2603.20602](https://doi.org/10.48550/arXiv.2603.20602).

**T ciblé** : [texte](https://arxiv.org/html/2603.20602v1), §2.3, §3.1–3.2, hypothèse 1 de §4.1 ; résumé des expériences 1D. Filtre neuronal agissant sur les coefficients singuliers. Recouvrement **fort** avec une revendication générique de « filtre spectral neuronal interprétable ». Garanties de convergence et d’invariance non auditées intégralement ; ne pas les transférer au problème de Cauchy, exponentiellement mal posé.

#### R04 — Réseaux pour Cauchy Laplace et biharmonique

**Xin Li, Fajie Wang, Renhao Wang, Sheng-Dong Zhao, Daigeng Yang.** *Fundamental solution neural networks for solving inverse Cauchy problems for the Laplace and biharmonic equations.* **Computers & Mathematics with Applications 201 (2026), 1–17**, DOI [10.1016/j.camwa.2025.09.032](https://doi.org/10.1016/j.camwa.2025.09.032). [Notice éditeur](https://www.sciencedirect.com/science/article/pii/S0898122125004110). [Prépublication des auteurs](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5236948), déposée le 30/04/2025, avec une variante du titre.

**R limité** : métadonnées/résumé indexés de l’éditeur et dépôt des auteurs ; ouvertures directes bloquées par HTTP 403. Méthode à solutions fondamentales. Antériorité suffisante pour exclure « premier réseau pour Cauchy Laplace ». Détails des bruits, coûts et code officiel non vérifiés. Le suffixe 2025 du DOI n’est pas l’année du volume 2026.

#### R05 — RBFNN pour le même problème physique

**F. Mostajeran, S. M. Hosseini.** *Radial basis function neural network (RBFNN) approximation of Cauchy inverse problems of the Laplace equation.* **Computers & Mathematics with Applications 141 (2023), 129–144**, DOI [10.1016/j.camwa.2023.04.026](https://doi.org/10.1016/j.camwa.2023.04.026). [Notice éditeur](https://www.sciencedirect.com/science/article/pii/S0898122123001694).

**R limité** : titre, auteurs, année, pagination et résumé indexés ; accès direct au texte bloqué. Architecture à fonctions radiales ; recouvrement physique certain, supériorité numérique non vérifiée. Ne pas prétendre reproduire cette méthode sans obtenir sa spécification complète. Pas de dépôt officiel confirmé.

#### R06 — Choix des paramètres sans solution vraie

**Michael J. Byrne, Rosemary A. Renaut.** *Learning spectral windowing parameters for regularization using unbiased predictive risk and generalized cross validation techniques for multiple data sets.* **Inverse Problems and Imaging 17(4), 841–869 (2023)**, DOI [10.3934/ipi.2023006](https://doi.org/10.3934/ipi.2023006). [Éditeur](https://www.aimsciences.org/article/doi/10.3934/ipi.2023006), publication en ligne 11/03/2023.

**R avec figures/légendes** : résumé, dates, légendes des figures 3–7. UPRE/GCV multifenêtres, validation sur jeux séparés en défloutage. Motive des sélections fondées sur les observations ; n’atteste pas que GCV ou notre risque prédictif sélectionne toujours la meilleure reconstruction de Cauchy. Démonstrations non examinées.

#### R07 — Stabilité des algorithmes déroulés

**Emilie Chouzenoux, Cécile Della Valle, Jean-Christophe Pesquet.** *Stability Bounds for the Unfolded Forward-Backward Algorithm.* **Journal of Mathematical Imaging and Vision 67, article 48 (2025)**, DOI [10.1007/s10851-025-01258-6](https://doi.org/10.1007/s10851-025-01258-6). [arXiv:2412.17888](https://arxiv.org/abs/2412.17888), dépôt 23/12/2024, v2 01/10/2025.

**T ciblé** : [texte v2](https://arxiv.org/html/2412.17888v2), §1, modèle (1)–(9), définition de la récurrence proximale. Opérateur linéaire connu, pénalité convexe. Motive B et une borne élémentaire par composition ; les bornes fines du §4 ne sont pas reproduites ici. La notice signale un recouvrement textuel avec arXiv:2105.15044 : ne pas compter les deux comme deux contributions indépendantes.

#### R08 — Fondation de LISTA

**Karol Gregor, Yann LeCun.** *Learning Fast Approximations of Sparse Coding.* **ICML 2010**, pp. 399–406. [PDF des proceedings](https://icml.cc/2010/papers/449.pdf).

**T** : PDF pp. 1–5, formulation parcimonieuse, ISTA et construction LISTA, algorithmes 1–4. Dérouler et apprendre une approximation à profondeur fixée est classique. Pour B, garder l’opérateur physique et apprendre seulement quelques paramètres est une adaptation restreinte ; ce n’est pas l’invention du déroulement. Les expériences ICML ne sont pas des expériences EEG.

#### R09 — Antériorité EEG directe

**Ke Liu, Hang Jiang, Hu Yang, Jun Zhang, Zhenghui Gu, Zhuliang Yu, Yu Zhang, Bin Xiao, Wei Wu.** *ADMM-ESINet: A Deep Unrolling Network for EEG Extended Source Imaging.* **IEEE Journal of Biomedical and Health Informatics 29(10), 7262–7273 (octobre 2025)**, DOI [10.1109/JBHI.2025.3568648](https://doi.org/10.1109/JBHI.2025.3568648). [Notice bibliographique et résumé original](https://pubmed.ncbi.nlm.nih.gov/40343816/).

**R + README** : publication et auteurs vérifiés dans la notice ; texte IEEE non accessible. [Dépôt associé](https://github.com/hangj-cache/ADMM-ESINet) lu : génération annoncée via Brainstorm/MATLAB, sources à haute résolution et liens de données. Exécution, licence de réutilisation et coût non audités. Notre B serait une étude réduite, **pas une reproduction complète** d’ADMM-ESINet ni une validation clinique.

#### R10 — Robin, stabilité et approche bayésienne

**Aksel K. Rasmussen, Fanny Seizilles, Mark Girolami, Ieva Kazlauskaite.** *The Bayesian Approach to Inverse Robin Problems.* **SIAM/ASA Journal on Uncertainty Quantification 12(3), 1050–1084 (2024)**, DOI [10.1137/23M1620624](https://doi.org/10.1137/23M1620624). [Archive institutionnelle DTU](https://orbit.dtu.dk/en/publications/the-bayesian-approach-to-inverse-robin-problems/) ; [arXiv:2311.17542](https://arxiv.org/abs/2311.17542), 29/11/2023.

**R et formulation introductive** : résumé, dates et premiers modèles Laplace/Stokes du [texte](https://arxiv.org/html/2311.17542v1). Coefficient caché, données de Cauchy accessibles ; taux dépendant de classes de régularité. Les résultats de consistance postérieure ne sont pas transférés au substitut GP de C. Les preuves et la chaîne MCMC ne sont pas auditées.

#### R11 — Robin : référence numérique récente solide

**Erik Burman, Siyu Cen, Bangti Jin, Zhi Zhou.** *Numerical approximation and analysis of the inverse Robin problem using the Kohn-Vogelius method.* **Mathematics of Computation, accepté/en cours de publication selon UCL (2026)**, DOI [10.1090/mcom/4222](https://doi.org/10.1090/mcom/4222). [UCL et manuscrit accepté](https://discovery.ucl.ac.uk/id/eprint/10229703/) ; [prépublication du 09/06/2025](https://arxiv.org/abs/2506.07370).

**T ciblé** : résumé, passages du §2 autour de l’hypothèse 2.1 et du théorème 2.1, PDF pp. 6–9. Coefficient par morceaux, approximation éléments finis et hypothèses sur les données bruitées ; analyse d’erreur conditionnelle. Ne pas revendiquer une première analyse numérique de Robin. Pagination du volume final non confirmée ; version acceptée explicitement identifiée.

#### R12 — Origine de KMF

**V. A. Kozlov, V. G. Maz’ya, A. V. Fomin.** *An iterative method for solving the Cauchy problem for elliptic equations.* **1991**, Zh. Vychisl. Mat. Mat. Fiz. 31(1), 64–74 ; traduction U.S.S.R. Computational Mathematics and Mathematical Physics 31(1), 45–52. [Archive Math-Net](https://www.mathnet.ru/php/archive.phtml?jrnid=zvmmf&option_lang=eng&paperid=3145&wshow=paper).

**M** : référence originale confirmée ; texte original non lu. Le schéma effectivement examiné est celui de C2, §4.4.2–4.4.4. DOI non vérifié, donc non ajouté. Référence de provenance, sans attribuer à cette mission la lecture de la démonstration de 1991.

#### R13 — Référence ML du substitut de C

**Carl Edward Rasmussen, Christopher K. I. Williams.** *Gaussian Processes for Machine Learning.* **MIT Press, 2006**, ISBN 0-262-18253-X. [Site des auteurs](https://gaussianprocess.org/gpml/) ; [chapitre 2, Regression](https://gaussianprocess.org/gpml/chapters/RW2.pdf).

**T ciblé** : §2.1 et passages autour de l’algorithme 2.1, p. imprimée 19 / PDF du chapitre p. 13. Régression GP, factorisation de Cholesky et coût cubique en nombre de points d’apprentissage. Référence méthodologique fondatrice pour un substitut du modèle direct, pas étude du Robin proposé. La variance prédictive du GP n’est pas une incertitude certifiée sur le coefficient inverse.

#### R14 — Prolongement récent de l’axe de l’encadrant

**Jamal Daoudi, Chakir Tajani.** *Genetic Algorithm-Based Optimization Approach for Solving a Class of Inverse Problems with Tikhonov Regularization.* **WSEAS Transactions on Mathematics 22, 842–853 (2023)**, DOI [10.37394/23206.2023.22.92](https://doi.org/10.37394/23206.2023.22.92). [Éditeur](https://www.wseas.com/journals/articles.php?id=8592).

**R** : notice, résumé et dates. Complétion Laplace formulée avec Tikhonov, optimisation génétique. Confirme un lien thématique récent avec l’encadrement. Un algorithme génétique n’est pas assimilé à un apprentissage statistique. Pour notre objectif quadratique convexe, une résolution linéaire constitue la référence plus directe ; aucune comparaison de performances avec cet article n’est prétendue.
