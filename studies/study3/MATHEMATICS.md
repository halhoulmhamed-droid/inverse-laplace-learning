# Construction et audit mathématique — étude 3

11 octobre 2026. Écrit avant implémentation et test final ; revue indépendante demandée avant campagne.

## Coordonnées, centrage, échelle et bruit
Le signal est le vecteur c des 64 premiers coefficients cosinus orthonormaux. La vérité synthétique comporte 128 modes indépendants ; les modes omis et la différence entre générateur analytique et direct FD sont une erreur de modèle, pas un bruit que l'on prétend connaître. A est l'opérateur FD nominal (64×16 intervalles, H=0,25), restreint au masque puis divisé par 8. L'observation est y=M(B_H c_128 + σ ε)/8, ε~N(0,I_64), indépendante des coefficients. Le bruit connu après masquage a covariance νI_m, ν=σ²/64. σ est absolu et n'est pas réestimé sur le test.

À partir de N signaux d'apprentissage, μ=(1/N)Σc_i et Cemp=(1/(N−1))Σ(c_i−μ)(c_i−μ)ᵀ. Les vérités synthétiques d'apprentissage sont autorisées ; aucune vérité de test n'est utilisée. La covariance concerne le signal, non les observations ni les blocs croisés de l'étude 2.

Poser Q=diag(1+(πk)²), k=0..63 ; la constante k=0 a une pénalité 1. Si tr(Cemp)>0, s=tr(Cemp)/tr(Q⁻¹)>0 et CH1=sQ⁻¹. Ainsi CH1 et Cemp ont la même trace, estimée uniquement sur l'apprentissage. CH1 impose une décroissance H1 dans ces coordonnées, et n'est pas la vraie covariance nominale, dont les variances décroissent comme (1+k²)^(-2,5). Pas d'universalité de cette cible. Si la trace est nulle ou non finie, l'ajustement est invalide et déclaré, sans remplacement caché.

Cγ=(1−γ)Cemp+γCH1, 0≤γ≤1. Le modèle fondé sur la covariance utilise
Rγ,α=Cγ Aᵀ(A Cγ Aᵀ+βI_m)⁻¹, β=α ν>0,
ĉ=μ+Rγ,α(y−Aμ).
On ajuste α sur validation : β est la pénalisation effective, pas une affirmation que le vrai bruit vaut α fois le bruit physique.

## Positivité, solvabilité et forme variationnelle
Cemp est symétrique PSD par construction, même quand N−1<64. CH1 est SPD si s>0 ; Cγ est PSD sur [0,1] et SPD pour γ>0. Pour tout v≠0,
vᵀ(A Cγ Aᵀ+βI)v=||Cγ^(1/2)Aᵀv||²+β||v||²≥β||v||²>0.
Le système d'observation est donc SPD, y compris pour γ=0 et une covariance de rang déficient. Résolution par Cholesky, pas d'inversion explicite ni de « jitter » non déclaré. Toute panne numérique est conservée. Une très faible valeur propre négative à la précision machine de Cemp n'est pas interprétée comme une covariance physique négative ; les résidus et tolérances des tests sont enregistrés.

Pour Cγ SPD, l'identité de push-through donne Rγ,α=(AᵀA+β Cγ⁻¹)⁻¹Aᵀ. Il s'agit de l'unique minimiseur de ||A x−(y−Aμ)||²+β xᵀCγ⁻¹x, avec ĉ=μ+x. Pour Cemp singulière, cette écriture avec inverse ordinaire est interdite ; la formule dans l'espace d'observation reste valide et produit un incrément dans Range(Cemp). Elle est aussi obtenue en paramétrant x=Cemp^(1/2)z et pénalisant ||z||².

## Correspondances exactes et familles distinctes
- γ=0 est exactement M0 : plug-in de Cemp dans le modèle nominal, même μ, même A, même α. Ce n'est pas en général le LMMSE empirique de l'étude 2.
- γ=1 donne R=(AᵀA+(β/s)Q)⁻¹Aᵀ. Ce H1 apparié (H1m) est centré sur Aμ. Dans la convention de l'étude 2 Qn=Q/Qmax, le paramètre équivalent est a=β Qmax/s, et non α lui-même. Le H1 historique peut avoir une autre grille et est centré sur μy empirique ; ses reconstructions ne sont donc pas nécessairement égales.
- L'empirique original (E) conserve exactement μ+ Ccy(Cyy+λ ty I)⁻¹(y−μy), ty=tr(Cyy)/m, Ccy et Cyy calculées sur les couples bruités. Les différences de blocs et de centre comptent : une covariance signal avec A nominal n'est pas égale à une covariance jointe empirique bruitée. Une égalité requiert notamment μy=Aμ, Ccy=Cemp Aᵀ et Cyy+λtyI=A Cemp Aᵀ+βI, conditions non supposées ici.
- Diag remplace Cemp par diag(diag(Cemp)) **dans la base cosinus**. Ce n'est pas D de l'étude 2, qui diagonalise Cbz et Czz dans la base singulière retenue.
- H1 historique (H1) conserve la formule affine de l'étude 2, μ+(AᵀA+aQn)⁻¹Aᵀ(y−μy). Les composants existants restent inchangés.

## Attribution de l'effet de mélange
La méthode S choisit γ et α ; M0 choisit α sur la même grille. La différence primaire compare des procédures réglées : S dispose de 140 candidats contre 20 pour M0, différence annoncée et non présentée comme budget identique. Deux diagnostics sans nouveau réglage utilisent M0 et H1m au **même α sélectionné pour S**. S contre M0 au α commun isole strictement le changement de covariance, à μ, A, bruit et β fixes. La comparaison H1m/H1 historique ne doit pas être attribuée au seul changement de covariance (centrage et grilles différents).

Il n'existe pas de garantie que γ>0 améliore le risque test : le mélange réduit certaines erreurs d'estimation mais peut introduire un biais de cible, notamment sous changement de distribution. La PSD et la solvabilité ne sont pas des preuves de robustesse statistique.

## Antériorité
La contraction de covariance vers une cible structurée est classique (Ledoit–Wolf), et l'apprentissage de régularisateurs quadratiques fondés sur les moments du signal est établi (Alberti et al.). Cette expérience étudie une cible H1 et un protocole de reconstruction particuliers ; elle ne reproduit pas la règle analytique de Ledoit–Wolf et ne revendique pas d'invention.
