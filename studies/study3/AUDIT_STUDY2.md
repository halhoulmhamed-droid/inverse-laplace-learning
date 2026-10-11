# Audit préalable des archives de l’étude 2

Référence : `1f11fea35d5367a355901c3866d21f5b3b21a5b0`, vérifiée localement et sur GitHub le 11 octobre 2026. Audit indépendant effectué avant le gel de l’étude 3, sans relance de la campagne antérieure.

Aucune anomalie substantielle détectée dans les éléments ciblés. Trois empreintes du protocole figé, cinq empreintes de sources et 150 empreintes d’archives concordent. Le gel `20a663d` précède le commit scientifique `57ab54e`. Les partitions PCG64/SeedSequence sont disjointes ; les tailles imbriquées et objets de test communs aux cinq entraînements sont intentionnels.

Recompositions à partir des archives : 7 800 statistiques (écart nul), 40 000 valeurs bootstrap (écart maximal 5,56e−17), 105 modèles projetés (écart maximal 5,301e−12), 468 pertes scalaires avec un transfert indépendant par intégrales de sinus (écart maximal 2,013e−15). Les 800 objets vrais sauvegardés sont retrouvés à l’arrondi numérique près. Aucune panne ni entrée non finie.

Le contraste F/D de l’étude 2 satisfait ses critères dans les cinq jeux. Il porte sur la diagonalisation de Cbz et Czz en base singulière r=12, et non sur Ccc cosinus. Les identités projetées sont valides à lambda commun et moyenne identique. Les résultats sont conditionnels et limités au protocole ; les comparaisons à H1 y sont descriptives.

L’estimateur empirique Ccy(Cyy+λtyI)⁻¹ avec centre μy n’est pas automatiquement le plug-in Cemp Aᵀ(A Cemp Aᵀ+βI)⁻¹ avec centre Aμ. Cette distinction est incorporée dans les formules et comparateurs de l’étude 3. Les coefficients cosinus indépendants n’établissent aucune corrélation physique nouvelle.

Les scripts de recomposition et journaux détaillés de cet audit restent privés. Une assertion initiale d’égalité bit à bit a été remplacée par une tolérance adaptée à l’arrondi des puissances scalaires/vectorisées (2,78e−17), sans modification des archives ni incidence scientifique. L’audit n’est pas une reproduction exhaustive de chaque étape de l’étude 2.
