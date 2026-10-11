# Précision mathématique après campagne — 11 octobre 2026

La revue indépendante finale relève une portée logique trop forte dans la phrase « Une égalité requiert notamment » du document `MATHEMATICS.md` figé. Les trois identités de moments/centre qui y suivent sont des **conditions suffisantes**, et non nécessaires en général, pour égaler les deux estimateurs affines M0 et E à paramètres correspondants. Des compensations dans les produits matriciels, ou un écart de centre dans le noyau de l’opérateur de reconstruction, peuvent également produire le même estimateur.

Pour écrire la condition nécessaire et suffisante au niveau des applications affines (même moyenne signal μ), il suffit de vérifier simultanément `R_M0 = R_E` et `R_M0 Aμ = R_E μ_y`. Les trois identités de moments constituent une manière suffisante de l’obtenir, sans être la seule.

Cette clarification a été écrite **après** la campagne et avant publication. Le document gelé, son empreinte et son commit restent intacts. Les familles implémentées, les grilles, les données, le protocole, les extrémités γ=0/1, les contrôles numériques et la conclusion ne changent pas. Le rapport adopte la formulation corrigée. Aucune expérience n’a été relancée pour cette correction rédactionnelle.
