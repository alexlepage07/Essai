# Mon essai de maîtrise

## Distribution d’un accroissement d’un processus de renouvellement avec récompense dans un contexte de dépendance

Sous la direction de:
- Hélène Cossette, directrice de recherche
- Étienne Marceau, codirecteur de recherche

2021


### Résumé:
Les processus de renouvellement avec récompenses sont très utiles notamment, en théorie
de la abilité (p. ex. Pal and Murthy (2003)) et dans la théorie de la ruine. Récemment,
dans la littérature, plusieurs auteurs se sont intéressé à relaxer l'hypothèse d'indépendance
qui existe généralement entre les variables des temps inter-sinistres et de la sévérité. Ces
chercheurs ont développé diérentes mesures de risque an d'étudier le comportement d'un
tel modèle. Mentionnons les deux premiers moments, la transformée de Laplace-Stieltjes, la
fonction de pénalité de Gerber-Shiu, etc. Cependant, aucune méthodologie n'a encore été
développée pour approximer la distribution d'un accroissement de ce processus aléatoire hormis
Marceau (2009) qui s'est intéressé au cas où les variables aléatoires sont discrètes. Le présent
essai vise donc à adapter les travaux de ce dernier lorsque les variables aléatoires sont continues.
An de modéliser la dépendance qui existe entre les temps inter-sinistres et les récompenses, la
théorie des copules est utilisée. En complément à cela, un algorithme est proposé pour eectuer
de la simulation Monte-Carlo. De plus, on regarde les propriétés asymptotiques en contexte
de dépendance et la démonstration de l'ordonnancement de ce processus selon le niveau de
dépendance est eectuée.


### Abstract:
Renewal processes with rewards are particularly useful in reliability theory (eg Pal and Murthy
(2003)) and in ruin theory. Recently, in the literature, several authors have been interested in
relaxing the independence hypothesis which generally exists between inter-disaster time and
severity variables. These researchers have developed dierent risk measures in order to study
the behavior of such a model. Let us mention the rst two moments, the Laplace-Stieltjes
transform, the Gerber-Shiu penalty function, etc. However, no methodology has yet been developed
to approximate the distribution of an increase in this random process except Marceau
(2009) which was interested in the case where the random variables are discretes. This essay
therefore aims to adapt the latter's work when the random variables are continuous. In order
to model the dependence that exists between inter-disaster times and rewards, copula theory
is used. In addition to this, an algorithm is proposed to perform Monte-Carlo simulation.
Moreover, we look at the asymptotic properties in a dependency context and the ordering of
this process according to the dependency level is demonstrated.


## Contenu du répertoire Git:
- '111144776.pdf': l'essai en format pdf
- 'Approximation_distribution_St.R': le document R où la distribution d'un accroissement d'un processus de renouvellement avec récompense est approximé conformément à la Section 2.2.3 de l'essai.
- 'Fonctions_essai.R': L'essentiel des fonctions utilisées et l'appel des paquetages se fait ici.
- 'Melange_Erlang.R': Ici, on fait l'approximation de distribution en utilisant un mélange d'Erlang pour produire les résultats de la Section 2.3 de l'essai.
- 'Simulation_Parcours_St.R': Ici, on simule le *parcours* d'un processus de renouvellement avec récompense.
- 'Simulation_processus_renouvellement_recompense.R': Ici, on présente l'algorithme de simulation présenté à la Section 2.1 de l'essai.
