rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Fonctions_essai.R")


#---- Paramètres ----

# Paramètres de la loi de W
shape_W = 1.5
scale_W = 5 / gamma(1 + 1 / shape_W)

Fw <- function(x) pweibull(x, shape_W, scale_W)

# Paramètres de la loi de J
m <- 10             # La longueur du vecteur de support de J est de 2^m
j <- 0:(2 ^ m - 1)  # Vecteur de support de la v.a. J.

prob_J <- 1/2              # Paramètre de la loi de J.
Fj <- 1 - (1 - prob_J) ^ j # On pose que J suit une loi Géométrique.

# Autres paramètres
scale_Erlang <- 1          # Paramètre d'échelle de la loi Erlang
aa <- c(-10, -5, 0, 5, 10) # Paramètre de dépendance de la copule de Frank
t_ <- c(1:5 * 10, 1000)    # La durée du processus S(t)
vecteur_support <- seq(0, 2^m, by=0.01) # Support sur lequel on désire calculer 
                           # les masses de probabilité du processus S(t).

#---- Calcul des distributions de probabilités ----

# Calcul des masses de probabilités du processus M(t)
fm <-
    lapply(t_, function(k)
        sapply(aa, function(a)
            Approx_MassesProbabilite_Mt(
                Fj,
                Fw,
                t_ = k,
                distribution_copule_frank,
                aa = a
            )))

# Distribution du processus S(t)
Fs <-
    lapply(1:length(t_), function(k)
        sapply(1:length(aa), function(i)
            Approx_Distribution_St_Erlang(
                fm[[k]][, i],
                vecteur_support,
                scale_Erlang)))

# Distribution asymptotique de S(t)/t
Fs_asympt <-
    lapply(1:length(t_), function(k)
        sapply(1:length(aa), function(i)
            Approx_Distribution_St_Erlang(
                fm[[k]][, i],
                vecteur_support = seq(0, 2, by = 0.0001) * t_[k],
                scale_Erlang)))


#---- Affichage des distributions de probabilités ----
print_distribution(
    Fs,
    t_,
    aa,
    vecteur_support,
    nb_show = 25,
    digits = 4,
    tbl_latex = FALSE
)


#---- Affichage des graphiques de distribution ----

# Distribution du processus S(t)
for (k in 1:length(t_)) {
    print_distribution_plots(Fs[[k]], t_[k], aa, vecteur_support)
}

# Distribution asymptotique de S(t)/t
for (k in 1:length(t_)) {
    print_AsymptoticDistribution_plots(Fs_asympt[[k]], t_[k], aa)
}


#---- Calcul de l'espérance ----

nn <- 2^m
h.w <- round((t_ + 1) / nn, 4)

tbl_esperances <-
    sapply(1:length(t_), function(k)
        sapply(1:length(aa), function(i)
            esperance_erlang(fm[[k]][, i], scale_Erlang)
            )
        ) %>% t()

tbl_esperances <- cbind(h.w, tbl_esperances)

dimnames(tbl_esperances) <- list(t_, c("h.w", aa))


tbl_esperances_asympt <- cbind(h.w, tbl_esperances[,-1] / t_)


xtable(tbl_esperances, digits=3)
xtable(tbl_esperances_asympt, digits=4)
