rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Fonctions_essai.R")


#---- Paramètres ----

shape_W <- 1.1
scale_W <- 1 / gamma(1 + 1/shape_W)

Qw <- function(u) qweibull(u, shape_W, scale_W)

Qx <- function(u) qexp(u, 1)

t_ <- 5                  # Durée du processus
aa <- c(-1, 0, 1)        # Paramètre de dépendance de la copule fgm
kappa <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99) # Seuils à considérer pour les mesures
                         # de risque


#---- Simulation ----

set.seed(2020)
# S_ <- sapply(aa, function(theta) simul_Proc_Renouv(theta))

S_ <- sapply(aa, function(theta)
    simul_Proc_Renouv(
        Qw,
        Qx,
        t_,
        Copule_inverse = Copule_inverse_fgm,
        aa = theta,
        n_sim = 1e+5
    ))


#---- Illustration graphique ----
plot(c(0, 20), c(0, 1),
    type = "n",
    ylab = TeX("$F_{S(t)}(x)$"),
    xlab = TeX("$x$")
)
sapply(1:3, function(j)
    lines(ecdf(S_[, j]), col = j))

legend("bottomright", 
       legend = c(TeX("$\\theta$ =-1 "),
                  TeX("$\\theta$ = 0 "),
                  TeX("$\\theta$ = 1 ")),
       col = 1:3,
       lty = 1, cex =1.2,
       box.lty=0)


#---- Calcul des mesures de risques VaR et TVaR ----
cdf_simul <- sapply(1:3,function(j)
    ecdf(S_[,j])(S_[,j]))

VaR <- sapply(1:3, function(j)
    sapply(kappa, function(k)
        VaR_empirique(k, S_[, j], cdf_simul[, j])))

TVaR <- sapply(1:3, function(j)
    sapply(VaR[,j], function(k)
        TVaR_empirique(k, S_[, j], cdf_simul[, j])))


xtable(rbind(kappa, "VaR"=t(VaR), "TVaR"=t(TVaR)),
       digits = 2, include.rownames=F)
