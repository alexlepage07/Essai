rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Fonctions_essai.R")


#---- Paramètres ----

# Paramètres des lois marginales.
shape_W <- 0.5
scale_W <- 1 / gamma(1 + 1 / shape_W)

Fw.exact <- function(x) pweibull(x, shape_W, scale_W)
Qw <- function(u) qweibull(u, shape_W, scale_W)

Fx.exact <- function(x) pexp(x, 2)
Qx <- function(u) qexp(u, 2)

# Hyperparamètres de l'approximation`:`
m <- 9     # Longueur des vecteurs de masses de probabilité
h.x <- 1/5  # Pas de discrétisation (h) de la v.a. X
t_ <- 2     # Durée du processus (t)
aa <- 1     # Paramètre de dépendance de la copule fgm

# Hyperparamètres de la simulation:
nb_sim <- 5e+4 # Nombre d'observations simulées pour fins de comparaison.

# Seuils utilisés pour calculer les VaR et TVaR:
kappas <- c(0, 0.5, 0.9, 0.99, 0.999)


#---- Calcul des approximation de densités de probabilités ----

densite_St_lower <- approx_densite_St(
   Fw = Fw.exact,
   Fx = Fx.exact,
   h.x = 1/5,
   t_ = 2,
   m = m,
   aa = 1,
   method = 'lower'
)

densite_St_upper <- approx_densite_St(
   Fw = Fw.exact,
   Fx = Fx.exact,
   h.x = 1/5,
   t_ = 2,
   m = m,
   aa = 1,
   method = 'upper'
)


#---- Simulation ----
set.seed(2020)
St_simul <- simul_Proc_Renouv(
   Qw,
   Qx,
   t_,
   Copule_inverse_fgm,
   aa,
   nb_sim)


#---- Illustration graphique ----

s <- (0:(2 ^ m - 1)) * h.x
cdf_St_simul <- sapply(s, function(j) mean(St_simul <= j))
cdf_St_lower <-  cumsum(densite_St_lower)
cdf_St_upper <-  cumsum(densite_St_upper)

matplot(
   x = s,
   y = cbind(cdf_St_simul, cdf_St_lower, cdf_St_upper),
   type = "s",
   xlab = "x",
   ylab = TeX("$F_{S(t)}(x)$"),
   col = 1:3,
   lty = 1,
   xlim = c(0, min(which(cdf_St_lower > 0.99)) * h.x)
)

legend(
   x = "bottomright",
   legend = c("Simulation",
              "Approximation lower",
              "Approximation upper"),
   col = 1:3,
   lty = 1,
   cex = 0.8,
   box.lty = 1
)


#---- Tableau des mesures de risques ----

print_mesuresRisque_tbl(
   vecteur_support = s, 
   densite_St_lower, 
   densite_St_upper,
   St_simul,
   kappas, 
   print_latex = FALSE
   )
