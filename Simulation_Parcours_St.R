rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Fonctions_essai.R")


#---- Paramètres ----

shape_W <- 1.1
scale_W <- 1 / gamma(1 + 1/shape_W)

Qw <- function(u) qweibull(u, shape_W, scale_W)

Qx <- function(u) qexp(u, 1)

t_ <- 10              # Durée du processus
aa <- c(-1,0,1)       # Paramètre de dépendance de la copule fgm


#---- Simulation ----

XT <- lapply(aa, function(theta)
    simul_Proc_Renouv(
        Qw,
        Qx,
        t_,
        Copule_inverse = Copule_inverse_fgm,
        aa = theta,
        n_sim = 1,
        return_path = TRUE,
        seed=TRUE)
    )


S <- lapply(1:3, function(i) {
    XT[[i]] %>% as_tibble() %>% transmute(Si = cumsum(Xi), Ti)
})
S[[1]] <- rbind(c(0,0), S[[1]], c(as.numeric(S[[1]][nrow(S[[1]]), 1]), t_))
S[[2]] <- rbind(c(0,0), S[[2]], c(as.numeric(S[[2]][nrow(S[[2]]), 1]), t_))
S[[3]] <- rbind(c(0,0), S[[3]], c(as.numeric(S[[3]][nrow(S[[3]]), 1]), t_))


ggplot() + geom_step(
    aes(x = Ti, y = Si, colour='red'),
    data = as.data.frame(S[[1]])
) + geom_step(
    aes(x = Ti, y = Si, colour='green'),
    data = as.data.frame(S[[2]])
) + geom_step(
    aes(x = Ti, y = Si, colour='blue'),
    data = as.data.frame(S[[3]])
) + scale_color_discrete(
    name=TeX("$\\theta$"),
    labels=aa
) + xlab(TeX("$T_i$")
) + ylab(TeX("$S_i$"))
