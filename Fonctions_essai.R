library(latex2exp)
library(xtable)
library(utils)
library(tidyverse)


#---- Calcul des probabilités conjointes de W et X ----

distribution_copule_fgm <- function(Fw, Fx, aa) {
   #' Fonction qui calcule la distribution conjointe des v.a. W et X à l'aide de
   #' la copule de Farlie–Gumbel–Morgenstern (FGM).
   #' 
   #' @param Fw : Vecteur des probabilités cumulées (cdf) de la variable W discrétisée
   #' @param Fx: Vecteur des probabilités cumulées (cdf) de la variable X discrétisée
   #' @param aa : paramètre de dépendance de la copule FGM
   #'
   #' @return Matrice des probabilités cumulées (cdf) conjointes de W et X.
   #' En rangée, on retrouve l'axe associé à la v.a. W
   #' En colonnes, l'axe est associé à la v.a. X
   
   nn <- max(length(Fw), length(Fx))
   
   mat.Fw <- matrix(Fw, nn, nn)
   mat.Fx <- t(matrix(Fx, nn, nn))
   
   mat.Fwx <- mat.Fw * mat.Fx * (1 + aa * (1 - mat.Fw) * (1 - mat.Fx))
   
   return(mat.Fwx)
}

distribution_copule_frank <- function(Fw, Fj, aa) {
   #' Fonction qui calcule la distribution conjointe des v.a. W et J à l'aide de
   #' la copule de Frank.
   #' 
   #' @param Fw : Vecteur des probabilités cumulées (cdf) de la variable W discrétisée
   #' @param Fj: Vecteur des probabilités cumulées (cdf) de la variable J discrétisée
   #' @param aa : paramètre de dépendance de la copule de Frank
   #'
   #' @return Matrice des probabilités cumulées (cdf) conjointes de W et J.
   #' En rangée, on retrouve l'axe associé à la v.a. W
   #' En colonnes, l'axe est associé à la v.a. J
   
   nn <- max(length(Fw), length(Fj))
   
   mat.Fw <- matrix(Fw, nn, nn)
   mat.Fj <- t(matrix(Fj, nn, nn))
   
   if (aa == 0) {
      mat.Fwj <- mat.Fw * mat.Fj
   } else {
      mat.Fwj <-
         (-1 / aa) * log(1 + (exp(-aa * mat.Fw) - 1) * (exp(-aa * mat.Fj) - 1) /
                            (exp(-aa) - 1))
   }
   
   return(mat.Fwj)
}

calcul_probabilites_conjointe <- function(Fwx) {
   #' Fonction qui calcule la matrice des masses de probabilité conjointe
   #' 
   #' @param Fwx : Matrice carrée contenant les probabilités cumulées (cdf) 
   #' conjointes
   #' 
   #' @return Une matrice carrée contenant les masses de probabilité conjointes.
   #' En rangée, on retrouve l'axe associé à la v.a. W
   #' En colonnes, l'axe est associé à la v.a. X
   
   n <- dim(Fwx)[1]
   
   mat.fwx <- matrix(0, n, n)
   
   mat.fwx[1, 1]  <- Fwx[1, 1]
   mat.fwx[1, -1] <- Fwx[1, -1] - Fwx[1, -n]
   mat.fwx[-1, 1] <- Fwx[-1, 1] - Fwx[-n, 1]
   
   mat.fwx[-1, -1] <- 
      Fwx[-1, -1] + Fwx[-n, -n] -
      Fwx[-1, -n] - Fwx[-n, -1]
   
   return(mat.fwx)
}


#---- Approximation de la densité de S(t) ----

approx_densite_St <- function(Fw, Fx, h.x, aa, t_, m=9,
                              method=c("lower", "upper")) {
   #' Fonction qui effectue l'approximation de la densité de probabilité d'un
   #' processus de renouvellement avec récompensenses et dont le lien de 
   #' dépendance unissant les v.a. des temps inter-occurences et des 
   #' "récompenses" est modélisé avec une copule de Farlie–Gumbel–Morgenstern 
   #' (FGM).
   #' 
   #' @param Fw (function): Fonction de répartition de la v.a. des temps 
   #' inter-occurences.
   #' @param Fx (function): Fonction de répartition de la v.a. de la sévérité.
   #' @param m (int): La longueur du vecteur de probabilité qui sera retourné
   #' sera de 2^m. Idéalement, prendre 9 ou 10.
   #' @param aa (float): Paramètre de dépendance de la copule FGM. -1 <= aa <= 1
   #' @param t_ (int): Durée du processus de renouvellement
   #' 
   #' @return Un vecteur de taille (2^m x 1) contenant les masses de 
   #' probabilités associées au processus de renouvellement avec récompenses 
   #' S(t).
   #' 
   method <- match.arg(method)
   
   nn <- 2^m
   
   h.w <- (t_ + 1) / nn
   
   if(method=="lower") {
      
      Fx <- c(0, Fx(1:(nn - 1) * h.x))
      Fw <- Fw(1:nn * h.w)
      
   } else {
      
      Fx <- Fx(1:nn * h.x)
      Fw <- c(0, Fw(1:(nn - 1) * h.w))
      
   }
   
   Fwx <- distribution_copule_fgm(Fw, Fx, aa)
   f_WX <- calcul_probabilites_conjointe(Fwx)
   
   densite_St <- function(){

      k.max <- floor(t_ / h.w + 1)
      pb <- txtProgressBar(min=2, max=k.max, style = 3)
      
      phi_S <- matrix(nrow = k.max, ncol = nn)
      phi_S[1,] <- rep(1, nn)
      
      for (k in 2:k.max) {
         
         phi_S[k,] <- fft(c(1 - Fw[k], rep(0, nn - 1)))
         
         for (j in 2:k) {
            
            phi_S[k, ] <- phi_S[k, ] +
               fft(f_WX[j, ]) * phi_S[k - j + 1, ]
         }
         
         phi_S[k, ] <- phi_S[k, ] + fft(f_WX[1, ]) * phi_S[k, ]
         
         phi_S[k, ] <- nn * phi_S[k, ] / sum(Re(fft(phi_S[k, ], T)))
         
         setTxtProgressBar(pb, k)
      }
      
      
      return(Re(fft(phi_S[k.max, ], inverse = T)) / nn)
   }
   
   return(densite_St())
}


#---- Simulation du processus S(t) ----

Copule_inverse_fgm <- function(aa, v1, u2) {
   # Copule fgm inversée
   
   b <- (1 + aa * (1 - 2 * u2))
   c <- aa * (1 - 2 * u2)
   
   return(2 * v1 / (sqrt(b ^ 2 - 4 * v1 * c) + b))
}

simul_Proc_Renouv <- function(Qw, Qx, t_, Copule_inverse, aa, n_sim=5e+4,
                              return_path=FALSE, seed=FALSE){
   #' Fonction qui permet de simuler des réalisations d'un processus de 
   #' renouvellement avec récompenses.
   #' 
   #' @param Qw (func): Fonction quantile de la v.a. des temps inter-occurences
   #' @param Wx (func): Fonction quantile de la v.a. de l'amplitude des 
   #' "récompenses"
   #' @param t_ (int): Durée du processus
   #' @param Copule_inverse (func) : Fonction associée à l'expression de la 
   #' copule conditionnelle inversée utilisée pour modéliser la dépendance.
   #' @param aa (float): le paramètre de dépendance de la copule
   #' @param n_sim (int) : Nombre de réalisations désirées pour la simulation.
   #' @param return_path (bool): Si TRUE, retourne une matrice contenant les 
   #' temps d'occurence avec la sévérité de chaque occurence.
   #' Si FALSE, retoure une réalisation du processus S(t).
   #' @param seed (bool): Doit-on inclure un ancrage de simulation ?  
   #' 
   #' @return Un vecteur de longueur (n_sim x 1) représentant des réalisations
   #'  simulées du processus S(t).
   
   simul_St <- function() {
      
      if (seed) 
         set.seed(653)
      
      if (n_sim > 2)
         pb <- txtProgressBar(min=1, max=n_sim, style=3)
      
      f <- function(i) {
         Ti <- 0
         Wi <- list()
         U2 <- numeric()
         
         while (Ti <= t_) {
            U2 <- append(U2, runif(1))
            Wi <- append(Wi, Qw(tail(U2, 1)))
            Ti <- as.numeric(tail(Wi, 1)) + Ti
         }
         
         U2 <- as.numeric(U2[-length(U2)])
         Wi <- as.numeric(Wi[-length(Wi)])
         Nt <- length(Wi)
         
         if (Nt == 0) {
            
            St <- 0
            
         } else{
            
            V <- runif(Nt)

            if (aa == 0) {
               U1 <- V
            } else{
               U1 <- mapply(Copule_inverse, rep(aa, Nt), V, U2)
            }
            
            Xi <- Qx(U1)
            
            St <- sum(Xi)
         }
         
         if (n_sim > 2)
            setTxtProgressBar(pb, i)
         
         if (return_path)
            return(cbind("Xi"=Xi, "Ti"=cumsum(Wi)))
        
         return(St)
      }
      
      if (n_sim > 1)
         return(sapply(1:n_sim, function(i) f(i)))
   
      return(f(1))
   }
   
   return(simul_St())
}


#---- Mesures de risques théoriques ----

pmf_theorique <- function(s, h.x, densite_St_approx){
   #' Fonction qui trouve la masse de probabilité associée à une observation s 
   #' selon les approximations. 
   #' 
   #' @param s (float): une observation dont on veut évaluer la probabilité.
   #' @param h.x (float): le pas de discrétisation utilisée pour faire les 
   #' approximations
   #' @param densite_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St.

   return(densite_St_approx[floor(s / h.x) + 1])
}

cdf_theorique <- function(s, h.x, cdf_St_approx){
   #' Fonction qui calcule la probabilité cumulée associée à une observation s 
   #' selon les approximations.
   #' 
   #' @param s (float): une observation dont on veut évaluer la probabilité.
   #' @param h.x (float): le pas de discrétisation utilisée pour faire les 
   #' approximations
   #' @param cdf_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St dont on a réalisé une somme cummulative (cumsum).
   
   return(cdf_St_approx[floor(s / h.x) + 1])
}

VaR_theorique <- function(kappa, vecteur_support, cdf_St_approx){
   #' Fonction qui calcule la valeur au risque (VaR) associée à un seuil kappa  
   #' selon les approximations.
   #' 
   #' @param kappa (float): un seuil tel que 0 <= kappa <= 1
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param cdf_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St dont on a réalisé une somme cummulative (cumsum).
   
   return(min(vecteur_support[cdf_St_approx >= kappa]))
}

stop_loss_theorique <- function(seuil, vecteur_support, h.x, densite_St_approx){
   #' Fonction qui calcule l'espérance d'excédent de seuil (Stop-loss) selon 
   #' les probabilités approximées.
   #' 
   #' @param seuil (float): un seuil supérieur à 0.
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations calculées.
   #' @param cdf_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St dont on a réalisé une somme cummulative (cumsum).
   
   xx <- vecteur_support[vecteur_support > seuil]
   probs <- sapply(xx, function(x) pmf_theorique(x, h.x, densite_St_approx))
   return(sum((xx - seuil) * probs))
}

TVaR_theorique <- function(kappa, vecteur_support, h.x, densite_St_approx,
                           cdf_St_approx=cumsum(densite_St_approx)){
   #' Fonction qui calcule la Tail Value at Risk (TVaR) selon les
   #' probabilités approximées.
   #' 
   #' @param kappa (float): un seuil tel que 0 <= kappa <= 1
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param h.x (float): le pas de discrétisation utilisée pour faire les 
   #' approximations
   #' @param densite_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St.
   #' @param cdf_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St dont on a réalisé une somme cummulative (cumsum).
   
   VaR <- VaR_theorique(kappa, vecteur_support, cdf_St_approx)
   sl <- stop_loss_theorique(VaR, vecteur_support, h.x, densite_St_approx)
   return(sl / (1 - kappa) + VaR)
}

esperance_theorique <- function(vecteur_support, h.x, densite_St_approx) {
   #' Fonction qui calcule l'espérance selon les probabilités approximées.
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param h.x (float): le pas de discrétisation utilisée pour faire les 
   #' approximations
   #' @param densite_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St.
   
   return(sum(
      vecteur_support * pmf_theorique(vecteur_support, h.x, densite_St_approx)
   ))
}

variance_theorique <- function(vecteur_support, h.x, densite_St_approx) {
   #' Fonction qui calcule l'espérance selon les probabilités approximées.
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param h.x (float): le pas de discrétisation utilisée pour faire les 
   #' approximations
   #' @param densite_St_approx (floats): Le vecteur résultant de la fonction
   #' approx_densite_St.
   
   e1 <- esperance_theorique(vecteur_support, h.x, densite_St_approx)
   e2 <- sum(
      vecteur_support^2 * pmf_theorique(vecteur_support, h.x, densite_St_approx)
   )

   return(e2 - e1^2)
}

calcul_mesuresRisque_theoriques <- 
   function(vecteur_support, densite_St_approx, kappas) {
      #' Fonction qui calcule l'espérance, la variance ainsi que les VaR et TVaR
      #' pour différentes valeurs de seuils kappa. Les calculs sont faits à 
      #' partir des probabilités approximées.
      #' 
      #' @param vecteur_support (floats): le vecteur du support sur lequel les 
      #' approximations sont calculées.
      #' @param densite_St_approx (floats): Le vecteur résultant de la fonction
      #' approx_densite_St.
      #' @param kappas (float): un vecteur de seuils kappa tels que 
      #' 0 <= kappa <= 1
      #' 
      #' @return une liste incluant toutes les mesures de risques calculées.
      
      h.x <- diff(vecteur_support[1:2])
      
      cdf_St_approx <- cumsum(densite_St_approx)

      
      mesures_risque <- list()
      
      mesures_risque$esperance <- 
         esperance_theorique(vecteur_support, h.x, densite_St_approx)
      
      mesures_risque$variance <- 
         variance_theorique(vecteur_support, h.x, densite_St_approx)
      
      mesures_risque$VaR <- sapply(kappas, function(k)
            VaR_theorique(k, vecteur_support, cdf_St_approx)
            )
      
      mesures_risque$TVaR <- sapply(kappas, function(k)
         TVaR_theorique(k, vecteur_support, h.x,
                        densite_St_approx, cdf_St_approx)
      )
      
      return(mesures_risque)
   }


#---- Mesures de risques observées -----

VaR_empirique <- function(kappa, vecteur_support, cdf_St_simul){
   #' Fonction qui calcule la Valeur au risque (VaR) empiriquement.
   #'
   #' @param kappa (float): un seuil tel que 0 <= kappa <= 1.
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param cdf_St_simul (floats): Le vecteur des probabilités cumulées 
   #' calculées sur les résultats de simulation avec la fonction ecdf.
   
   return(min(vecteur_support[cdf_St_simul >= kappa]))
}

TVaR_empirique <- function(VaR, St_simul, cdf_St_simul){
   #' Fonction qui calcule la Tail Value at Risk (TVaR) empiriquement.
   #'
   #' @param VaR (float): Valeur au risque
   #' @param St_simul (floats): Le vecteur résultant de la fonction
   #' simul_Proc_Renouv.
   #' @param cdf_St_simul (floats): Le vecteur des probabilités cumulées 
   #' calculées sur les résultats de simulation avec la fonction ecdf.
   
   return(mean(St_simul[St_simul >= VaR]))
}

calcul_mesuresRisque_simul <- function(vecteur_support, St_simul, kappas) {
   #' Fonction qui calcule l'espérance, la variance ainsi que les VaR et TVaR
   #' pour différentes valeurs de seuils kappa. Les calculs sont faits à partir
   #' des réalisations simulées.
   #' 
   #' @param vecteur_support (floats): le vecteur du support sur lequel les 
   #' approximations sont calculées.
   #' @param St_simul (floats): Le vecteur résultant de la fonction
   #' simul_Proc_Renouv
   #' @param kappas (float): un vecteur de seuils kappa tels que 
   #' 0 <= kappa <= 1
   #' 
   #' @return une liste incluant toutes les mesures de risques calculées.
   
   cdf_St_simul <- ecdf(St_simul)(vecteur_support)
   
   
   mesures_risque <- list()
   
   mesures_risque$esperance <- mean(St_simul)
   
   mesures_risque$variance <- var(St_simul)
   
   mesures_risque$VaR <- sapply(kappas, function(k)
      VaR_empirique(k, vecteur_support, cdf_St_simul)
   )
   
   mesures_risque$TVaR <- sapply(mesures_risque$VaR, function(VaR)
      TVaR_empirique(VaR, St_simul, cdf_St_simul)
   )
   
   return(mesures_risque)
}


#---- Produire le tableau des mesures de risque ----

print_mesuresRisque_tbl <- function(vecteur_support, 
                                    densite_St_lower, densite_St_upper,
                                    St_simul, kappas, print_latex=F) {
   #' Fonction qui produit un tableau résumant les mesures de risque pour les
   #' données simulées et avec les probabilités approximées.
   #'
   #' @param vecteur_support (floats): le vecteur du support sur lequel les
   #' approximations sont calculées.
   #' @param densite_St_lower (floats): Le vecteur résultant de la fonction
   #' approx_densite_St pour la méthode lower.
   #' @param densite_St_upper (floats): Le vecteur résultant de la fonction
   #' approx_densite_St pour la méthode upper
   #' @param St_simul (floats): Le vecteur résultant de la fonction
   #' simul_Proc_Renouv
   #' @param kappas (float): un vecteur de seuils kappa tels que
   #' 0 <= kappa <= 1
      
   mesures_lower <- calcul_mesuresRisque_theoriques(
      vecteur_support, densite_St_lower, kappas)
   
   mesures_upper <- calcul_mesuresRisque_theoriques(
      vecteur_support, densite_St_upper, kappas)
      
   mesures_simul <- calcul_mesuresRisque_simul(
      vecteur_support, St_simul, kappas)
   
   
   tbl_moments <- rbind(c(mesures_lower$esperance, mesures_lower$variance),
                        c(mesures_simul$esperance, mesures_simul$variance),
                        c(mesures_upper$esperance, mesures_upper$variance))
   dimnames(tbl_moments) <-  list(c("Lower", "Simulation", "Upper"),
                                  c("Esperance", "Variance"))
   
   
   tbl_VaR <- rbind("Lower" = mesures_lower$VaR,
                    "Simulation" = mesures_simul$VaR,
                    "Upper" = mesures_upper$VaR)
   
   tbl_TVaR <- rbind("Lower" = mesures_lower$TVaR,
                     "Simulation" = mesures_simul$TVaR,
                     "Upper" = mesures_upper$TVaR)
   
   colnames(tbl_VaR) <- kappas
   colnames(tbl_TVaR) <- kappas
   
   
   if (print_latex) {
      
      xtable(cbind(tbl_moments, tbl_TVaR, tbl_VaR), digits = 2)
      
   } else {
      
      print(tbl_moments)
      cat(rep('-', 25), fill=T)
      print(tbl_TVaR)
      cat(rep('-', 25), fill=T)
      print(tbl_VaR)
   }
   
}


#---- Mélanges d'Erlang ----

Approx_MassesProbabilite_Mt <- function(Fj, Fw, t_, distribution_copule, aa) {
   #' Fonction qui approxime les masses de probabilités associées au processus
   #' M(t).
   #' 
   #' @param Fj (float): Vecteur des probabilités cumulées (cdf) de la loi
   #'  de J.
   #' @param Fw (function): Fonction de répartition de la loi de W.
   #' @param t_ (int): Durée du processus de renouvellement M(t).
   #' @param distribution_copule (function): fonction de répartition conjointe
   #'  de W et de J calculée à partir d'une copule.
   #' @param aa (float): Paramètre de dépendance de la copule
   #' 
   #' @return Le vecteur des masses de probabilités associées au processus M(t).
   
   n <- length(Fj)
   h.w <- (t_ + 1) / n # Pas de discrétisation de W
   
   Fw_l <- c(0, Fw((1:(n - 1)) * h.w)) # méthode Lower
   
   Fwj <- distribution_copule(Fw_l, Fj, aa) # Distribution conjointe
   fwj <- calcul_probabilites_conjointe(Fwj) # Masses de probabilités conjointes
   
   
   k.max <- floor(t_/h.w) + 1
   
   pb <- txtProgressBar(min=2, max=k.max, style = 3)
   
   phi_M <- matrix(nrow = k.max, ncol = n)
   phi_M[1,] <- rep(1, n)
   
   for (j in 2:k.max) {
      
      phi_M[j,] <- fft(c(1 - Fw_l[j], rep(0, n - 1)))
      
      for (w in 2:j) {
         
         phi_M[j,] <- phi_M[j,] + fft(fwj[w,]) * phi_M[j - w + 1,]
         
      }
      
      setTxtProgressBar(pb, j)
   }
   
   return(Re(fft(phi_M[k.max, ], inverse = TRUE)) / n)
}

Approx_Distribution_St_Erlang <- function(fm, vecteur_support, scale_Erlang){
   #' Fonction qui approxime la distribution de probabilités associées au 
   #' processus S(t).
   #' 
   #' @param fm (float): Vecteur des masses de probabilité du processus M(t)
   #' calculé avec la fonction Approx_MassesProbabilite_Mt.
   #' @param vecteur_support (float): Vecteur des valeurs représentant le 
   #' support du processus S(t)
   #' @param scale_Erlang (float): Paramètre d'échelle de la loi Erlang
   #' 
   #' @return Le vecteur des probabilités cumulées (cdf) associées au 
   #' processus S(t).
   
   m.max <- which.max(cumsum(fm)) - 1
   
   pb <- txtProgressBar(min = 1, max = m.max, style = 3)
   
   Fs <- rep(fm[1], length(vecteur_support))
   
   for (m in 1:m.max) {
      
      Fs <- Fs + fm[m + 1] * pgamma(vecteur_support, m, scale_Erlang)
      
      setTxtProgressBar(pb, m)
   }
   
   return(Fs)
}

esperance_erlang <- function(fm, scale_Erlang) {
   #' Fonction qui calcule l'espérance du processus S(t).
   #'
   #' @param fm (vecteur): Vecteur des masses de probabilité du processus M(t)
   #' calculé avec la fonction Approx_MassesProbabilite_Mt.
   #' @param scale_Erlang (float): Paramètre d'échelle de la loi Erlang
   
   m <- 0:(which.max(cumsum(fm)) - 1)
   
   esperance_Mt <- sum(m * fm)
   
   return(esperance_Mt / scale_Erlang)
}

print_distribution <- function(Fs, t_, aa, vecteur_support, nb_show=50,
                               digits=4, tbl_latex=F) {
   #' Fonction qui affiche les distributions de probabilités pour chaque 
   #' scénarios testés sous forme de tableaux.
   #' 
   #' @param Fs (float): Matrice des probabilités cumulées (cdf) de S(t) 
   #' approximées avec la fonction Approx_Distribution_St_Erlang, où 
   #' chaque colonne présente une force de dépendance (aa) différente. 
   #' @param t_ (int): Durée du processus.
   #' @param aa (float): vecteur du paramètre de dépendance entre W et J.
   #' @param vecteur_support (vecteur): Vecteur des valeurs représentant le 
   #' support du processus S(t).
   #' @param nb_show (int): Nb d'observations à afficher pour chaque scénario.
   #' @param digits (int): Nombre de décimales à afficher
   #' @param tbl_latex (bool): Est-ce que le tableau doit être affiché en 
   #' format LaTeX ?
   
   nb_data <- length(vecteur_support)
   
   indices_toshow <- round(nb_data / nb_show) * (0:(nb_show - 1)) + 1
   x <- vecteur_support[indices_toshow]
   
   for (k in 1:length(t_)) {
      cat(fill=T)
      cat(rep('-', 15))
      cat(sprintf(" t = %d ", t_[k]))
      cat(rep('-', 15), fill=T)
      cat(fill=T)
      
      tbl_resultats <- Fs[[k]][indices_toshow,]
      tbl_resultats <- round(tbl_resultats, digits)
      
      dimnames(tbl_resultats) <- list(x, aa)
      
      if (tbl_latex) {
         tbl_resultats <- xtable(tbl_resultats, digits=digits)
      }
      
      print(tbl_resultats)
   }
}

print_distribution_plots <- function(Fs, t_, aa, vecteur_support, 
                                     xlim = c(0, 75)) {
   #' Fonction qui affiche les graphiques des distributions de probabilités
   #' du processus S(t) pour tous les niveaux de dépendance et toutes les 
   #' durées du processus désirés.
   #' 
   #' @param Fs (float): Matrice des probabilités cumulées (cdf) de S(t) 
   #' approximées avec la fonction Approx_Distribution_St_Erlang, où 
   #' chaque colonne présente une force de dépendance (aa) différente. 
   #' @param t_ (int): Durée du processus.
   #' @param aa (float): vecteur du paramètre de dépendance entre W et J.
   #' @param vecteur_support (vecteur): Vecteur des valeurs représentant le 
   #' support du processus S(t).
   #' @param xlim (int): limite du support de S(t) pour les fins des 
   #' graphiques.
   
   matplot(
      vecteur_support,
      Fs,
      type = "l",
      ylim = c(0, 1),
      xlim = xlim,
      ylab = TeX(sprintf("$F_{S_{%d}}(x)$", t_)),
      xlab = "x",
      col = 1:5,
      lty = 1:5
   )
   legend(
      "bottomright",
      legend = c(TeX("$\\alpha = $-10 "),
                 TeX("$\\alpha = $-5 "),
                 TeX("Independance"),
                 TeX("$\\alpha = $ 5 "),
                 TeX("$\\alpha = $ 10 ")),
      col = 1:5,
      lty = 1:5,
      cex = 0.8,
      box.lty = 1)
   
}

print_AsymptoticDistribution_plots <- 
   function(Fs, t_, aa, 
            vecteur_support = seq(0, 2, by=0.0001) * t_,
            xlim = c(0, 75)) {
   #' Fonction qui affiche les graphiques des distributions asymptotiques
   #' du processus S(t)/t pour tous les niveaux de dépendance et toutes les 
   #' durées du processus désirés.
   #' 
   #' @param Fs (float): Matrice des probabilités cumulées (cdf) de S(t) 
   #' approximées avec la fonction Approx_Distribution_St_Erlang, où 
   #' chaque colonne présente une force de dépendance (aa) différente. 
   #' @param t_ (int): Durée du processus.
   #' @param aa (float): vecteur du paramètre de dépendance entre W et J.
   #' @param vecteur_support (vecteur): Vecteur des valeurs représentant le 
   #' support du processus S(t).
   #' @param xlim (int): limite du support de S(t) pour les fins des 
   #' graphiques.
   
      matplot(
         vecteur_support / t_,
         Fs,
         type = "l",
         ylim = c(0, 1),
         xlim = c(0, 2),
         ylab = TeX(sprintf("$F_{S_{%d}/%d}(x)$", t_, t_)),
         xlab = "x",
         col = 1:5,
         lty = 1:5
      )
      
      lines(rep(0.4, 11), (0:10) / 10, type = "l", lty = 6, col = 6)
      
      legend(
         "bottomright",
         legend = c(
            TeX("$\\alpha = $-10 "),
            TeX("$\\alpha = $-5 "),
            TeX("Independance"),
            TeX("$\\alpha = $ 5 "),
            TeX("$\\alpha = $ 10 "),
            TeX("$x=E\\[X\\]/E\\[W\\]$")
         ),
         col = 1:6,
         lty = 1:6,
         cex = 0.8,
         box.lty = 1
      )
   }
