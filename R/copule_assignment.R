# ============================================================
# ASSIGNMENT STATISTICA DEI MERCATI FINANZIARI
# Esercizio 1 sulle copule
# Marginali Variance-Gamma + Copula Gaussiana
# ============================================================
#
# NOTA SUI PERCORSI: lo script va eseguito con la cartella di lavoro
# impostata sulla radice del progetto (copula-portfolio-simulation/).
# In RStudio: Session > Set Working Directory > To Source File Location...
# oppure piu' semplicemente apri il file .Rproj nella radice del progetto.

# --- Pacchetti -------------------------------------------------
library(VarianceGamma)   # stima e funzioni della distribuzione Variance-Gamma
library(mvtnorm)         # simulazione dalla normale multivariata
library(writexl)         # esportazione dei risultati in formato Excel
# nota: writexl va installato una sola volta con install.packages("writexl")

# --- Importazione dati ----------------------------------------
# uso il file CSV con i prezzi di chiusura scaricato da Yahoo Finance tramite Python
if (!file.exists("data/prezzi_chiusura.csv")) {
  stop("File data/prezzi_chiusura.csv non trovato. Imposta la working directory sulla radice del progetto.")
}

prezzi <- read.csv("data/prezzi_chiusura.csv", header = TRUE, sep = ",")

# separo le date dalla matrice numerica dei prezzi
date <- as.Date(prezzi$Date)
prezzi_matrice <- as.matrix(prezzi[, -1])

# passo dai prezzi ai log-rendimenti, una colonna per asset
log_ret <- apply(prezzi_matrice, 2, function(x) diff(log(x)))

k <- ncol(log_ret)   # numero di asset
n <- nrow(log_ret)   # numero di osservazioni

# --- Stima delle marginali Variance-Gamma ---------------------
# per ogni asset stimo i 4 parametri (vgC, sigma, theta, nu) via massima
# verosimiglianza.
#
# Accorgimento numerico: i log-rendimenti giornalieri sono numeri molto
# piccoli (ordine 0.01) e su questa scala l'ottimizzatore di vgFit converge
# male. Stimo allora sui rendimenti moltiplicati per un fattore di scala e
# riporto poi indietro i parametri. La Variance-Gamma e' chiusa per
# trasformazioni di scala: se X ~ VG(vgC, sigma, theta, nu) allora
# c*X ~ VG(c*vgC, c*sigma, c*theta, nu), mentre nu resta invariato.
SCALE <- 100

par_vg <- matrix(NA, nrow = k, ncol = 4)
colnames(par_vg) <- c("vgC", "sigma", "theta", "nu")
rownames(par_vg) <- colnames(log_ret)

for (i in 1:k) {
  fit <- vgFit(log_ret[, i] * SCALE)
  p <- fit$param
  # riporto i primi tre parametri alla scala originale; nu e' scale-invariante
  par_vg[i, ] <- c(p[1] / SCALE, p[2] / SCALE, p[3] / SCALE, p[4])
}

par_vg

# --- Diagnostica delle marginali stimate ----------------------
# verifico che le VG stimate riproducano i momenti empirici dei dati:
# media e deviazione standard implicite dai parametri vs quelle campionarie.
# Per la parametrizzazione del pacchetto: E[X] = vgC + theta e
# Var[X] = sigma^2 + nu * theta^2. Il rapporto fra sd stimata e sd empirica
# deve essere vicino a 1 se il fit e' buono.
check_marg <- data.frame(
  asset       = colnames(log_ret),
  media_emp   = colMeans(log_ret),
  media_VG    = par_vg[, "vgC"] + par_vg[, "theta"],
  sd_emp      = apply(log_ret, 2, sd),
  sd_VG       = sqrt(par_vg[, "sigma"]^2 + par_vg[, "nu"] * par_vg[, "theta"]^2),
  row.names   = NULL
)
check_marg$rapporto_sd <- check_marg$sd_VG / check_marg$sd_emp
print(check_marg)

# --- Stima della copula Gaussiana -----------------------------
# trasformo ogni serie di rendimenti nelle pseudo-osservazioni della copula
# applicando la funzione di ripartizione Variance-Gamma stimata (teorema di
# trasformazione integrale): U_i = F_i(X_i) ha distribuzione uniforme.
# Porto poi le uniformi sulla scala normale standard con Phi^-1.
# La matrice di correlazione di queste pseudo-normali e' il parametro Sigma
# della copula Gaussiana.
U <- matrix(NA, nrow = n, ncol = k)
colnames(U) <- colnames(log_ret)

for (i in 1:k) {
  U[, i] <- pvg(log_ret[, i], param = par_vg[i, ])
}

# qui evito che qnorm lavori esattamente su 0 o 1
U <- pmin(pmax(U, 1e-6), 1 - 1e-6)

Z <- qnorm(U)

Sigma <- cor(Z)
colnames(Sigma) <- colnames(log_ret)
rownames(Sigma) <- colnames(log_ret)

Sigma

# --- Simulazione della distribuzione congiunta ----------------
# ricostruisco la distribuzione congiunta stimata seguendo il metodo di
# inversione: dalla normale multivariata alla copula Gaussiana, dalla copula
# alle marginali Variance-Gamma.
set.seed(123)
B <- 10000

# passo 1: simulo dalla normale multivariata con matrice di correlazione Sigma
X <- rmvnorm(B, sigma = Sigma)

# passo 2: trasformo in osservazioni dalla copula Gaussiana
Usim <- pnorm(X)

# passo 3: aggiungo le marginali Variance-Gamma stimate
Rsim <- matrix(NA, nrow = B, ncol = k)
colnames(Rsim) <- colnames(log_ret)

for (i in 1:k) {
  Rsim[, i] <- qvg(Usim[, i], param = par_vg[i, ])
}

# --- Distribuzione del rendimento di portafoglio equi-pesato --
# pesi pari a 1/k su ciascun asset; il rendimento di portafoglio e' la media
# dei rendimenti simulati riga per riga.
pesi <- rep(1 / k, k)
rend_pf <- as.numeric(Rsim %*% pesi)

# sintesi della distribuzione simulata e VaR al 99% e 95%
media_pf <- mean(rend_pf)
sd_pf    <- sd(rend_pf)
VaR      <- quantile(rend_pf, probs = c(0.01, 0.05))

media_pf
sd_pf
VaR

# --- Grafico 1: istogramma della distribuzione simulata -------
png("output/istogramma_portafoglio.png", width = 1000, height = 700)

hist(rend_pf, breaks = 50, probability = TRUE,
     main = "Distribuzione simulata del rendimento di portafoglio equi-pesato",
     xlab = "rendimento giornaliero del portafoglio",
     col = "lightgray", border = "white")

dev.off()

# --- Grafico 2: confronto con la normale ----------------------
# sovrappongo alla densita' simulata una normale con la stessa media e la
# stessa deviazione standard, per mettere in evidenza code piu' spesse e
# asimmetria rispetto al modello gaussiano.
x_norm <- seq(min(rend_pf), max(rend_pf), length.out = 1000)
y_norm <- dnorm(x_norm, mean = media_pf, sd = sd_pf)

png("output/confronto_normale.png", width = 1000, height = 700)

plot(density(rend_pf), col = "blue", lwd = 2,
     main = "Densita' simulata del portafoglio vs Normale",
     xlab = "rendimento giornaliero del portafoglio")

lines(x_norm, y_norm, col = "red", lwd = 2)

legend("topright", legend = c("densita' simulata", "normale di confronto"),
       col = c("blue", "red"), lwd = 2)

dev.off()

# --- Esportazione dei risultati in Excel ----------------------
# quattro fogli: parametri delle marginali, diagnostica del fit, matrice
# della copula, sintesi del portafoglio (media, deviazione standard, VaR).
par_vg_df <- data.frame(asset = rownames(par_vg), par_vg, row.names = NULL)

Sigma_df <- data.frame(asset = rownames(Sigma), Sigma, row.names = NULL)

sintesi_df <- data.frame(
  statistica = c("media", "deviazione standard", "VaR 99%", "VaR 95%"),
  valore     = c(media_pf, sd_pf, VaR[1], VaR[2])
)

write_xlsx(
  list(
    parametri_VG        = par_vg_df,
    diagnostica_fit     = check_marg,
    matrice_copula      = Sigma_df,
    sintesi_portafoglio = sintesi_df
  ),
  path = "output/risultati_copule.xlsx"
)