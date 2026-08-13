# Assignment Statistica dei mercati finanziari - Punto 2

library(VarianceGamma)
library(mvtnorm)
library(writexl)


# Importo i prezzi di chiusura giornalieri (scaricati da Yahoo Finance con uno script Python che ho fatto a parte)

prezzi <- read.csv("data/prezzi_chiusura.csv", header = TRUE, sep = ",")

date <- as.Date(prezzi$Date)
prezzi_chiusura <- as.matrix(prezzi[, -1])

# Calcolo i log-rendimenti giornalieri

rendimenti <- apply(prezzi_chiusura, 2, function(x) diff(log(x)))

k <- ncol(rendimenti)
n <- nrow(rendimenti)


# Stimo per ogni asset i 4 parametri della Variance-Gamma con massima verosimiglianza

parametri_VG <- matrix(NA, nrow = k, ncol = 4)
colnames(parametri_VG) <- c("c", "sigma", "theta", "nu")
rownames(parametri_VG) <- colnames(rendimenti)

for (i in 1:k) {
  fit <- vgFit(rendimenti[, i] * 100)
  p <- fit$param
  parametri_VG[i, ] <- c(p[1]/100, p[2]/100, p[3]/100, p[4])
}

parametri_VG


# Controllo se le VG stimate riproducono media e sd dei dati.

controllo_fit <- data.frame(
  asset = colnames(rendimenti),
  media_emp = colMeans(rendimenti),
  media_VG = parametri_VG[, "c"] + parametri_VG[, "theta"],
  sd_emp = apply(rendimenti, 2, sd),
  sd_VG = sqrt(parametri_VG[, "sigma"]^2 + parametri_VG[, "nu"] * parametri_VG[, "theta"]^2),
  row.names = NULL
)
controllo_fit$rapporto_sd <- controllo_fit$sd_VG / controllo_fit$sd_emp
print(controllo_fit)


# Stimo la copula gaussiana (teorema di trasformazione integrale)

uniformi <- matrix(NA, nrow = n, ncol = k)
colnames(uniformi) <- colnames(rendimenti)

for (i in 1:k) {
  uniformi[, i] <- pvg(rendimenti[, i], param = parametri_VG[i, ])
}

# Schiaccio gli estremi 0/1 per evitare che qnorm restituisca +/- Inf.

uniformi <- pmin(pmax(uniformi, 1e-6), 1 - 1e-6)

normali_standard <- qnorm(uniformi)

Sigma <- cor(normali_standard)
colnames(Sigma) <- colnames(rendimenti)
rownames(Sigma) <- colnames(rendimenti)

Sigma


# Simulo dalla distribuzione congiunta col metodo di inversione:
# (1) normale multivariata con correlazione Sigma,
# (2) la passo nella Phi per ottenere uniformi dalla copula,
# (3) le inverto con qvg per applicare le marginali VG stimate.

set.seed(123)
B <- 10000

normali_simulate <- rmvnorm(B, sigma = Sigma)
uniformi_simulate <- pnorm(normali_simulate)

rendimenti_simulati <- matrix(NA, nrow = B, ncol = k)
colnames(rendimenti_simulati) <- colnames(rendimenti)

for (i in 1:k) {
  rendimenti_simulati[, i] <- qvg(uniformi_simulate[, i], param = parametri_VG[i, ])
}


# Portafoglio equi-pesato.
# I rendimenti simulati sono log-rendimenti; li riconverto prima in
# rendimenti semplici, che si aggregano linearmente con i pesi di portafoglio.

pesi <- rep(1/k, k)
rendimenti_semplici_simulati <- exp(rendimenti_simulati) - 1
rendimento_portafoglio <- as.numeric(rendimenti_semplici_simulati %*% pesi)

media_portafoglio <- mean(rendimento_portafoglio)
sd_portafoglio <- sd(rendimento_portafoglio)
VaR <- quantile(rendimento_portafoglio, probs = c(0.01, 0.05))

# Skewness e kurtosis: una normale ha skewness=0 e kurtosis=3, se i miei
# valori si discostano da quelli il portafoglio non e' gaussiano.

scarti_standardizzati <- (rendimento_portafoglio - media_portafoglio) / sd_portafoglio
skewness_portafoglio <- mean(scarti_standardizzati^3)
kurtosis_portafoglio <- mean(scarti_standardizzati^4)

media_portafoglio
sd_portafoglio
VaR
skewness_portafoglio
kurtosis_portafoglio


# Creo la cartella di output se non esiste

if (!dir.exists("results/copula")) dir.create("results/copula", recursive = TRUE)


# Istogramma della distribuzione simulata del portafoglio

png("results/copula/istogramma_portafoglio.png", width = 1000, height = 700)
hist(rendimento_portafoglio, breaks = 50, probability = TRUE,
     main = "Distribuzione simulata del rendimento di portafoglio equi-pesato",
     xlab = "rendimento giornaliero del portafoglio",
     col = "lightgray", border = "white")
dev.off()


# Confronto con una normale di stessa media e sd, per vedere quanto sono  piu' grasse le code rispetto al modello gaussiano

griglia_x <- seq(min(rendimento_portafoglio), max(rendimento_portafoglio), length.out = 1000)
densita_normale <- dnorm(griglia_x, mean = media_portafoglio, sd = sd_portafoglio)

png("results/copula/confronto_normale.png", width = 1000, height = 700)
plot(density(rendimento_portafoglio), col = "blue", lwd = 2,
     main = "Densita' simulata del portafoglio vs Normale",
     xlab = "rendimento giornaliero del portafoglio")
lines(griglia_x, densita_normale, col = "red", lwd = 2)
legend("topright", legend = c("densita' simulata", "normale di confronto"),
       col = c("blue", "red"), lwd = 2)
dev.off()


# Esporto tutto in un Excel con quattro fogli.

parametri_VG_df <- data.frame(asset = rownames(parametri_VG), parametri_VG, row.names = NULL)
Sigma_df <- data.frame(asset = rownames(Sigma), Sigma, row.names = NULL)

sintesi_df <- data.frame(
  statistica = c("media", "deviazione standard", "VaR 99%", "VaR 95%", "skewness", "kurtosis"),
  valore = c(media_portafoglio, sd_portafoglio, VaR[1], VaR[2], skewness_portafoglio, kurtosis_portafoglio)
)

write_xlsx(
  list(
    parametri_VG = parametri_VG_df,
    controllo_fit = controllo_fit,
    matrice_copula = Sigma_df,
    sintesi_portafoglio = sintesi_df
  ),
  path = "results/copula/risultati_copule.xlsx"
)
