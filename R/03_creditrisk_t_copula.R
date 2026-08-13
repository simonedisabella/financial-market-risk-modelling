# ============================================================
# COSTRUZIONE INPUT CREDITRISK+ CON 20 FASCE
# Gamma-Poisson
# ============================================================
# Descrizione iniziale: 
# il seguente codice usa come riferimento principale l'esercizio in excel sul modello 
# CR + realizzato per il corso di Credit Risk dello scorso anno. 
# (Si veda Lezione16_GP_svolto). Il codice cerca di essere il più fedele possibile 
# all'approccio originale, in particolare per la definizione dei parametri per 
# le fasce.
# Fanno seguito differenze sull'approccio utilizzato data l'introduzione delle 
# copule come strumento utile per cogliere la dipendenza tra i default

# Pulizia ambiente
rm(list = ls())
# Pacchetto necessario
library(mvtnorm)
library(writexl)   # esportazione risultati in Excel

# Cartella di output
output_dir <- "results/credit-risk"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# ------------------------------------------------------------
# 1. Parametri generali
# ------------------------------------------------------------

k <- 20          # numero di fasce
j <- 1:k         # indice delle fasce

# ------------------------------------------------------------
# 2. Scelta di E[N_j]
# ------------------------------------------------------------
# Il numero atteso di default è assunto decrescente
# al crescere della fascia.
# NB: come nell'esercizio excel supponiamo di conoscere a priori valore atteso 
# e varianza per ciascuna fascia, immaginando una loro funzione "deterministica"
# data l'assenza di dati reali a partire da cui sviluppare l'analisi

E_N_j <- 4 * exp(-0.12 * (j - 1)) + 0.15

# ------------------------------------------------------------
# 3. Scelta di Var(N_j)
# ------------------------------------------------------------
# Per una Gamma-Poisson deve valere:
# Var(N_j) > E[N_j]

# imponiamo in modo deterministico una relazione che leghi valore atteso e 
# varianza.
# Assumiamo che le fasce più alte, pur avendo meno default attesi, abbiano 
# una variabilità relativa maggiore, dato il minor numero di soggetti nella fascia

overdispersion_j <- 1.25 + 0.04 * j

V_N_j <- E_N_j * overdispersion_j

# ------------------------------------------------------------
# 4. Calcolo dei parametri con il metodo dei momenti
# ------------------------------------------------------------

# Utilizziamo il metodo dei momenti per stimare i parametri alpha e beta della gamma.
alpha_j <- E_N_j^2 / (V_N_j - E_N_j)
beta_j <- (V_N_j - E_N_j) / E_N_j

# ------------------------------------------------------------
# 5. Creazione tabella input
# ------------------------------------------------------------

input_table <- data.frame(
  j = j,
  E_N_j = E_N_j,
  V_N_j = V_N_j,
  alpha_j = alpha_j,
  beta_j = beta_j
)

# Arrotondamento per visualizzazione
input_table_round <- data.frame(
  j = input_table$j,
  E_N_j = round(input_table$E_N_j, 4),
  V_N_j = round(input_table$V_N_j, 4),
  alpha_j = round(input_table$alpha_j, 4),
  beta_j = round(input_table$beta_j, 6)
)

print(input_table_round)

# ------------------------------------------------------------
# 6. Controlli di coerenza
# ------------------------------------------------------------

# Check di controllo volto a verificare che la varianza sia maggiore del valore
# atteso in ogni fascia

check_variance <- all(V_N_j > E_N_j)

cat("\nControllo Var(N_j) > E[N_j]:", check_variance, "\n")

# ------------------------------------------------------------
# 7. Costruzione delle distribuzioni marginali di N_j
# ------------------------------------------------------------

# come nel caso dell'esercizio di riferimento si suppone un numero massimo di 
# default di 25, oltre il quale viene fatta un'approssimazione. Si considerano
# quindi trascurabili oltre a 25 default per ogni fascia, data la ridotta probabilità.
# NB: l'approccio è naturalmente un'approssimazione, la distribuzione teorica 
# andrebbe da 0 a infinito

n_max <- 25
n_values <- 0:n_max

prob_N <- matrix(NA, nrow = length(n_values), ncol = k)

for (jj in 1:k) {
  prob_N[, jj] <- dnbinom(
    x = n_values,
    size = alpha_j[jj],
    prob = 1 / (1 + beta_j[jj])
  )
}

colnames(prob_N) <- paste0("P_N_", 1:k)
rownames(prob_N) <- n_values

prob_N_table <- data.frame(
  n = n_values,
  prob_N
)

print(round(prob_N_table, 8))

# ============================================================
# VISUALIZZAZIONE DISTRIBUZIONI MARGINALI PER FASCIA
# ============================================================
# Analisi iniziale delle distribuzioni marginali di ogni fascia per valutare 
# differenze all'aumentare del rischio secondo la struttura data originariamente

# Pulizia dispositivi grafici
graphics.off()
par(mfrow = c(1, 1))

# Palette per le 20 fasce
cols_all <- colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
)(k)

# ------------------------------------------------------------
# 1. Griglia 4x5: distribuzione per fascia
# ------------------------------------------------------------
png(
  filename = file.path(output_dir, "01_distribuzioni_marginali_per_fascia.png"),
  width = 1600,
  height = 1200
)
old_par <- par(no.readonly = TRUE)

par(
  mfrow = c(4, 5),
  mar   = c(1.8, 1.8, 2.2, 0.5),
  oma   = c(3, 3, 4, 0.5)
)

for (jj in 1:k) {
  
  bp <- barplot(
    prob_N[, jj],
    names.arg = FALSE,
    col       = cols_all[jj],
    border    = "white",
    main      = paste0(
      "Fascia ", jj,
      "\nE=", round(E_N_j[jj], 2),
      "  α=", round(alpha_j[jj], 2),
      "  β=", round(beta_j[jj], 2)
    ),
    cex.main  = 0.65,
    axes      = FALSE,
    ylim      = c(0, max(prob_N[, jj]) * 1.20)
  )
  
  axis(
    side = 1,
    at = bp[seq(1, length(bp), by = 5)],
    labels = n_values[seq(1, length(n_values), by = 5)],
    cex.axis = 0.45,
    lwd = 0,
    lwd.ticks = 0.5
  )
  
  axis(
    side = 2,
    labels = FALSE,
    lwd = 0,
    lwd.ticks = 0.4
  )
  
  box(col = "gray80", lwd = 0.6)
}

mtext(
  "Distribuzioni marginali di N_j per fascia",
  outer = TRUE,
  cex   = 1.25,
  font  = 2
)

mtext(
  "n",
  side = 1,
  outer = TRUE,
  line = 1.2,
  cex = 0.9
)

mtext(
  expression(P(N[j] == n)),
  side = 2,
  outer = TRUE,
  line = 1.2,
  cex = 0.9
)

par(old_par)
dev.off()
# ------------------------------------------------------------
# 2. Confronto fasce selezionate
# ------------------------------------------------------------
png(
  filename = file.path(output_dir, "02_confronto_marginali_fasce_selezionate.png"),
  width = 1200,
  height = 700
)
par(mfrow = c(1, 1))
par(mar = c(4.5, 4.5, 3.5, 1))

fasce_plot <- c(1, 5, 10, 15, 20)
cols_plot  <- cols_all[fasce_plot]

plot(
  n_values,
  prob_N[, fasce_plot[1]],
  type = "n",
  xlim = c(0, max(n_values)),
  ylim = c(0, max(prob_N[, fasce_plot]) * 1.15),
  xlab = "n  (numero di default)",
  ylab = expression(P(N[j] == n)),
  main = "Confronto tra distribuzioni marginali selezionate",
  cex.lab  = 1.00,
  cex.main = 1.10
)

grid(col = "gray90", lty = 1)

for (idx in seq_along(fasce_plot)) {
  
  jj <- fasce_plot[idx]
  
  lines(
    n_values,
    prob_N[, jj],
    type = "b",
    pch  = 19,
    cex  = 0.65,
    col  = cols_plot[idx],
    lwd  = 2
  )
}

legend(
  "topright",
  legend = paste0(
    "Fascia ", fasce_plot,
    " | E=", round(E_N_j[fasce_plot], 2),
    " | β=", round(beta_j[fasce_plot], 2)
  ),
  col    = cols_plot,
  lty    = 1,
  pch    = 19,
  lwd    = 2,
  bty    = "n",
  cex    = 0.80
)
dev.off()
# ------------------------------------------------------------
# 3. Heatmap probabilità P(N_j = n)
# ------------------------------------------------------------
png(
  filename = file.path(output_dir, "03_heatmap_probabilita_marginali.png"),
  width = 1300,
  height = 750
)
par(mfrow = c(1, 1))
par(mar = c(4.5, 4.5, 3.5, 1))

z_mat <- t(prob_N)

heat_cols <- colorRampPalette(
  c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B")
)(200)

image(
  x    = 1:k,
  y    = n_values,
  z    = z_mat,
  col  = heat_cols,
  xlab = "Fascia j",
  ylab = "n  (numero di default)",
  main = expression("Heatmap delle probabilità marginali  " * P(N[j] == n)),
  xaxt = "n",
  yaxt = "n",
  useRaster = TRUE
)

axis(
  side = 1,
  at = 1:k,
  labels = 1:k,
  cex.axis = 0.75
)

axis(
  side = 2,
  at = seq(0, max(n_values), by = 5),
  labels = seq(0, max(n_values), by = 5),
  las = 1,
  cex.axis = 0.80
)

# Griglia leggera
abline(
  v = seq(1.5, k - 0.5, by = 1),
  col = rgb(1, 1, 1, 0.30),
  lwd = 0.4
)

abline(
  h = seq(2.5, max(n_values) - 0.5, by = 5),
  col = rgb(1, 1, 1, 0.40),
  lwd = 0.5
)

box(col = "gray40")
dev.off()

# ------------------------------------------------------------
# 4. Riepilogo numerico a console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("RIEPILOGO DISTRIBUZIONI MARGINALI\n")
cat("========================================\n")

cat(sprintf(
  "%-8s %-8s %-10s %-8s %-8s\n",
  "Fascia", "E[N_j]", "Var[N_j]", "alpha_j", "beta_j"
))

cat(rep("-", 50), "\n", sep = "")

for (jj in 1:k) {
  cat(sprintf(
    "%-8d %-8.3f %-10.3f %-8.3f %-8.3f\n",
    jj,
    E_N_j[jj],
    V_N_j[jj],
    alpha_j[jj],
    beta_j[jj]
  ))
}

cat("========================================\n")

# Questo blocco ci fa comprendere anche graficamente la logica delle distribuzioni
# costruite per le fasce: quanto più ci si sposta verso le fasce più sicure quanto più le distribuzioni
# divengono J-Shaped (massa concentrata su pochi default). 
# Le meno sicure assumono invece la tipica forma campanulare asimmetrica.

## ============================================================
# 8. Simulazione Monte Carlo con copula t di Student
# ============================================================
# 4 steps:
# 1. Simulazione di variabili dipendenti da una t multivariata
# 2. Trasformazione in uniformi (Teorema di trasformazione integrale di probabilità)
# 3. Trasformazione delle uniformi in marginali Gamma Poisson (Metodo d'inversione)
# 4. Calcolo della perdita aggregata

# fissato un seed per garantire la riproducibilità dell'esperimento
set.seed(123)

# Numero di simulazioni Monte Carlo
# Al diminuire del numero di simulazioni naturalmente si hanno risultati meno stabili
# viceversa all'aumentare delle simulazioni. 
B <- 50000

# Gradi di libertà della copula t di Student
nu <- 3

# Livello di confidenza per VaR e TCE
alpha_level <- 0.99

# Costante base di perdita
L <- 1000

# Perdita associata a ciascuna fascia: calcolata come perdita rappresentativa per
# fascia calcolata tramite banding
J <- 1:k
loss_j <- L * J

# Griglia di valori di Rho analizzata al fine di apprezzare il cambiamento delle
# misure di rischio al variare del parametro.
# La comonotonia viene simulata separatamente rispetto alla copula t.
# Infatti, rho = 1 renderebbe la matrice di equi-correlazione singolare.
# Il caso comonotono è quindi trattato esternamente come benchmark limite di dipendenza perfetta,
# ottenuto imponendo la stessa uniforme comune a tutte le marginali.
# Partiamo da 0 poichè l'interesse principale è valutare la dipendenza positiva tra default

rho_values <- seq(0, 0.9, by = 0.1)

# Tabella vuota per salvare VaR, TCE e statistiche sintetiche
results_rho <- data.frame(
  scenario = character(),
  rho = numeric(),
  VaR = numeric(),
  TCE = numeric(),
  mean_loss = numeric(),
  sd_loss = numeric()
)

# Lista per salvare l'intera distribuzione simulata delle perdite aggregate Z
# per ciascun valore di rho
Z_list <- list()

# ------------------------------------------------------------
# 8.1 Simulazione al variare di rho
# ------------------------------------------------------------
# Applicazione inversa del teorema di Sklar: costruzione della distribuzione 
# congiunta a partire dalla copula scelta.

for (rho in rho_values) {
  
  cat("\nSimulazione con rho =", rho, "\n")
  
  # Matrice di equi-correlazione
  Corr <- matrix(rho, nrow = k, ncol = k)
  diag(Corr) <- 1
  
  # Passo 1: simulazione dalla t multivariata
  X <- rmvt(B, sigma = Corr, df = nu)
  
  # Passo 2: trasformazione in uniformi, pt = probabilità cumulata t di student
  U <- pt(X, df = nu)
  
  # Passo 3: trasformazione nelle marginali Gamma-Poisson
  N_sim <- matrix(NA, nrow = B, ncol = k)
  
  # blocco chiave: Inversione della CDF della NegBinomiale sulle uniformi.
  # Scelta della NegBinomiale giustificata dal programma di Credit Risk, in particolare
  # nella sezione sul miscuglio Gamma-Poisson.
  # Trasformazione delle uniformi dipendenti (copula t) nelle marginali Gamma-Poisson 
  # desiderate, preservando la struttura di dipendenza
  
  for (i in 1:k) {
    N_sim[, i] <- qnbinom(
      p = U[, i],
      size = alpha_j[i],
      prob = 1 / (1 + beta_j[i])
    )
  }
  
  # Perdita aggregata: calcolata come perdita rappresentativa x conteggio dei default
  Z <- as.vector(N_sim %*% loss_j)
  
  # Salvataggio della distribuzione completa delle perdite aggregate
  Z_list[[paste0("rho = ", rho)]] <- Z
  
  # VaR e TCE
  VaR_rho <- as.numeric(quantile(Z, probs = alpha_level, type = 1))
  TCE_rho <- mean(Z[Z >= VaR_rho])
  
  # Salvataggio risultati sintetici
  results_rho <- rbind(
    results_rho,
    data.frame(
      scenario = paste0("rho = ", rho),
      rho = rho,
      VaR = VaR_rho,
      TCE = TCE_rho,
      mean_loss = mean(Z),
      sd_loss = sd(Z)
    )
  )
}

# ============================================================
# 9. Caso di comonotonia
# ============================================================

cat("\nSimulazione caso comonotono\n")

# Limite superiore di Frechet
# Nella comonotonia tutte le marginali sono guidate
# dalla stessa variabile uniforme
# runif(B) genera un vettore di lunghezza B di variabili casuali uniformi tra 0 e 1

U_common <- runif(B)

N_comono <- matrix(NA, nrow = B, ncol = k)

for (i in 1:k) {
  N_comono[, i] <- qnbinom(
    p = U_common,
    size = alpha_j[i],
    prob = 1 / (1 + beta_j[i])
  )
}

# Perdita aggregata nel caso comonotono
Z_comono <- as.vector(N_comono %*% loss_j)

# VaR e TCE nel caso comonotono
VaR_comono <- as.numeric(quantile(Z_comono, probs = alpha_level, type = 1))
TCE_comono <- mean(Z_comono[Z_comono >= VaR_comono])

result_comono <- data.frame(
  scenario = "Comonotonia",
  rho = NA,
  VaR = VaR_comono,
  TCE = TCE_comono,
  mean_loss = mean(Z_comono),
  sd_loss = sd(Z_comono)
)

# ============================================================
# 10. Tabella finale dei risultati
# ============================================================

# Unione dei risultati ottenuti con il primo ciclo e quelli ottenuti con il caso
# comonotono

results <- rbind(results_rho, result_comono)
print(results)

# ============================================================
# 11. istogrammi delle perdite aggregate
# ============================================================
# Aggiungo anche la comonotonia alla lista

Z_list[["Comonotonia"]] <- Z_comono

# Ordine corretto degli scenari
scenario_names <- c(
  paste0("rho = ", rho_values),
  "Comonotonia"
)

# Livello di confidenza
alpha <- 0.99

png(
  filename = file.path(output_dir, "04_istogrammi_perdite_aggregate_tutti_scenari.png"),
  width = 1400,
  height = 1000
)

# Layout grafico: 4 righe x 3 colonne
par(
  mfrow = c(4, 3),
  mar = c(4, 4, 3, 1),
  oma = c(0, 0, 3, 0)
)

# Istogramma per ciascuno scenario
for (sc in scenario_names) {
  
  Z_sc <- Z_list[[sc]]
  
  # Calcolo VaR e TCE dello scenario
  VaR_sc <- as.numeric(quantile(Z_sc, probs = alpha, type = 1))
  TCE_sc <- mean(Z_sc[Z_sc >= VaR_sc])
  
  hist(
    Z_sc,
    breaks = 60,
    main = sc,
    xlab = "Perdita aggregata Z",
    col = "lightblue",
    border = "white"
  )
  
  # Linea VaR
  abline(
    v = VaR_sc,
    col = "red",
    lwd = 2,
    lty = 2
  )
  
  # Linea TCE
  abline(
    v = TCE_sc,
    col = "darkblue",
    lwd = 2,
    lty = 3
  )
  
  legend(
    "topright",
    legend = c("VaR 99%", "TCE 99%"),
    col = c("red", "darkblue"),
    lwd = 2,
    lty = c(2, 3),
    cex = 0.7,
    bty = "n"
  )
}

# Titolo generale
mtext(
  "Distribuzione delle perdite aggregate per tutti i livelli di dipendenza",
  outer = TRUE,
  cex = 1.4,
  font = 2
)
dev.off()
# ============================================================
# GRAFICI VaR e TCE
# ============================================================

# ------------------------------------------------------------
# 2. Dataset per includere la comonotonia come rho = 1
# ------------------------------------------------------------

results_static <- rbind(
  results_rho[, c("scenario", "rho", "VaR", "TCE", "mean_loss", "sd_loss")],
  data.frame(
    scenario  = "Comonotonia",
    rho       = 1,
    VaR       = VaR_comono,
    TCE       = TCE_comono,
    mean_loss = mean(Z_comono),
    sd_loss   = sd(Z_comono)
  )
)

# ============================================================
# 3. Grafico VaR con comonotonia in rho = 1
# ============================================================

png(
  filename = file.path(output_dir, "05_var_al_variare_della_dipendenza.png"),
  width = 1100,
  height = 700
)

par(mar = c(5, 5.5, 4, 2))

ylim_var <- range(results_static$VaR) * c(0.97, 1.03)

plot(
  results_rho$rho,
  results_rho$VaR,
  type = "n",
  xlim = c(0, 1),
  ylim = ylim_var,
  xlab = expression(rho),
  ylab = "VaR",
  main = "VaR al variare della dipendenza",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 1.1,
  cex.main = 1.2
)

axis(
  side = 1,
  at = seq(0, 1, by = 0.1),
  labels = seq(0, 1, by = 0.1),
  cex.axis = 0.85
)

axis(
  side = 2,
  at = pretty(ylim_var, n = 5),
  labels = formatC(pretty(ylim_var, n = 5), format = "f", digits = 0, big.mark = ","),
  las = 1,
  cex.axis = 0.85
)

grid(col = "gray90", lty = 1)

lines(
  results_rho$rho,
  results_rho$VaR,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#1F4E79"
)

segments(
  x0 = max(results_rho$rho),
  y0 = results_rho$VaR[which.max(results_rho$rho)],
  x1 = 1,
  y1 = VaR_comono,
  lty = 2,
  lwd = 1.8,
  col = "#B2182B"
)

points(
  1,
  VaR_comono,
  pch = 17,
  cex = 1.5,
  col = "#B2182B"
)

text(
  x = 1,
  y = VaR_comono,
  labels = "Comonotonia",
  pos = 3,
  cex = 0.85,
  col = "#B2182B"
)

legend(
  "topleft",
  legend = c("Copula t", "Comonotonia"),
  col = c("#1F4E79", "#B2182B"),
  lty = c(1, 2),
  pch = c(19, 17),
  lwd = c(2, 1.8),
  bty = "n",
  cex = 0.9
)

box(col = "gray50")
dev.off()

# ============================================================
# 4. Grafico TCE con comonotonia in rho = 1
# ============================================================

png(
  filename = file.path(output_dir, "06_tce_al_variare_della_dipendenza.png"),
  width = 1100,
  height = 700
)
par(mar = c(5, 5.5, 4, 2))

ylim_tce <- range(results_static$TCE) * c(0.97, 1.03)

plot(
  results_rho$rho,
  results_rho$TCE,
  type = "n",
  xlim = c(0, 1),
  ylim = ylim_tce,
  xlab = expression(rho),
  ylab = "TCE",
  main = "TCE al variare della dipendenza",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 1.1,
  cex.main = 1.2
)

axis(
  side = 1,
  at = seq(0, 1, by = 0.1),
  labels = seq(0, 1, by = 0.1),
  cex.axis = 0.85
)

axis(
  side = 2,
  at = pretty(ylim_tce, n = 5),
  labels = formatC(pretty(ylim_tce, n = 5), format = "f", digits = 0, big.mark = ","),
  las = 1,
  cex.axis = 0.85
)

grid(col = "gray90", lty = 1)

lines(
  results_rho$rho,
  results_rho$TCE,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#1F4E79"
)

segments(
  x0 = max(results_rho$rho),
  y0 = results_rho$TCE[which.max(results_rho$rho)],
  x1 = 1,
  y1 = TCE_comono,
  lty = 2,
  lwd = 1.8,
  col = "#B2182B"
)

points(
  1,
  TCE_comono,
  pch = 17,
  cex = 1.5,
  col = "#B2182B"
)

text(
  x = 1,
  y = TCE_comono,
  labels = "Comonotonia",
  pos = 3,
  cex = 0.85,
  col = "#B2182B"
)

legend(
  "topleft",
  legend = c("Copula t", "Comonotonia"),
  col = c("#1F4E79", "#B2182B"),
  lty = c(1, 2),
  pch = c(19, 17),
  lwd = c(2, 1.8),
  bty = "n",
  cex = 0.9
)

box(col = "gray50")
dev.off()

# ============================================================
# 5. Barplot VaR e TCE affiancate
# ============================================================

png(
  filename = file.path(output_dir, "07_barplot_var_tce_dipendenza.png"),
  width = 1300,
  height = 750
)

par(mar = c(6, 5.5, 4, 2))

bar_data <- rbind(
  VaR = results_static$VaR,
  TCE = results_static$TCE
)

colnames(bar_data) <- c(
  paste0("rho=", results_rho$rho),
  "Comonotonia"
)

barplot(
  bar_data,
  beside = TRUE,
  col = c("#1F4E79", "#B2182B"),
  border = "white",
  ylim = c(0, max(bar_data) * 1.15),
  ylab = "Misura di rischio",
  main = "VaR e TCE al variare della dipendenza",
  las = 2,
  cex.names = 0.75,
  cex.axis = 0.85
)

grid(
  nx = NA,
  ny = NULL,
  col = "gray90",
  lty = 1
)

legend(
  "topleft",
  legend = c("VaR", "TCE"),
  fill = c("#1F4E79", "#B2182B"),
  bty = "n",
  cex = 0.9
)

box(col = "gray50")
dev.off()

# ============================================================
# 12. Esportazione risultati in Excel
# ============================================================

# Tabella input delle fasce
input_table_excel <- data.frame(
  fascia = j,
  E_N_j = E_N_j,
  V_N_j = V_N_j,
  overdispersion_j = overdispersion_j,
  alpha_j = alpha_j,
  beta_j = beta_j,
  loss_j = loss_j
)

# Tabella delle probabilità marginali
prob_N_excel <- data.frame(
  n = n_values,
  prob_N
)

# Risultati principali VaR, TCE, media e deviazione standard
results_excel <- results

# Versione con comonotonia trattata graficamente come rho = 1
results_static_excel <- results_static

# Parametri generali della simulazione
simulation_settings <- data.frame(
  parametro = c(
    "numero_fasce",
    "numero_simulazioni_MC",
    "gradi_liberta_copula_t",
    "livello_confidenza",
    "perdita_base_L",
    "n_max_marginali",
    "seed"
  ),
  valore = c(
    k,
    B,
    nu,
    alpha_level,
    L,
    n_max,
    123
  )
)

# Sintesi delle distribuzioni delle perdite aggregate
loss_distribution_summary <- data.frame()

for (sc in scenario_names) {
  
  Z_sc <- Z_list[[sc]]
  
  loss_distribution_summary <- rbind(
    loss_distribution_summary,
    data.frame(
      scenario = sc,
      min_loss = min(Z_sc),
      q_50 = as.numeric(quantile(Z_sc, 0.50, type = 1)),
      q_95 = as.numeric(quantile(Z_sc, 0.95, type = 1)),
      q_99 = as.numeric(quantile(Z_sc, 0.99, type = 1)),
      q_995 = as.numeric(quantile(Z_sc, 0.995, type = 1)),
      max_loss = max(Z_sc),
      mean_loss = mean(Z_sc),
      sd_loss = sd(Z_sc)
    )
  )
}

# Esportazione Excel
write_xlsx(
  list(
    input_fasce = input_table_excel,
    probabilita_marginali = prob_N_excel,
    risultati_var_tce = results_excel,
    risultati_grafici = results_static_excel,
    sintesi_distribuzioni_perdite = loss_distribution_summary,
    impostazioni_simulazione = simulation_settings
  ),
  path = file.path(output_dir, "Risultati_creditrisk_gamma_poisson.xlsx")
)
# ============================================================
# 13. Apertura esterna dei grafici salvati
# ============================================================

# Funzione per aprire un file con il visualizzatore esterno del sistema operativo
open_external <- function(file_path) {
  
  file_path <- normalizePath(file_path, winslash = "/", mustWork = TRUE)
  
  if (.Platform$OS.type == "windows") {
    
    shell.exec(file_path)
    
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    
    system2("open", shQuote(file_path))
    
  } else {
    
    system2("xdg-open", shQuote(file_path))
  }
}

# Lista dei grafici PNG salvati nella cartella di output
plot_files <- list.files(
  path = output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

plot_files <- sort(plot_files)

# Apertura esterna dei grafici
if (length(plot_files) == 0) {
  
  warning("Nessun grafico PNG trovato nella cartella di output.")
  
} else {
  
  cat("\nApertura esterna dei grafici salvati...\n")
  
  for (f in plot_files) {
    open_external(f)
    Sys.sleep(0.5)
  }
}