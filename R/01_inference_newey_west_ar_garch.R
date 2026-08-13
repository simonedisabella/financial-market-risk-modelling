# ASSIGNMENT ESERCIZIO 1 - INFERENZA
# Librerie necessarie
library(readxl)
library(fGarch)


#########################################################################################
# PUNTO 1
# Funzioni R per stime Newey-West della varianza asintotica:
# 1) dello stimatore dello scarto quadratico medio
# 2) dello stimatore del VaR
#########################################################################################
# 1) Funzione Newey-West
NeweyWest_varianza_lungo_periodo <- function(serie_osservata_x, numero_lag_NeweyWest_m = NULL) {
  # x: serie osservata
  # m: numero massimo di ritardi Newey-West
  numero_osservazioni_n <- length(serie_osservata_x)
  # se non specificato il numero dei lag è scelto automaticamente 
  if (is.null(numero_lag_NeweyWest_m)) {
    # Rule of thumb di Stock-Watson
    # floor = arrotonda un numero per difetto all'intero più vicino.
    numero_lag_NeweyWest_m <- floor(4 * (numero_osservazioni_n / 100)^(2 / 9))
  }
  numero_lag_NeweyWest_m <- min(numero_lag_NeweyWest_m, numero_osservazioni_n - 1)
  media_serie_osservata_mu <- mean(serie_osservata_x)
  # Autocovarianza campionaria a lag 0: coincide con la varianza campionaria della serie perché la serie viene confrontata con sé stessa senza ritardo temporale.
  varianza_stimata_v <- mean((serie_osservata_x - media_serie_osservata_mu)^2)
  # Autocovarianze ai lag successivi con pesi Bartlett
  if (numero_lag_NeweyWest_m >= 1) {
    for (j in 1:numero_lag_NeweyWest_m) {
      autocovarianza_j <- sum((serie_osservata_x[(j + 1):numero_osservazioni_n] - media_serie_osservata_mu) * (serie_osservata_x[1:(numero_osservazioni_n - j)] - media_serie_osservata_mu)) / numero_osservazioni_n
      # Peso Bartlett Newey-West: w_j = 1 - j / (m + 1), con m massimo lag incluso
      peso_Bartlett_j <- 1 - j / (numero_lag_NeweyWest_m + 1)
      varianza_stimata_v <- varianza_stimata_v + 2 * peso_Bartlett_j * autocovarianza_j
    }
  } 
  return(varianza_stimata_v)
}


# 2) Newey-West per lo stimatore dello scarto quadratico medio
NeweyWest_scarto_quadratico_medio <- function(serie_osservata_x, numero_lag_NeweyWest_m = NULL) {
  # x: serie storica dei rendimenti_simulati_AR1
  # m: numero massimo di ritardi Newey-West
  numero_osservazioni_n <- length(serie_osservata_x)
  media_campionaria <- mean(serie_osservata_x)
  # Varianza campionaria corretta: denominatore n - 1
  varianza_campionaria_corretta <- sum((serie_osservata_x - media_campionaria)^2) / (numero_osservazioni_n - 1)
  # Stimatore dello scarto quadratico medio
  scarto_quadratico_medio_campionario <- sqrt(varianza_campionaria_corretta)
  # Serie trasformata: scarti quadratici
  scarti_quadratici_dalla_media_campionaria_y <- (serie_osservata_x - media_campionaria)^2
  # Varianza di lungo periodo Newey-West degli scarti quadratici
  varianza_scarti_quadratici_y <- NeweyWest_varianza_lungo_periodo(scarti_quadratici_dalla_media_campionaria_y, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
  # Metodo delta:
  # g(v) = sqrt(v)
  # g'(v) = 1 / (2 sqrt(v))
  # AVAR(sqrt(varianza_campionaria_corretta)) = AVAR(varianza_campionaria_corretta) / (4 * var)           AVAR = varianza asintotica
  varianza_asintotica_scarto_quadratico_medio <- varianza_scarti_quadratici_y / (4 * varianza_campionaria_corretta)
  # Varianza stimata dello stimatore dello Scarto Quadratico Medio (SQM)
  # Divido la varianza asintotica per n perché lo stimatore ha convergenza a velocità sqrt(n).
  varianza_stimata_stimatore_SQM <- varianza_asintotica_scarto_quadratico_medio / numero_osservazioni_n
  # Errore standard stimato SQM
  standard_error_stimato_SQM <- sqrt(varianza_stimata_stimatore_SQM)
  risultati <- list(
    numero_osservazioni_n = numero_osservazioni_n,
    # ifelse = sceglie il primo valore se la condizione è TRUE, sceglie il secondo valore se la condizione è FALSE.
    # floor = arrotonda un numero per difetto all'intero più vicino.
    numero_lag_NeweyWest_m = ifelse(is.null(numero_lag_NeweyWest_m), floor(4 * (numero_osservazioni_n / 100)^(2 / 9)), numero_lag_NeweyWest_m),
    media_campionaria = media_campionaria,
    varianza_campionaria_corretta = varianza_campionaria_corretta,
    scarto_quadratico_medio_campionario = scarto_quadratico_medio_campionario,
    varianza_scarti_quadratici_y = varianza_scarti_quadratici_y,
    varianza_asintotica_scarto_quadratico_medio = varianza_asintotica_scarto_quadratico_medio,
    varianza_stimata_stimatore_SQM = varianza_stimata_stimatore_SQM,
    standard_error_stimato_SQM = standard_error_stimato_SQM
  )
  return(risultati)
}


# 3) Newey-West per lo stimatore del VaR
NeweyWest_VaR <- function(serie_osservata_x, livello_quantile_p = 0.05, numero_lag_NeweyWest_m = NULL) {
  # x: serie storica dei rendimenti_simulati_AR1
  # p: livello del quantile
  # m: numero massimo di ritardi Newey-West
  numero_osservazioni_n <- length(serie_osservata_x)
  # Quantile campionario
  quantile_campionario <- as.numeric(quantile(serie_osservata_x, probs = livello_quantile_p))
  # VaR lasciato con il segno del quantile dei rendimenti.
  VaR_stimato <- quantile_campionario
  # Serie indicatrice usata per la rappresentazione asintotica del quantile.
  # La media di questa serie corrisponde alla frequenza empirica delle osservazioni nella coda sinistra.
  indicatore <- serie_osservata_x <= quantile_campionario
  # Varianza di lungo periodo Newey-West dell'indicatore.
  # Serve per tenere conto dell'eventuale dipendenza seriale negli eventi di coda.
  varianza_indicatore_quantile <- NeweyWest_varianza_lungo_periodo(indicatore, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
  # Stima kernel non parametrica della densità della serie osservata.
  # La densità è necessaria perché la varianza asintotica del quantile dipende dalla densità nel punto del quantile.  
  # density = stima la densità di una variabile continua a partire dai dati campionari.
  densita_kernel_stimata <- density(serie_osservata_x)
  # Densità stimata con metodo kernel, valutata nel punto del quantile campionario tramite interpolazione.
  # approx = interpola linearmente tra i punti disponibili; qui serve per ottenere la densità stimata esattamente nel punto del quantile campionario.
  densita_kernel_stimata_quantile <- approx(densita_kernel_stimata$x, densita_kernel_stimata$y, xout = quantile_campionario)$y
  # Metodo delta per il quantile: la varianza asintotica del quantile stimato è la varianza di lungo periodo dell'indicatore divisa per densità f(q)^2.
  # AVAR(q_hat) = LRV(1{x <= q}) / f(q)^2
  varianza_asintotica_VaR <- varianza_indicatore_quantile / (densita_kernel_stimata_quantile^2)
  # Divido per n perché la varianza dello stimatore è la varianza asintotica divisa per la numerosità campionaria.
  var_stimatore_VaR <- varianza_asintotica_VaR / numero_osservazioni_n
  # Errore standard stimato
  standard_error_VaR <- sqrt(var_stimatore_VaR)
  risultati <- list(
    numero_osservazioni_n = numero_osservazioni_n,
    # ifelse = sceglie il primo valore se la condizione è TRUE, sceglie il secondo valore se la condizione è FALSE.
    # floor = arrotonda un numero per difetto all'intero più vicino.
    numero_lag_NeweyWest_m = ifelse(is.null(numero_lag_NeweyWest_m), floor(4 * (numero_osservazioni_n / 100)^(2 / 9)), numero_lag_NeweyWest_m),
    livello_quantile_p = livello_quantile_p,
    quantile_campionario = quantile_campionario,
    VaR_stimato = VaR_stimato,
    densita_kernel_stimata_quantile = densita_kernel_stimata_quantile,
    varianza_indicatore_quantile = varianza_indicatore_quantile,
    varianza_asintotica_VaR = varianza_asintotica_VaR,
    var_stimatore_VaR = var_stimatore_VaR,
    standard_error_VaR = standard_error_VaR
  )
  return(risultati)
}


#########################################################################################
# PUNTO 2
# Analisi di una serie storica di rendimenti azionari
#########################################################################################
# INTERVALLI DI CONFIDENZA PER LA SIMULAZIONE
# 1) Intervallo di confidenza per il valore atteso.
# Se correggi_dipendenza = TRUE uso Newey-West del punto 1.
# Se correggi_dipendenza = FALSE uso la formula iid, cioè ignoro erroneamente la dipendenza seriale.
intervallo_confidenza_media_simulazione <- function(serie_osservata_x, correggi_dipendenza, livello_significativita, numero_lag_NeweyWest_m) {
  numero_osservazioni_n <- length(serie_osservata_x)
  media_campionaria <- mean(serie_osservata_x)
  # qnorm = restituisce il quantile della distribuzione normale standard.
  quantile_normale_bilaterale <- qnorm(1 - livello_significativita / 2)
  if (correggi_dipendenza == TRUE) {
    varianza_asintotica_media <- NeweyWest_varianza_lungo_periodo(serie_osservata_x, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
  } else {
    # Formula iid: considero solo la varianza a lag 0 e ignoro le autocovarianze.
    varianza_asintotica_media <- mean((serie_osservata_x - media_campionaria)^2)
  }
  standard_error_media <- sqrt(varianza_asintotica_media / numero_osservazioni_n)
  intervallo_di_confidenza <- c(
    limite_inferiore = media_campionaria - quantile_normale_bilaterale * standard_error_media,
    limite_superiore = media_campionaria + quantile_normale_bilaterale * standard_error_media,
    ampiezza = 2 * quantile_normale_bilaterale * standard_error_media
  )
  return(intervallo_di_confidenza)
}


# 2) Intervallo di confidenza per lo scarto quadratico medio.
# Se correggi_dipendenza = TRUE uso la funzione Newey-West del punto 1.
# Se correggi_dipendenza = FALSE uso la formula iid, cioè ignoro erroneamente la dipendenza seriale.
intervallo_confidenza_SQM_simulazione <- function(serie_osservata_x, correggi_dipendenza, livello_significativita, numero_lag_NeweyWest_m) {
  numero_osservazioni_n <- length(serie_osservata_x)
  media_campionaria <- mean(serie_osservata_x)
  scarti_quadratici <- (serie_osservata_x - media_campionaria)^2
  varianza_campionaria_corretta <- sum(scarti_quadratici) / (numero_osservazioni_n - 1)
  scarto_quadratico_medio_campionario <- sqrt(varianza_campionaria_corretta)
  # qnorm = restituisce il quantile della distribuzione normale standard.
  quantile_normale_bilaterale <- qnorm(1 - livello_significativita / 2)
  if (correggi_dipendenza == TRUE) {
    risultato_NeweyWest_SQM <- NeweyWest_scarto_quadratico_medio(serie_osservata_x, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    standard_error_SQM <- risultato_NeweyWest_SQM$standard_error_stimato_SQM
  } else {
    # Formula iid: uso solo la varianza degli scarti quadratici e ignoro le autocovarianze.
    varianza_scarti_quadratici_iid <- mean((scarti_quadratici - mean(scarti_quadratici))^2)
    varianza_asintotica_SQM_iid <- varianza_scarti_quadratici_iid / (4 * varianza_campionaria_corretta)
    standard_error_SQM <- sqrt(varianza_asintotica_SQM_iid / numero_osservazioni_n)
  }
  intervallo_di_confidenza <- c(
    limite_inferiore = scarto_quadratico_medio_campionario - quantile_normale_bilaterale * standard_error_SQM,
    limite_superiore = scarto_quadratico_medio_campionario + quantile_normale_bilaterale * standard_error_SQM,
    ampiezza = 2 * quantile_normale_bilaterale * standard_error_SQM
  )
  return(intervallo_di_confidenza)
}


# 3) Intervallo di confidenza per il VaR.
# Il VaR viene espresso come: VaR = quantile al livello p.
# Se correggi_dipendenza = TRUE uso la funzione Newey-West del punto 1.
# Se correggi_dipendenza = FALSE uso la formula iid, cioè ignoro erroneamente la dipendenza seriale.
intervallo_confidenza_VaR_simulazione <- function(serie_osservata_x, livello_quantile_p, correggi_dipendenza, livello_significativita, numero_lag_NeweyWest_m) {
  numero_osservazioni_n <- length(serie_osservata_x)
  quantile_normale_bilaterale <- qnorm(1 - livello_significativita / 2)
  # Quantile campionario della coda sinistra.
  quantile_campionario <- as.numeric(quantile(serie_osservata_x, probs = livello_quantile_p))
  # VaR lasciato con il segno del quantile dei rendimenti.
  VaR_stimato <- quantile_campionario
  if (correggi_dipendenza == TRUE) {
    risultato_NeweyWest_VaR <- NeweyWest_VaR(serie_osservata_x, livello_quantile_p = livello_quantile_p, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    standard_error_VaR <- risultato_NeweyWest_VaR$standard_error_VaR
  } else {
    # Formula iid: uso solo la varianza dell'indicatore del quantile e ignoro le autocovarianze.
    indicatore <- serie_osservata_x <= quantile_campionario
    varianza_indicatore_iid <- mean((indicatore - mean(indicatore))^2)
    # Stima kernel della densità e valutazione della densità nel punto del quantile campionario.
    densita_kernel_stimata <- density(serie_osservata_x)
    densita_kernel_stimata_quantile <- approx(densita_kernel_stimata$x, densita_kernel_stimata$y, xout = quantile_campionario)$y
    varianza_asintotica_VaR_iid <- varianza_indicatore_iid / (densita_kernel_stimata_quantile^2)
    standard_error_VaR <- sqrt(varianza_asintotica_VaR_iid / numero_osservazioni_n)
  }
  intervallo_di_confidenza <- c(
    limite_inferiore = VaR_stimato - quantile_normale_bilaterale * standard_error_VaR,
    limite_superiore = VaR_stimato + quantile_normale_bilaterale * standard_error_VaR,
    ampiezza = 2 * quantile_normale_bilaterale * standard_error_VaR
  )
  return(intervallo_di_confidenza)
}


#########################################################################################
# PUNTO 3
# Esperimento Monte Carlo con simulazione da AR(1)-GARCH(1,1) con innovazioni non gaussiane.
# Confrontare copertura effettiva e ampiezza degli intervalli di confidenza per:
# 1) valore atteso
# 2) scarto quadratico medio
# 3) VaR
# Distinzione tra:
# a) considerando la dipendenza seriale tramite stima Newey-West
# b) ignorando erroneamente la dipendenza seriale, cioè trattando i dati come iid
#########################################################################################
# FUNZIONE DI SIMULAZIONE AR(1)-GARCH(1,1)
#   rendimento_t  = mu + phi * rendimento_t-1 + shock_t 
#   shock_t = sqrt(varianza_condizionata_t) * epsilon_t
#   varianza_condizionata_t = omega + alpha * (shock_t-1)^2 + beta * varianza_condizionata_t-1
# dove epsilon_t è non gaussiana e standardizzata con media 0 e varianza 1.
simula_ar1_garch11 <- function(n, mu, phi, omega, alpha, beta, innovazione, gradi_di_liberta, da_escludere) {
  if (abs(phi) >= 1) {stop("Per stazionarietà AR(1) deve valere valore assoluto di phi < 1.")}
  if (omega <= 0) {stop("Omega deve essere positivo.")}
  if (alpha < 0 || beta < 0) {stop("Alpha e beta devono essere non negativi.")}
  if (alpha + beta >= 1) {stop("Per stazionarietà GARCH deve valere alpha + beta < 1.")}
  N <- n + da_escludere
  # numeric = crea un vettore numerico; se le si passa una lunghezza, inizializza un vettore di zeri.  
  rendimento <- numeric(N)
  shock <- numeric(N)
  varianza_condizionata <- numeric(N)
  # epsilon = componente casuale di base 
  innovazioni_standardizzate_epsilon <- numeric(N)
  
  # Valori iniziali (uso la media e la varianza incondizionata teoriche).
  rendimento[1] <- mu / (1 - phi)
  shock[1] <- 0
  varianza_condizionata[1] <- omega / (1 - alpha - beta)
  for (i in 2:N) {
    # Innovazioni t-Student standardizzate.
    # La t-Student ha già media 0, ma ha una varianza pari a gradi_di_liberta / (gradi_di_liberta - 2).
    # Divido per il suo scarto quadratico medio teorico così l'innovazione finale ha media 0 e varianza 1.
    # Lo scopo è quello di simulare innovazioni non gaussiane con code pesanti mantenendo però la stessa scala di varianza del modello GARCH.
    if (innovazione == "t") {
      # Genero una singola innovazione t-Student con gradi_di_liberta.
      # rt = genera numeri casuali da una distribuzione t-Student.
      innovazioni_standardizzate_epsilon[i] <- rt(1, df = gradi_di_liberta) / sqrt(gradi_di_liberta / (gradi_di_liberta - 2))
    }
    
    # Innovazioni chi-quadro standardizzate.
    # La chi-quadro ha media pari a gradi_di_liberta e ha varianza pari a 2 * gradi_di_liberta.
    # Sottraggo la sua media teorica e divido per il suo scarto quadratico medio teorico, così ottengo un'innovazione standardizzata con media 0 e varianza 1.
    # Lo scopo è quello di simulare innovazioni non gaussiane e asimmetriche, mantenendo però la stessa scala di varianza del modello GARCH.
    if (innovazione == "chisq") {
      # Genero una singola variabile chi-quadro con gradi_di_liberta.
      # rchisq = genera numeri casuali da una distribuzione chi-quadro.
      innovazione_standardizzate_chi_quadro_z <- rchisq(1, df = gradi_di_liberta)
      innovazioni_standardizzate_epsilon[i] <- (innovazione_standardizzate_chi_quadro_z - gradi_di_liberta) / sqrt(2 * gradi_di_liberta)
    }
    
    # Ricorsione GARCH(1,1).
    # La varianza al tempo i dipende da:
    # - omega: livello costante della varianza;
    # - alpha * shock[i - 1]^2: effetto dello shock quadratico del periodo precedente;
    # - beta * varianza_condizionata[i - 1]: persistenza della varianza passata.
    # Lo scopo è generare volatilità variabile nel tempo e clustering della volatilità.
    varianza_condizionata[i] <- omega + alpha * shock[i - 1]^2 + beta * varianza_condizionata[i - 1]
    # Shock eteroschedastico.
    # Moltiplico l'innovazione standardizzata (con media 0 e varianza 1) per la deviazione standard condizionata; in questo modo lo shock ha media condizionata 0 e varianza condizionata pari a varianza_condizionata[i].
    shock[i] <- sqrt(varianza_condizionata[i]) * innovazioni_standardizzate_epsilon[i]
    # Ricorsione AR(1).
    # Il rendimento al tempo i dipende da:
    # - mu: termine costante;
    # - phi * rendimento[i - 1]: dipendenza lineare dal rendimento precedente;
    # - shock[i]: nuova componente casuale eteroschedastica.
    # Lo scopo è simulare una serie con dipendenza nella media e volatilità condizionata variabile nel tempo.
    rendimento[i] <- mu + phi * rendimento[i - 1] + shock[i]
  }
  # Elimino i valori da escludere (indicato con da_escludere) e restituisco solo le osservazioni finali.
  return(rendimento[(da_escludere + 1):N])
}


# CALCOLO DEI VERI VALORI AR(1)-GARCH(1,1)
# Calcolo i valori veri del processo simulato per valutare la copertura degli intervalli di confidenza.
# La copertura risponde alla domanda: "Quante volte l'intervallo contiene il vero valore del parametro?"
# Restituisce:
# - media vera
# - scarto quadratico medio vero
# - VaR vero approssimato tramite simulazione lunga
parametri_veri_ar1_garch11 <- function(mu, phi, omega, alpha, beta, innovazione, gradi_di_liberta, livello_quantile_p, numero_osservazioni_per_VaR_vero, da_escludere_per_VaR_vero) {
  # Calcolo la serie_lunga per il VaR vero.
  serie_lunga <- simula_ar1_garch11(n = numero_osservazioni_per_VaR_vero, mu = mu, phi = phi, omega = omega, alpha = alpha, beta = beta, innovazione = innovazione, gradi_di_liberta = gradi_di_liberta, da_escludere = da_escludere_per_VaR_vero)
  # Media incondizionata teorica del processo AR(1).
  media_vera <- mu / (1 - phi)
  # Varianza incondizionata teorica dello shock GARCH.
  varianza_incondizionata_shock <- omega / (1 - alpha - beta)
  # Varianza incondizionata teorica del rendimento AR(1).
  varianza_incondizionata_rendimento <- varianza_incondizionata_shock / (1 - phi^2)
  # Scarto quadratico medio vero del rendimento.
  scarto_quadratico_medio_vero <- sqrt(varianza_incondizionata_rendimento)
  # VaR vero approssimato tramite simulazione lunga: VaR = quantile_p.
  VaR_vero <- as.numeric(quantile(serie_lunga, probs = livello_quantile_p))
  return(list(media_vera = media_vera, scarto_quadratico_medio_vero = scarto_quadratico_medio_vero, VaR_vero = VaR_vero))
}


# MONTECARLO PER UNO SCENARIO
# Uno scenario è identificato da:
# - ampiezza campionaria;
# - parametri AR(1)-GARCH(1,1): mu, phi, omega, alpha, beta;
# - tipo di innovazione non gaussiana;
# - gradi di libertà dell'innovazione.
# Per ogni replica Monte Carlo:
# 1) simula una serie AR(1)-GARCH(1,1);
# 2) costruisce gli intervalli di confidenza per media, SQM e VaR;
# 3) costruisce sia gli Intervalli di Confidenza corretti con Newey-West sia gli Intervalli di Confidenza iid;
# 4) verifica se ciascun Intervallo di Confidenza contiene il valore vero del parametro;
# 5) salva l'ampiezza di ciascun Intervallo di Confidenza.
# Alla fine restituisce:
# - copertura effettiva degli Intervalli di Confidenza;
# - ampiezza media degli Intervalli di Confidenza.
esperimento_montecarlo_uno_scenario <- function(numero_simulazioni_montecarlo, ampiezza_campionaria, mu, phi, omega, alpha, beta, innovazione, gradi_di_liberta, da_escludere, livello_quantile_p, livello_significativita, numero_lag_NeweyWest_m, numero_osservazioni_per_VaR_vero, da_escludere_per_VaR_vero) {
  # Per valutare la copertura devo sapere quali sono i veri valori di: media, scarto quadratico medio, VaR.
  valori_veri <- parametri_veri_ar1_garch11(mu = mu, phi = phi, omega = omega, alpha = alpha, beta = beta, innovazione = innovazione, gradi_di_liberta = gradi_di_liberta, livello_quantile_p = livello_quantile_p, numero_osservazioni_per_VaR_vero = numero_osservazioni_per_VaR_vero, da_escludere_per_VaR_vero = da_escludere_per_VaR_vero)
  
  # Vettore in cui salvare i risultati di MonteCarlo (lunghezza pari al numero di simulazioni).
  # In copertura (frequenza con cui il valore vero del parametro appartiene all'intervallo di confidenza) salvo:
  # - 1 se l'intervallo contiene il valore vero;
  # - 0 se l'intervallo non contiene il valore vero.
  # In ampiezza salvo:
  # - ampiezza = limite superiore - limite inferiore.
  copertura_media_Newey_West <- numeric(numero_simulazioni_montecarlo)
  copertura_media_iid <- numeric(numero_simulazioni_montecarlo)
  copertura_scarto_quadratico_medio_Newey_West <- numeric(numero_simulazioni_montecarlo)
  copertura_scarto_quadratico_medio_iid <- numeric(numero_simulazioni_montecarlo)
  copertura_VaR_Newey_West <- numeric(numero_simulazioni_montecarlo)
  copertura_VaR_iid <- numeric(numero_simulazioni_montecarlo)
  ampiezza_media_Newey_West <- numeric(numero_simulazioni_montecarlo)
  ampiezza_media_iid <- numeric(numero_simulazioni_montecarlo)
  ampiezza_scarto_quadratico_medio_Newey_West <- numeric(numero_simulazioni_montecarlo)
  ampiezza_scarto_quadratico_medio_iid <- numeric(numero_simulazioni_montecarlo)
  ampiezza_VaR_Newey_West <- numeric(numero_simulazioni_montecarlo)
  ampiezza_VaR_iid <- numeric(numero_simulazioni_montecarlo)
  
  # Ripeto l'esperimento numero_simulazioni_montecarlo volte. Ad ogni iterazione b:
  # - simulo una nuova serie;
  # - calcolo gli intervalli;
  # - verifico la copertura;
  # - salvo l'ampiezza.
  for (b in 1:numero_simulazioni_montecarlo) {
    # Simulazione di una serie AR(1)-GARCH(1,1)
    serie_simulata <- simula_ar1_garch11(n = ampiezza_campionaria, mu = mu, phi = phi, omega = omega, alpha = alpha, beta = beta, innovazione = innovazione, gradi_di_liberta = gradi_di_liberta, da_escludere = da_escludere)
    
    # Intervalli di confidenza per la media
    # - uno correggendo la dipendenza seriale con Newey-West;
    # - uno ignorando erroneamente la dipendenza seriale trattando la serie come iid.
    intervallo_media_Newey_West <- intervallo_confidenza_media_simulazione(serie_osservata_x = serie_simulata, correggi_dipendenza = TRUE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    intervallo_media_iid <- intervallo_confidenza_media_simulazione(serie_osservata_x = serie_simulata, correggi_dipendenza = FALSE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    
    # Intervalli di confidenza per lo scarto quadratico medio
    # - uno correggendo la dipendenza seriale con Newey-West;
    # - uno ignorando erroneamente la dipendenza seriale trattando la serie come iid.
    intervallo_scarto_quadratico_medio_Newey_West <- intervallo_confidenza_SQM_simulazione(serie_osservata_x = serie_simulata, correggi_dipendenza = TRUE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    intervallo_scarto_quadratico_medio_iid <- intervallo_confidenza_SQM_simulazione(serie_osservata_x = serie_simulata, correggi_dipendenza = FALSE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    
    # Intervalli di confidenza per il VaR
    # - uno correggendo la dipendenza seriale con Newey-West;
    # - uno ignorando erroneamente la dipendenza seriale trattando la serie come iid.
    intervallo_VaR_Newey_West <- intervallo_confidenza_VaR_simulazione(serie_osservata_x = serie_simulata, livello_quantile_p = livello_quantile_p, correggi_dipendenza = TRUE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    intervallo_VaR_iid <- intervallo_confidenza_VaR_simulazione(serie_osservata_x = serie_simulata, livello_quantile_p = livello_quantile_p, correggi_dipendenza = FALSE, livello_significativita = livello_significativita, numero_lag_NeweyWest_m = numero_lag_NeweyWest_m)
    
    # Verifica copertura degli intervalli
    # Un intervallo copre il valore vero se il valore_vero è compreso tra limite_inferiore e limite_superiore
    copertura_media_Newey_West[b] <- (intervallo_media_Newey_West["limite_inferiore"] <= valori_veri$media_vera && valori_veri$media_vera <= intervallo_media_Newey_West["limite_superiore"])
    copertura_media_iid[b] <- (intervallo_media_iid["limite_inferiore"] <= valori_veri$media_vera && valori_veri$media_vera <= intervallo_media_iid["limite_superiore"])
    copertura_scarto_quadratico_medio_Newey_West[b] <- (intervallo_scarto_quadratico_medio_Newey_West["limite_inferiore"] <= valori_veri$scarto_quadratico_medio_vero && valori_veri$scarto_quadratico_medio_vero <= intervallo_scarto_quadratico_medio_Newey_West["limite_superiore"])
    copertura_scarto_quadratico_medio_iid[b] <- (intervallo_scarto_quadratico_medio_iid["limite_inferiore"] <= valori_veri$scarto_quadratico_medio_vero && valori_veri$scarto_quadratico_medio_vero <= intervallo_scarto_quadratico_medio_iid["limite_superiore"])
    copertura_VaR_Newey_West[b] <- (intervallo_VaR_Newey_West["limite_inferiore"] <= valori_veri$VaR_vero && valori_veri$VaR_vero <= intervallo_VaR_Newey_West["limite_superiore"])
    copertura_VaR_iid[b] <- (intervallo_VaR_iid["limite_inferiore"] <= valori_veri$VaR_vero && valori_veri$VaR_vero <= intervallo_VaR_iid["limite_superiore"])
    
    # Salvataggio dell'ampiezza degli intervalli
    ampiezza_media_Newey_West[b] <- intervallo_media_Newey_West["ampiezza"]
    ampiezza_media_iid[b] <- intervallo_media_iid["ampiezza"]
    ampiezza_scarto_quadratico_medio_Newey_West[b] <- intervallo_scarto_quadratico_medio_Newey_West["ampiezza"]
    ampiezza_scarto_quadratico_medio_iid[b] <- intervallo_scarto_quadratico_medio_iid["ampiezza"]
    ampiezza_VaR_Newey_West[b] <- intervallo_VaR_Newey_West["ampiezza"]
    ampiezza_VaR_iid[b] <- intervallo_VaR_iid["ampiezza"]
  }
  # Sintesi dei risultati
  # La copertura effettiva è la media degli indicatori di copertura.
  # L'ampiezza media è la media delle ampiezze ottenute nelle diverse repliche.
  risultati_scenario <- data.frame(
    # Parametri dello scenario simulato
    numero_simulazioni_montecarlo = numero_simulazioni_montecarlo,
    ampiezza_campionaria = ampiezza_campionaria,
    mu = mu,
    phi = phi,
    omega = omega,
    alpha = alpha,
    beta = beta,
    innovazione = innovazione,
    gradi_di_liberta = gradi_di_liberta,
    livello_quantile_p = livello_quantile_p,
    livello_significativita = livello_significativita,
    # Valori veri (approssimati come veri) dei parametri stimati
    media_vera = valori_veri$media_vera,
    scarto_quadratico_medio_vero = valori_veri$scarto_quadratico_medio_vero,
    VaR_vero = valori_veri$VaR_vero,
    # Risultati per la media
    copertura_media_Newey_West = mean(copertura_media_Newey_West),
    copertura_media_iid = mean(copertura_media_iid),
    ampiezza_media_Newey_West = mean(ampiezza_media_Newey_West),
    ampiezza_media_iid = mean(ampiezza_media_iid),
    # Risultati per lo scarto quadratico medio
    copertura_scarto_quadratico_medio_Newey_West = mean(copertura_scarto_quadratico_medio_Newey_West),
    copertura_scarto_quadratico_medio_iid = mean(copertura_scarto_quadratico_medio_iid),
    ampiezza_scarto_quadratico_medio_Newey_West = mean(ampiezza_scarto_quadratico_medio_Newey_West),
    ampiezza_scarto_quadratico_medio_iid = mean(ampiezza_scarto_quadratico_medio_iid),
    # Risultati per il VaR
    copertura_VaR_Newey_West = mean(copertura_VaR_Newey_West),
    copertura_VaR_iid = mean(copertura_VaR_iid),
    ampiezza_VaR_Newey_West = mean(ampiezza_VaR_Newey_West),
    ampiezza_VaR_iid = mean(ampiezza_VaR_iid)
  )
  # La funzione restituisce una tabella con una sola riga.
  # Ogni riga corrisponde a uno scenario Monte Carlo.
  return(risultati_scenario)
}


# ABBREVIAZIONE DELLE ETICHETTE DEGLI SCENARI PER GRAFICO
abbrevia_etichette_scenario <- function(etichette) {
  # gsub(pattern, replacement, x) = sostituisce una stringa di testo con un'altra.
  # - pattern = testo da cercare.
  # - replacement = testo da inserire al posto del testo cercato.
  # - x = vettore di stringhe su cui effettuare la sostituzione.
  etichette <- gsub("Scenario base Amazon", "Base Amazon", etichette)
  etichette <- gsub("Ampiezza campionaria: n = ", "n = ", etichette)
  etichette <- gsub("Innovazione t-Student", "t base", etichette)
  etichette <- gsub("Innovazione t-Student code molto pesanti", "t code molto pes.", etichette)
  etichette <- gsub("Innovazione t-Student code poco pesanti", "t code poco pes.", etichette)
  etichette <- gsub("Innovazione t-Student vicino Normale", "t vicino Normale", etichette)
  etichette <- gsub("Innovazione Chi-Quadro", "Chi2 base", etichette)
  etichette <- gsub("Innovazione Chi-Quadro alta asimmetria", "Chi2 alta asim.", etichette)
  etichette <- gsub("Innovazione Chi-Quadro bassa asimmetria", "Chi2 bassa asim.", etichette)
  etichette <- gsub("Innovazione Chi-Quadro vicino Normale", "Chi2 vicino Normale", etichette)
  etichette <- gsub("GARCH persistenza bassa", "GARCH pers. bassa", etichette)
  etichette <- gsub("GARCH persistenza moderata", "GARCH pers. mod.", etichette)
  etichette <- gsub("GARCH bassa reattività shock", "GARCH bassa reatt.", etichette)
  etichette <- gsub("GARCH alta reattività shock", "GARCH alta reatt.", etichette)
  return(etichette)
}


# COSTRUZIONE GRAFICO A BARRE ORIZZONTALI
grafico_confronto_NeweyWest_iid <- function(tabella_grafico, parametro_scelto, titolo_grafico, etichetta_asse_x, aggiungi_linea_copertura = FALSE) {
  # Seleziono solo le righe relative al parametro scelto
  tabella_parametro <- tabella_grafico[tabella_grafico$parametro == parametro_scelto, ]
  # rbind = combina vettori o tabelle per riga (formato richiesto da barplot per disegnare barre affiancate).
  # Qui costruisco una matrice con una riga per i risultati Newey-West e una riga per i risultati i.i.d.
  matrice_barre <- rbind("Newey-West" = tabella_parametro$Newey_West, "i.i.d." = tabella_parametro$iid)
  # Se sto rappresentando la copertura fisso l'asse x tra 0 e 1 dato che è una frequenza. 
  if (aggiungi_linea_copertura == TRUE) {
    limite_asse_x <- 1
  } 
  # Se sto rappresentando l'ampiezza scelgo automaticamente il limite massimo dell'asse x.
  else {
    limite_asse_x <- max(matrice_barre, na.rm = TRUE) * 1.20
  }
  # layout = divide la finestra grafica in più aree.
  # Divido la finestra grafica in due parti:
  layout(matrix(c(1, 2), nrow = 2), heights = c(12, 1.4))
  # par(mar = ...) = imposta i margini del grafico: basso, sinistra, alto, destra.
  # Parte superiore: grafico
  par(mar = c(5, 14, 4, 2) + 0.1)
  # barplot = costruisce un grafico a barre.
  # - height = valori da rappresentare.
  # - beside = TRUE disegna le barre Newey-West e i.i.d. affiancate per ogni scenario.
  # - horiz = TRUE rende il grafico orizzontale più leggibile con molte etichette.
  # - names.arg = nomi degli scenari mostrati sull'asse verticale.
  barplot(
    height = matrice_barre,
    beside = TRUE,
    horiz = TRUE,
    names.arg = tabella_parametro$scenario_breve,
    las = 1,
    cex.names = 0.70,
    xlim = c(0, limite_asse_x),
    main = titolo_grafico,
    xlab = etichetta_asse_x,
    col = c("gray", "white"),
    border = "black"
  )
  if (aggiungi_linea_copertura == TRUE) {
    # abline = aggiunge una linea al grafico. Aggiungo la linea verticale della copertura teorica pari a 0.95.
    abline(v = 0.95, lty = 2, lwd = 1)
  }
  # Parte inferiore: legenda
  par(mar = c(0, 0, 0, 0))
  plot.new()
  if (aggiungi_linea_copertura == TRUE) {
    # legend = aggiunge la legenda al grafico.
    # - fill indica il colore dei simboli pieni; 
    # - lty e lwd servono per rappresentare la linea tratteggiata.
    # - horiz = TRUE dispone gli elementi della legenda in orizzontale.
    legend(
      "center",
      legend = c("= Newey-West", "= i.i.d.", "= Copertura teorica 0.95"),
      fill = c("gray", "white", NA),
      border = c("black", "black", NA),
      lty = c(NA, NA, 2),
      lwd = c(NA, NA, 1),
      horiz = TRUE,
      cex = 0.90
    )
  } else {
    legend(
      "center",
      legend = c("= Newey-West", "= i.i.d."),
      fill = c("gray", "white"),
      border = "black",
      horiz = TRUE,
      cex = 0.90
    )
  }
  # Ripristino il layout standard per i grafici successivi, in questo modo i grafici successivi non ereditano la divisione in due parti.
  layout(1)
}



#########################################################################################
# Caricamento dei dati azionari Amazon
# read_excel = legge un file Excel; sheet = indica il foglio del file da importare.
dati_azionari_Amazon <- read_excel("data/Price_Amazon.xlsx", sheet = "AMZN-US")
# Ordino in ordine temporale crescente (i dati sono ordinati in ordine temporale decrescente)
# as.Date = converte la variabile Date in formato data.
dati_azionari_Amazon$Date <- as.Date(dati_azionari_Amazon$Date)
# order = restituisce l'ordine degli indici che dispone le date in ordine crescente.
dati_azionari_Amazon <- dati_azionari_Amazon[order(dati_azionari_Amazon$Date), ]
# Estraggo i prezzi e li converto in formato numerico.
prezzi <- as.numeric(dati_azionari_Amazon$Price)
# is.na = individua i valori mancanti; !is.na seleziona solo i valori non mancanti.
prezzi <- prezzi[!is.na(prezzi)]
# Calcolo i rendimenti logaritmici.
rendimenti_logaritmici <- diff(log(prezzi))
rendimenti_logaritmici <- rendimenti_logaritmici[!is.na(rendimenti_logaritmici)]

#########################################################################################
# PUNTO 2 - STATISTICHE DESCRITTIVE e INTERVALLI DI CONFIDENZA SULLA SERIE STORICA AMAZON
#########################################################################################
#----------------------------------------------------------------------------------------
# 2A) Statistiche descrittive di base della serie storica.
#----------------------------------------------------------------------------------------
statistiche_descrittive <- list(
  numero_osservazioni_n = length(rendimenti_logaritmici),
  percentuale_rendimenti_positivi = mean(rendimenti_logaritmici >= 0),
  percentuale_rendimenti_negativi = mean(rendimenti_logaritmici < 0),
  media = mean(rendimenti_logaritmici),
  mediana = median(rendimenti_logaritmici),
  minimo = min(rendimenti_logaritmici),
  massimo = max(rendimenti_logaritmici),
  varianza = mean((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2),
  varianza_corretta = sum((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2)/(length(rendimenti_logaritmici)-1),
  scarto_quadratico_medio = sqrt(mean((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2)), 
  # quantile(x, probs) = calcola il quantile campionario della variabile x al livello indicato da probs.
  # - x = serie numerica su cui calcolare il quantile.
  # - probs = probabilità/livello del quantile, compreso tra 0 e 1.
  primo_quartile = quantile(rendimenti_logaritmici, 0.25),  
  terzo_quartile = quantile(rendimenti_logaritmici, 0.75),  
  VaR_5_rendimenti = quantile(rendimenti_logaritmici, 0.05),           
  asimmetria = mean(((rendimenti_logaritmici - mean(rendimenti_logaritmici)) / sqrt(mean((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2)))^3),
  curtosi = mean(((rendimenti_logaritmici - mean(rendimenti_logaritmici)) / sqrt(mean((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2)))^4),
  eccesso_curtosi = (mean(((rendimenti_logaritmici - mean(rendimenti_logaritmici)) / sqrt(mean((rendimenti_logaritmici - mean(rendimenti_logaritmici))^2)))^4))-3
)
statistiche_descrittive

if (!dir.exists("results/inference")) {dir.create("results/inference", recursive = TRUE)}
# Grafici diagnostici sui rendimenti storici.
x11()
plot(rendimenti_logaritmici, type = "l", main = "Rendimenti logaritmici Amazon", ylab = "rendimenti_t")
x11()
hist(rendimenti_logaritmici, breaks = 40, probability = TRUE, main = "Istogramma rendimenti", xlab = "rendimenti_t")
x11()
acf(rendimenti_logaritmici, main = "ACF rendimenti")
x11()
acf(rendimenti_logaritmici^2, main = "ACF rendimenti^2")
png(filename = "results/inference/Rendimenti_logaritmici_Amazon.png", width = 1400, height = 900)
plot(rendimenti_logaritmici, type = "l", main = "Rendimenti logaritmici Amazon", ylab = "rendimenti_t")
dev.off()
png(filename = "results/inference/Istogramma_rendimenti_Amazon.png", width = 1400, height = 900)
hist(rendimenti_logaritmici, breaks = 40, probability = TRUE, main = "Istogramma rendimenti", xlab = "rendimenti_t")
dev.off()
png(filename = "results/inference/ACF_rendimenti_Amazon.png", width = 1400, height = 900)
acf(rendimenti_logaritmici, main = "ACF rendimenti")
dev.off()
png(filename = "results/inference/ACF_rendimenti_quadrati_Amazon.png", width = 1400, height = 900)
acf(rendimenti_logaritmici^2, main = "ACF rendimenti^2")
dev.off()

#----------------------------------------------------------------------------------------
# 2B) Intervalli di confidenza
#----------------------------------------------------------------------------------------
# 2a) Intervallo di confidenza per il rendimento atteso
intervallo_media_Amazon_Newey_West <- intervallo_confidenza_media_simulazione(
  serie_osservata_x = rendimenti_logaritmici,
  correggi_dipendenza = TRUE,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL
)
intervallo_media_Amazon_Newey_West

rendimento_atteso_diverso_da_zero <- !(intervallo_media_Amazon_Newey_West["limite_inferiore"] <= 0 && 0 <= intervallo_media_Amazon_Newey_West["limite_superiore"])
# Il rendimento atteso è significativamente diverso da 0 (0 NON appartiene all'intervallo di confidenza bilaterale)?
rendimento_atteso_diverso_da_zero

varianza_asintotica_media_Amazon <- NeweyWest_varianza_lungo_periodo(serie_osservata_x = rendimenti_logaritmici, numero_lag_NeweyWest_m = NULL)
standard_error_media_Amazon <- sqrt(varianza_asintotica_media_Amazon / length(rendimenti_logaritmici))
limite_inferiore_unilaterale_media_Amazon <- statistiche_descrittive$media - qnorm(1 - 0.05) * standard_error_media_Amazon
rendimento_atteso_maggiore_di_zero <- limite_inferiore_unilaterale_media_Amazon > 0
# Il rendimento atteso è significativamente maggiore di 0 (H0: mu <= 0 contro H1: mu > 0, intervallo unilaterale inferiore: media - z_{1-alpha} * SE)?
rendimento_atteso_maggiore_di_zero

# 2b) Intervallo di confidenza per lo scarto quadratico medio
intervallo_scarto_quadratico_medio_Amazon_Newey_West <- intervallo_confidenza_SQM_simulazione(
  serie_osservata_x = rendimenti_logaritmici,
  correggi_dipendenza = TRUE,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL
)
intervallo_scarto_quadratico_medio_Amazon_Newey_West

# 2c) Intervallo di confidenza per il VaR
intervallo_VaR_Amazon_Newey_West <- intervallo_confidenza_VaR_simulazione(
  serie_osservata_x = rendimenti_logaritmici,
  livello_quantile_p = 0.05,
  correggi_dipendenza = TRUE,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL
)
intervallo_VaR_Amazon_Newey_West

#########################################################################################
# PUNTO 3 
#########################################################################################
# Stimo un modello AR(1)-GARCH(1,1) con innovazioni t-Student (AR(1) è come scrivere ARMA(1,0)).
# garchFit = funzione del pacchetto fGarch usata per stimare modelli GARCH.
# formula = specifica la dinamica del modello:
# - arma(1, 0) indica una componente autoregressiva AR(1) nella media;
# - garch(1, 1) indica una componente GARCH(1,1) nella varianza condizionata.
# data = serie storica su cui stimare il modello.
# cond.dist = "std" impone innovazioni t-Student standardizzate nella stima storica del modello base. Questa scelta serve solo per stimare i parametri empirici iniziali; negli scenari Monte Carlo il tipo di innovazione viene poi scelto separatamente tramite innovazione = "t" oppure innovazione = "chisq".
# trace = FALSE evita di stampare a video i dettagli iterativi della procedura di stima.
stima_storica_AR1_GARCH11 <- garchFit(
  formula = ~ arma(1, 0) + garch(1, 1),
  data = rendimenti_logaritmici,
  cond.dist = "std",
  trace = FALSE
)
# Estraggo i coefficienti stimati.
# coef = estrae dal modello stimato i coefficienti numerici stimati.
coefficienti_stimati <- coef(stima_storica_AR1_GARCH11)

# Assegno i nomi che voglio.
variabili <- c("mu", "ar1", "omega", "alpha1", "beta1", "shape")
mu_stimato <- as.numeric(coefficienti_stimati["mu"])
phi_stimato <- as.numeric(coefficienti_stimati["ar1"])
omega_stimato <- as.numeric(coefficienti_stimati["omega"])
alpha_stimato <- as.numeric(coefficienti_stimati["alpha1"])
beta_stimato <- as.numeric(coefficienti_stimati["beta1"])
gradi_liberta_t_student_stimati <- as.numeric(coefficienti_stimati["shape"])
cat("\nParametri stimati del modello AR(1)-GARCH(1,1) scenario base per la simulazione:\n")
cat("- mu stimato: intercetta/media condizionata del processo AR(1) =", mu_stimato, "\n")
cat("- phi stimato: coefficiente autoregressivo AR(1) (peso del rendimento precedente) =", phi_stimato, "\n")
cat("- omega stimato: termine costante dell'equazione GARCH(1,1) della varianza condizionata =", omega_stimato, "\n")
cat("- alpha stimato: coefficiente ARCH (peso dello shock quadratico precedente nella varianza condizionata) =", alpha_stimato, "\n")
cat("- beta stimato: coefficiente GARCH (peso della varianza condizionata precedente) =", beta_stimato, "\n")
cat("- gradi di libertà t-Student stimati: parametro della distribuzione non gaussiana delle innovazioni =", gradi_liberta_t_student_stimati, "\n")

#----------------------------------------------------------------------------------------
# 3A) ESPERIMENTO MONTE CARLO: SCENARIO BASE AMAZON
#----------------------------------------------------------------------------------------
# Uso i parametri stimati sulla serie Amazon come scenario base mantenendo la struttura empirica stimata sui rendimenti Amazon.
set.seed(123)
risultati_scenario_base_Amazon <- esperimento_montecarlo_uno_scenario(
  # 1000 repliche sono un buon compromesso tra stabilità dei risultati e tempo di calcolo.
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)
risultati_scenario_base_Amazon$scenario <- "Scenario base Amazon"
risultati_scenario_base_Amazon <- risultati_scenario_base_Amazon[ , c("scenario", setdiff(names(risultati_scenario_base_Amazon), "scenario"))]
risultati_scenario_base_Amazon

#----------------------------------------------------------------------------------------
# 3B) IMPATTO DELL'AMPIEZZA CAMPIONARIA
#----------------------------------------------------------------------------------------
risultati_ampiezza_250 <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = 250,
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_ampiezza_500 <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = 500,
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_ampiezza_1000 <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = 1000,
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_ampiezza_2000 <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = 2000,
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)
risultati_ampiezza_250$scenario <- "Ampiezza campionaria: n = 250"
risultati_ampiezza_500$scenario <- "Ampiezza campionaria: n = 500"
risultati_ampiezza_1000$scenario <- "Ampiezza campionaria: n = 1000"
risultati_ampiezza_2000$scenario <- "Ampiezza campionaria: n = 2000"
# rbind = combina vettori o tabelle per riga.
risultati_impatto_ampiezza_campionaria <- rbind(risultati_ampiezza_250, risultati_ampiezza_500, risultati_ampiezza_1000, risultati_ampiezza_2000)
# Riordino le colonne della tabella mettendo "scenario" come prima colonna.
# names(...) restituisce i nomi di tutte le colonne della tabella.
# setdiff(names(...), "scenario") restituisce tutti i nomi di colonna tranne "scenario".
# c("scenario", ...) ricostruisce l'ordine desiderato delle colonne.
# La virgola prima di c(...) indica che sto selezionando tutte le righe e solo le colonne nell'ordine specificato.
risultati_impatto_ampiezza_campionaria <- risultati_impatto_ampiezza_campionaria[ , c("scenario", setdiff(names(risultati_impatto_ampiezza_campionaria), "scenario"))]
risultati_impatto_ampiezza_campionaria

#----------------------------------------------------------------------------------------
# 3C) IMPATTO DEL TIPO E DELLA FORMA DISTRIBUTIVA DELLE INNOVAZIONI
#----------------------------------------------------------------------------------------
risultati_innovazione_tStudent_code_molto_pesanti <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = max(min(gradi_liberta_t_student_stimati / 2, 2.5), 2.05),
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_chi_quadro_alta_asimmetria <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "chisq",
  gradi_di_liberta = max(min(gradi_liberta_t_student_stimati / 2, 2.5), 2.05),
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_tStudent <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_chi_quadro <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "chisq",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_tStudent_code_meno_pesanti <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati*2,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_chi_quadro_bassa_asimmetria <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "chisq",
  gradi_di_liberta = gradi_liberta_t_student_stimati*2,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_tStudent_vicino_normale <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = 20,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

risultati_innovazione_chi_quadro_vicino_normale <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "chisq",
  gradi_di_liberta = 20,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)
risultati_innovazione_tStudent_code_molto_pesanti$scenario <- "Innovazione t-Student code molto pesanti"
risultati_innovazione_chi_quadro_alta_asimmetria$scenario <- "Innovazione Chi-Quadro alta asimmetria"
risultati_innovazione_tStudent$scenario <- "Innovazione t-Student"
risultati_innovazione_chi_quadro$scenario <- "Innovazione Chi-Quadro"
risultati_innovazione_tStudent_code_meno_pesanti$scenario <- "Innovazione t-Student code poco pesanti"
risultati_innovazione_chi_quadro_bassa_asimmetria$scenario <- "Innovazione Chi-Quadro bassa asimmetria"
risultati_innovazione_tStudent_vicino_normale$scenario <- "Innovazione t-Student vicino Normale"
risultati_innovazione_chi_quadro_vicino_normale$scenario <- "Innovazione Chi-Quadro vicino Normale"
# rbind = combina vettori o tabelle per riga.
risultati_impatto_tipo_innovazione <- rbind(risultati_innovazione_tStudent_code_molto_pesanti, risultati_innovazione_chi_quadro_alta_asimmetria , risultati_innovazione_tStudent, risultati_innovazione_chi_quadro, risultati_innovazione_tStudent_code_meno_pesanti , risultati_innovazione_chi_quadro_bassa_asimmetria, risultati_innovazione_tStudent_vicino_normale, risultati_innovazione_chi_quadro_vicino_normale)
# Riordino le colonne della tabella mettendo "scenario" come prima colonna.
# names(...) restituisce i nomi di tutte le colonne della tabella.
# setdiff(names(...), "scenario") restituisce tutti i nomi di colonna tranne "scenario".
# c("scenario", ...) ricostruisce l'ordine desiderato delle colonne.
# La virgola prima di c(...) indica che sto selezionando tutte le righe e solo le colonne nell'ordine specificato.
risultati_impatto_tipo_innovazione <- risultati_impatto_tipo_innovazione[ , c("scenario", setdiff(names(risultati_impatto_tipo_innovazione), "scenario"))]
risultati_impatto_tipo_innovazione

#----------------------------------------------------------------------------------------
# 3D) IMPATTO DEI PARAMETRI DEL PROCESSO
#----------------------------------------------------------------------------------------
# Scenario con media condizionata più bassa
risultati_mu_basso <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato/2,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con media condizionata più alta
risultati_mu_alto <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato * 2,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con dipendenza moderata nella media.
risultati_phi_moderato  <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = 0.10,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con alta dipendenza nella media.
risultati_phi_alto <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = 0.60,
  omega = omega_stimato,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con livello medio di volatilità più basso.
risultati_omega_basso <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato/2,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con livello medio di volatilità più alto.
risultati_omega_alto <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato * 2,
  alpha = alpha_stimato,
  beta = beta_stimato,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

#----------------------------------------------------------------------------------------
# Scenari GARCH: persistenza complessiva e composizione alpha/beta
# Nel modello GARCH(1,1):
# - alpha + beta misura la persistenza complessiva della volatilità;
# - alpha misura la reattività della volatilità agli shock recenti;
# - beta misura la persistenza lenta della volatilità passata.
# Il caso base Amazon ha già alpha + beta molto vicino a 1, quindi rappresenta già uno scenario di alta persistenza GARCH.
# Per confrontare gli scenari in modo pulito, calcolo la varianza incondizionata dello shock nello scenario base:
#     Var(shock_t) = omega / (1 - alpha - beta)
# e, quando modifico alpha e beta, ricalcolo omega in modo da mantenere costante questa varianza incondizionata.
#----------------------------------------------------------------------------------------

# Scenario con persistenza GARCH bassa (alpha + beta = 0.75)
# Mantengo la stessa proporzione alpha/(alpha+beta) osservata nel caso base.
risultati_garch_bassa_persistenza <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = (omega_stimato / (1 - (alpha_stimato + beta_stimato))) * (1 - 0.75),
  alpha = (alpha_stimato / (alpha_stimato + beta_stimato)) * 0.75,
  beta = (1 - (alpha_stimato / (alpha_stimato + beta_stimato))) * 0.75,
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con persistenza GARCH moderata (alpha + beta = 0.90)
# Mantengo la stessa proporzione alpha/(alpha+beta) osservata nel caso base.
risultati_garch_moderata_persistenza <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = ((omega_stimato / (1 - (alpha_stimato + beta_stimato))) * (1 - 0.90)),
  alpha = ((alpha_stimato / (alpha_stimato + beta_stimato)) * 0.90),
  beta = ((1 - (alpha_stimato / (alpha_stimato + beta_stimato))) * 0.90),
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con bassa reattività agli shock recenti.
# Mantengo costante alpha + beta al caso base (se riduco alpha allora aumento beta).
# In questo modo la persistenza complessiva resta invariata, ma do meno peso agli shock recenti e do più peso alla volatilità passata.
if (((alpha_stimato + beta_stimato) - 0.01) <= 0) {
  stop("Scenario non valido: beta <= 0.")
}
risultati_garch_bassa_reattivita <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = 0.01,
  beta = ((alpha_stimato + beta_stimato) - 0.01),
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)

# Scenario con alta reattività agli shock recenti.
# Mantengo costante alpha + beta al caso base (se aumento alpha allora riduco beta).
# In questo modo la persistenza complessiva resta invariata, ma la volatilità reagisce più intensamente agli shock recenti.
if ((min(0.10, (alpha_stimato + beta_stimato) - 0.01)) <= 0 || ((alpha_stimato + beta_stimato) - (min(0.10, (alpha_stimato + beta_stimato) - 0.01))) <= 0) {
  stop("Scenario non valido: alpha o beta non sono positivi.")
}
risultati_garch_alta_reattivita <- esperimento_montecarlo_uno_scenario(
  numero_simulazioni_montecarlo = 1000,
  ampiezza_campionaria = length(rendimenti_logaritmici),
  mu = mu_stimato,
  phi = phi_stimato,
  omega = omega_stimato,
  alpha = (min(0.10, (alpha_stimato + beta_stimato) - 0.01)),
  beta = ((alpha_stimato + beta_stimato) - (min(0.10, (alpha_stimato + beta_stimato) - 0.01))),
  innovazione = "t",
  gradi_di_liberta = gradi_liberta_t_student_stimati,
  da_escludere = 100,
  livello_quantile_p = 0.05,
  livello_significativita = 0.05,
  numero_lag_NeweyWest_m = NULL,
  numero_osservazioni_per_VaR_vero = 100000,
  da_escludere_per_VaR_vero = 1000
)
risultati_mu_basso$scenario <- "Mu basso"
risultati_mu_alto$scenario <- "Mu alto"
risultati_phi_moderato$scenario <- "Phi moderato"
risultati_phi_alto$scenario <- "Phi alto"
risultati_omega_basso$scenario <- "Omega basso"
risultati_omega_alto$scenario <- "Omega alto"
risultati_garch_bassa_persistenza$scenario <- "GARCH persistenza bassa"
risultati_garch_moderata_persistenza$scenario <- "GARCH persistenza moderata"
risultati_garch_bassa_reattivita$scenario <- "GARCH bassa reattività shock"
risultati_garch_alta_reattivita$scenario <- "GARCH alta reattività shock"
# rbind = combina vettori o tabelle per riga.
risultati_impatto_parametri_processo <- rbind(risultati_mu_basso, risultati_mu_alto, risultati_phi_moderato , risultati_phi_alto, risultati_omega_basso, risultati_omega_alto, risultati_garch_bassa_persistenza, risultati_garch_moderata_persistenza, risultati_garch_bassa_reattivita, risultati_garch_alta_reattivita)
# Riordino le colonne della tabella mettendo "scenario" come prima colonna.
# names(...) restituisce i nomi di tutte le colonne della tabella.
# setdiff(names(...), "scenario") restituisce tutti i nomi di colonna tranne "scenario".
# c("scenario", ...) ricostruisce l'ordine desiderato delle colonne.
# La virgola prima di c(...) indica che sto selezionando tutte le righe e solo le colonne nell'ordine specificato.
risultati_impatto_parametri_processo <- risultati_impatto_parametri_processo[ , c("scenario", setdiff(names(risultati_impatto_parametri_processo), "scenario"))]
risultati_impatto_parametri_processo

# TABELLA COMPLESSIVA DEI RISULTATI DEL PUNTO 3
# rbind = combina vettori o tabelle per riga.
risultati_punto_3_completi <- rbind(risultati_scenario_base_Amazon, risultati_impatto_ampiezza_campionaria, risultati_impatto_tipo_innovazione, risultati_impatto_parametri_processo)
# Riordino le colonne della tabella mettendo "scenario" come prima colonna.
# names(...) restituisce i nomi di tutte le colonne della tabella.
# setdiff(names(...), "scenario") restituisce tutti i nomi di colonna tranne "scenario".
# c("scenario", ...) ricostruisce l'ordine desiderato delle colonne.
# La virgola prima di c(...) indica che sto selezionando tutte le righe e solo le colonne nell'ordine specificato.
risultati_punto_3_completi <- risultati_punto_3_completi[ , c("scenario", setdiff(names(risultati_punto_3_completi), "scenario"))]
risultati_punto_3_completi

#----------------------------------------------------------------------------------------
# La tabella finale serve per confrontare:
# - copertura e ampiezza degli intervalli Newey-West;
# - copertura e ampiezza degli intervalli i.i.d.;
# - effetto dell'ampiezza campionaria;
# - effetto del tipo di innovazione;
# - effetto dei parametri AR(1)-GARCH(1,1).
#
# Aspettative:
# - aumentando l'ampiezza campionaria, gli intervalli diventano più stretti;
# - gli intervalli iid possono avere copertura peggiore quando la dipendenza seriale è rilevante;
# - innovazioni con code pesanti o asimmetriche rendono più difficile stimare correttamente il VaR;
# - maggiore dipendenza nella media, maggiore persistenza GARCH e maggiore reattività della volatilità agli shock recenti possono aumentare l'importanza della correzione Newey-West, soprattutto per SQM e VaR.
#----------------------------------------------------------------------------------------


#########################################################################################
# GRAFICI RIASSUNTIVI DEI RISULTATI MONTECARLO
# Serve per il confronto tra intervalli Newey-West e intervalli i.i.d. in:
# - copertura effettiva;
# - ampiezza media.
# Creo grafici separati per parametro:
# - Media
# - SQM
# - VaR
#########################################################################################
# Preparazione dati per grafici
coperture_montecarlo <- rbind(
  # data.frame = costruisce una tabella.
  # Qui creo  una per ciascun parametro: Media, SQM e VaR.
  # Ogni tabella contiene:
  # - scenario = nome dello scenario Monte Carlo;
  # - parametro = parametro a cui si riferisce la copertura;
  # - Newey_West = copertura effettiva ottenuta usando intervalli corretti con Newey-West;
  # - iid = copertura effettiva ottenuta usando intervalli che trattano erroneamente i dati come indipendenti.
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "Media",
    Newey_West = risultati_punto_3_completi$copertura_media_Newey_West,
    iid = risultati_punto_3_completi$copertura_media_iid
    ),
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "SQM",
    Newey_West = risultati_punto_3_completi$copertura_scarto_quadratico_medio_Newey_West,
    iid = risultati_punto_3_completi$copertura_scarto_quadratico_medio_iid
    ),
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "VaR",
    Newey_West = risultati_punto_3_completi$copertura_VaR_Newey_West,
    iid = risultati_punto_3_completi$copertura_VaR_iid
    )
)
ampiezze_montecarlo <- rbind(
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "Media",
    Newey_West = risultati_punto_3_completi$ampiezza_media_Newey_West,
    iid = risultati_punto_3_completi$ampiezza_media_iid
    ),
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "SQM",
    Newey_West = risultati_punto_3_completi$ampiezza_scarto_quadratico_medio_Newey_West,
    iid = risultati_punto_3_completi$ampiezza_scarto_quadratico_medio_iid
    ),
  data.frame(
    scenario = risultati_punto_3_completi$scenario,
    parametro = "VaR",
    Newey_West = risultati_punto_3_completi$ampiezza_VaR_Newey_West,
    iid = risultati_punto_3_completi$ampiezza_VaR_iid
    )
)
# Aggiungo le etichette abbreviate ai dataset grafici
coperture_montecarlo$scenario_breve <- abbrevia_etichette_scenario(coperture_montecarlo$scenario)
ampiezze_montecarlo$scenario_breve <- abbrevia_etichette_scenario(ampiezze_montecarlo$scenario)
#----------------------------------------------------------------------------------------
# GRAFICI DELLA COPERTURA EFFETTIVA
#----------------------------------------------------------------------------------------
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "Media", titolo_grafico = "Copertura effettiva - Media", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "SQM", titolo_grafico = "Copertura effettiva - Scarto Quadratico Medio", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "VaR", titolo_grafico = "Copertura effettiva - VaR", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
#----------------------------------------------------------------------------------------
# GRAFICI DELL'AMPIEZZA MEDIA
#----------------------------------------------------------------------------------------
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "Media", titolo_grafico = "Ampiezza media - Media", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "SQM", titolo_grafico = "Ampiezza media - Scarto Quadratico Medio", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
x11()
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "VaR", titolo_grafico = "Ampiezza media - VaR", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
png(filename = "results/inference/Copertura_effettiva_Media.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "Media", titolo_grafico = "Copertura effettiva - Media", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
dev.off()
png(filename = "results/inference/Copertura_effettiva_SQM.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "SQM", titolo_grafico = "Copertura effettiva - Scarto Quadratico Medio", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
dev.off()
png(filename = "results/inference/Copertura_effettiva_VaR.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = coperture_montecarlo, parametro_scelto = "VaR", titolo_grafico = "Copertura effettiva - VaR", etichetta_asse_x = "Copertura effettiva", aggiungi_linea_copertura = TRUE)
dev.off()
png(filename = "results/inference/Ampiezza_media_Media.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "Media", titolo_grafico = "Ampiezza media - Media", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
dev.off()
png(filename = "results/inference/Ampiezza_media_SQM.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "SQM", titolo_grafico = "Ampiezza media - Scarto Quadratico Medio", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
dev.off()
png(filename = "results/inference/Ampiezza_media_VaR.png", width = 1400, height = 1000)
grafico_confronto_NeweyWest_iid(tabella_grafico = ampiezze_montecarlo, parametro_scelto = "VaR", titolo_grafico = "Ampiezza media - VaR", etichetta_asse_x = "Ampiezza media", aggiungi_linea_copertura = FALSE)
dev.off()



# Riordiniamo i risultati in un file excel
library(openxlsx)
#########################################################################################
# ESPORTAZIONE DEI RISULTATI IN EXCEL
#########################################################################################
# Numero di osservazioni
numero_osservazioni_Amazon <- length(rendimenti_logaritmici)
# Lag Newey-West
numero_lag_NeweyWest_Amazon <- floor(4 * (numero_osservazioni_Amazon / 100)^(2 / 9))
numero_lag_NeweyWest_Amazon <- min(numero_lag_NeweyWest_Amazon, numero_osservazioni_Amazon - 1)
# Risultati Newey-West per SQM e VaR
risultato_NeweyWest_SQM_Amazon <- NeweyWest_scarto_quadratico_medio(serie_osservata_x = rendimenti_logaritmici, numero_lag_NeweyWest_m = NULL)
risultato_NeweyWest_VaR_Amazon <- NeweyWest_VaR(serie_osservata_x = rendimenti_logaritmici, livello_quantile_p = 0.05, numero_lag_NeweyWest_m = NULL)
# Statistica test e p-value per il rendimento atteso
t_stat_media_Amazon <- statistiche_descrittive$media / standard_error_media_Amazon
p_value_bilaterale_media_Amazon <- 2 * (1 - pnorm(abs(t_stat_media_Amazon)))
p_value_unilaterale_media_Amazon <- 1 - pnorm(t_stat_media_Amazon)
# Conclusioni test in formato testuale
conclusione_bilaterale_media_Amazon <- ifelse(
  rendimento_atteso_diverso_da_zero,
  "Rifiuto H0: il rendimento atteso è significativamente diverso da 0 al 5%",
  "Non rifiuto H0: il rendimento atteso non è significativamente diverso da 0 al 5%"
)
conclusione_unilaterale_media_Amazon <- ifelse(
  rendimento_atteso_maggiore_di_zero,
  "Rifiuto H0: il rendimento atteso è significativamente maggiore di 0 al 5%",
  "Non rifiuto H0: il rendimento atteso non è significativamente maggiore di 0 al 5%"
)
# Creo un nuovo file Excel.
file_risultati <- createWorkbook()
# Lista in cui salvo posizione e dimensioni delle tabelle scritte nei fogli Excel.
tabelle_da_formattare <- list()
# Funzione di supporto per registrare le dimensioni di ogni tabella esportata.
registra_tabella_excel <- function(nome_foglio, nome_tabella, tabella, riga_iniziale = 1, colonna_iniziale = 1) {
  tabelle_da_formattare[[nome_tabella]] <<- list(
    nome_foglio = nome_foglio,
    riga_iniziale = riga_iniziale,
    colonna_iniziale = colonna_iniziale,
    numero_righe = nrow(tabella) + 1,      # +1 perché considero anche la riga di intestazione
    numero_colonne = ncol(tabella)
  )
}
#----------------------------------------------------------------------------------------
# Foglio 1: Statistiche descrittive
addWorksheet(file_risultati, "Statistiche_descrittive")
tabella_statistiche_descrittive <- data.frame(statistica = names(statistiche_descrittive), valore = unlist(statistiche_descrittive))
writeData(wb = file_risultati, sheet = "Statistiche_descrittive", x = tabella_statistiche_descrittive)
registra_tabella_excel(nome_foglio = "Statistiche_descrittive", nome_tabella = "Statistiche_descrittive", tabella = tabella_statistiche_descrittive, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 2: Intervalli di confidenza e test sulla media
addWorksheet(file_risultati, "Intervallo_Confidenza")
tabella_intervalli_Amazon <- data.frame(
  parametro = c("Rendimento atteso", "Scarto Quadratico Medio", "VaR"),
  stima_puntuale = c(statistiche_descrittive$media, risultato_NeweyWest_SQM_Amazon$scarto_quadratico_medio_campionario, risultato_NeweyWest_VaR_Amazon$VaR_stimato),
  errore_standard_Newey_West = c(standard_error_media_Amazon, risultato_NeweyWest_SQM_Amazon$standard_error_stimato_SQM, risultato_NeweyWest_VaR_Amazon$standard_error_VaR),
  numero_lag_NeweyWest = c(numero_lag_NeweyWest_Amazon, risultato_NeweyWest_SQM_Amazon$numero_lag_NeweyWest_m, risultato_NeweyWest_VaR_Amazon$numero_lag_NeweyWest_m),
  limite_inferiore = c(intervallo_media_Amazon_Newey_West["limite_inferiore"], intervallo_scarto_quadratico_medio_Amazon_Newey_West["limite_inferiore"], intervallo_VaR_Amazon_Newey_West["limite_inferiore"]),
  limite_superiore = c(intervallo_media_Amazon_Newey_West["limite_superiore"], intervallo_scarto_quadratico_medio_Amazon_Newey_West["limite_superiore"], intervallo_VaR_Amazon_Newey_West["limite_superiore"]),
  ampiezza = c(intervallo_media_Amazon_Newey_West["ampiezza"], intervallo_scarto_quadratico_medio_Amazon_Newey_West["ampiezza"], intervallo_VaR_Amazon_Newey_West["ampiezza"])
)
writeData(wb = file_risultati, sheet = "Intervallo_Confidenza", x = tabella_intervalli_Amazon)
registra_tabella_excel(nome_foglio = "Intervallo_Confidenza", nome_tabella = "Intervallo_Confidenza", tabella = tabella_intervalli_Amazon, riga_iniziale = 1, colonna_iniziale = 1)
# Aggiungo anche gli esiti dei test sulla media.
tabella_test_media_Amazon <- data.frame(
  test = c(
    "Test bilaterale sulla media", 
    "Test unilaterale destro sulla media"
  ),
  ipotesi_nulla_H0 = c(
    "Rendimento atteso = 0", 
    "Rendimento atteso <= 0"
  ),
  ipotesi_alternativa_H1 = c(
    "Rendimento atteso diverso da 0", 
    "Rendimento atteso > 0"
  ),
  statistica_test_Newey_West = c(t_stat_media_Amazon, t_stat_media_Amazon),
  p_value = c(p_value_bilaterale_media_Amazon, p_value_unilaterale_media_Amazon),
  livello_significativita = c(0.05, 0.05),
  decisione = c(
    ifelse(rendimento_atteso_diverso_da_zero, "Rifiuto H0", "Non rifiuto H0"), 
    ifelse(rendimento_atteso_maggiore_di_zero, "Rifiuto H0", "Non rifiuto H0")
  ),
  conclusione = c(conclusione_bilaterale_media_Amazon, conclusione_unilaterale_media_Amazon)
)
writeData(wb = file_risultati, sheet = "Intervallo_Confidenza", x = tabella_test_media_Amazon, startRow = nrow(tabella_intervalli_Amazon) + 4)
registra_tabella_excel(nome_foglio = "Intervallo_Confidenza", nome_tabella = "Intervallo_Confidenza_Test_Media", tabella = tabella_test_media_Amazon, riga_iniziale = nrow(tabella_intervalli_Amazon) + 4, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 3: Scenario base Monte Carlo
addWorksheet(file_risultati, "MonteCarlo_Scenario_Base")
writeData(wb = file_risultati, sheet = "MonteCarlo_Scenario_Base", x = risultati_scenario_base_Amazon)
registra_tabella_excel(nome_foglio = "MonteCarlo_Scenario_Base", nome_tabella = "MonteCarlo_Scenario_Base", tabella = risultati_scenario_base_Amazon, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 4: Impatto ampiezza campionaria
addWorksheet(file_risultati, "MonteCarlo_Ampiezza_Campione")
writeData(wb = file_risultati, sheet = "MonteCarlo_Ampiezza_Campione", x = risultati_impatto_ampiezza_campionaria)
registra_tabella_excel(nome_foglio = "MonteCarlo_Ampiezza_Campione", nome_tabella = "MonteCarlo_Ampiezza_Campione", tabella = risultati_impatto_ampiezza_campionaria, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 5: Impatto tipo innovazione
addWorksheet(file_risultati, "MonteCarlo_Innovazioni")
writeData(wb = file_risultati, sheet = "MonteCarlo_Innovazioni", x = risultati_impatto_tipo_innovazione)
registra_tabella_excel(nome_foglio = "MonteCarlo_Innovazioni", nome_tabella = "MonteCarlo_Innovazioni", tabella = risultati_impatto_tipo_innovazione, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 6: Impatto parametri processo
addWorksheet(file_risultati, "MonteCarlo_Parametri")
writeData(wb = file_risultati, sheet = "MonteCarlo_Parametri", x = risultati_impatto_parametri_processo)
registra_tabella_excel(nome_foglio = "MonteCarlo_Parametri", nome_tabella = "MonteCarlo_Parametri", tabella = risultati_impatto_parametri_processo, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Foglio 7: Tutti i risultati del punto 3
addWorksheet(file_risultati, "MonteCarlo_Analisi_Completa")
writeData(wb = file_risultati, sheet = "MonteCarlo_Analisi_Completa", x = risultati_punto_3_completi)
registra_tabella_excel(nome_foglio = "MonteCarlo_Analisi_Completa", nome_tabella = "MonteCarlo_Analisi_Completa", tabella = risultati_punto_3_completi, riga_iniziale = 1, colonna_iniziale = 1)
#----------------------------------------------------------------------------------------
# Formattazione file
# Stile per l'intestazione
stile_intestazione <- createStyle(textDecoration = "bold", fgFill = "#D9EAF7", halign = "center", valign = "center", border = "Bottom", borderStyle = "thick")
# Stile per il corpo della tabella
stile_corpo_tabella <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
# Stili per il bordo esterno spesso della tabella.
stile_bordo_superiore_spesso <- createStyle(border = "Top", borderStyle = "thick")
stile_bordo_inferiore_spesso <- createStyle(border = "Bottom", borderStyle = "thick")
stile_bordo_sinistro_spesso <- createStyle(border = "Left", borderStyle = "thick")
stile_bordo_destro_spesso <- createStyle(border = "Right", borderStyle = "thick")
# Applico la formattazione a tutte le tabelle.
for (nome_tabella in names(tabelle_da_formattare)) {
  # Recupero le informazioni della tabella corrente.
  informazioni_tabella <- tabelle_da_formattare[[nome_tabella]]
  foglio <- informazioni_tabella$nome_foglio
  riga_iniziale <- informazioni_tabella$riga_iniziale
  colonna_iniziale <- informazioni_tabella$colonna_iniziale
  numero_righe_tabella <- informazioni_tabella$numero_righe
  numero_colonne_tabella <- informazioni_tabella$numero_colonne
  # Calcolo righe e colonne effettivamente occupate dalla tabella.
  righe_tabella <- riga_iniziale:(riga_iniziale + numero_righe_tabella - 1)
  colonne_tabella <- colonna_iniziale:(colonna_iniziale + numero_colonne_tabella - 1)
  riga_intestazione <- riga_iniziale
  riga_finale <- riga_iniziale + numero_righe_tabella - 1
  colonna_sinistra <- colonna_iniziale
  colonna_destra <- colonna_iniziale + numero_colonne_tabella - 1
  # - Bordo sottile su tutta la tabella
  addStyle(wb = file_risultati, sheet = foglio, style = stile_corpo_tabella, rows = righe_tabella, cols = colonne_tabella, gridExpand = TRUE, stack = TRUE)
  # - Stile intestazione con bordo inferiore spesso
  addStyle(wb = file_risultati, sheet = foglio, style = stile_intestazione, rows = riga_intestazione, cols = colonne_tabella, gridExpand = TRUE, stack = TRUE)
  # - Bordo esterno spesso della tabella
  #     Bordo superiore spesso.
  addStyle(wb = file_risultati, sheet = foglio, style = stile_bordo_superiore_spesso, rows = riga_iniziale, cols = colonne_tabella, gridExpand = TRUE, stack = TRUE)
  #     Bordo inferiore spesso.
  addStyle(wb = file_risultati, sheet = foglio, style = stile_bordo_inferiore_spesso, rows = riga_finale, cols = colonne_tabella, gridExpand = TRUE, stack = TRUE)
  #     Bordo sinistro spesso.
  addStyle(wb = file_risultati, sheet = foglio, style = stile_bordo_sinistro_spesso, rows = righe_tabella, cols = colonna_sinistra, gridExpand = TRUE, stack = TRUE)
  #     Bordo destro spesso.
  addStyle(wb = file_risultati, sheet = foglio, style = stile_bordo_destro_spesso, rows = righe_tabella, cols = colonna_destra, gridExpand = TRUE, stack = TRUE)
  # - Larghezza automatica delle colonne
  setColWidths(wb = file_risultati, sheet = foglio, cols = colonne_tabella, widths = "auto")
}
#----------------------------------------------------------------------------------------
# Salvataggio del file Excel
saveWorkbook(wb = file_risultati, file = "results/inference/Assignment_prima_parte_Inferenza_Amazon.xlsx", overwrite = TRUE)