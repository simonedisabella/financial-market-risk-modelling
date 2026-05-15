# Copula-based Portfolio Simulation

**Variance-Gamma marginals and Gaussian copula** — joint distribution
estimation and Monte Carlo simulation of an equally-weighted portfolio of
10 assets.

University assignment for *Statistica dei Mercati Finanziari*
(MSc, Univ. Milano-Bicocca).

\---

## What this project does

Given the historical daily prices of 10 heterogeneous assets, the project:

1. computes daily log-returns;
2. fits a **Variance-Gamma marginal distribution** to each asset by maximum
likelihood, with a numerical rescaling trick to ensure the optimiser
converges to a meaningful local maximum;
3. validates the marginal fits by comparing model-implied moments with
empirical ones;
4. estimates the **Gaussian copula** correlation matrix Σ from the
probability-integral-transformed data;
5. simulates 10 000 joint scenarios via the **inversion method**
(multivariate normal → Gaussian copula → Variance-Gamma marginals);
6. derives the simulated distribution of the equally-weighted portfolio
and reports its mean, standard deviation, and historical VaR at 95%
and 99% confidence levels.

## Methodological highlights

* **Sklar's theorem** is the backbone: the joint distribution is decomposed
into marginals + dependence structure, estimated separately and
recombined.
* **Numerical stability of the MLE.** Daily log-returns live on a scale
of \~0.01, where `vgFit()` from the `VarianceGamma` package converges
poorly. The code estimates on `log\_ret \* 100` and scales the parameters
back using the scale-closure property of the Variance-Gamma family
(`cX \~ VG(c·vgC, c·σ, c·θ, ν)`). Without this step, the implied
standard deviation was up to \~3× the empirical one for several assets.
* **Goodness-of-fit diagnostic.** A dedicated block prints, for every
asset, the ratio between the model-implied standard deviation and the
empirical one. After rescaling all ratios fall in `\[0.93, 0.99]`.
* **Sanity check on the portfolio.** The simulated portfolio is compared
with the historical equally-weighted portfolio: mean, sd, and 1%
quantile match in order of magnitude, confirming the model is coherent
with the data.
* **Known limitation: tail dependence.** The simulated 99%-VaR is slightly
less severe than the historical one. This is structural: the Gaussian
copula does not capture tail dependence (joint crashes). This is the
motivation for alternatives such as the t-copula or Archimedean
copulas (Clayton/Gumbel/Frank), which lie outside the scope of the
assignment.

## Repository structure

```
.
├── README.md
├── .gitignore
├── data/
│   └── prezzi\_chiusura.csv         daily close prices (Yahoo Finance, 2020-2025)
├── python/
│   ├── download\_prezzi.py          one-shot Yahoo Finance downloader
│   └── requirements.txt
├── R/
│   └── copule\_assignment.R         main analysis script
└── output/
    ├── istogramma\_portafoglio.png  simulated portfolio histogram
    ├── confronto\_normale.png       simulated density vs normal
    └── risultati\_copule.xlsx       4 sheets: VG params, fit diagnostic,
                                    copula matrix, portfolio summary
```

## Universe

10 heterogeneous assets to obtain a non-trivial dependence structure:

* **Tech**: AAPL, MSFT, GOOGL, NVDA
* **Financials**: JPM, GS
* **Energy**: XOM
* **Consumer staples**: KO
* **Gold ETF**: GLD
* **Long-duration Treasuries ETF**: TLT

The mix ensures both strong intra-sector dependence (e.g. GS-JPM ≈ 0.80)
and weak/negative cross-asset dependence (TLT vs equities ≈ -0.10 to -0.23).

## How to run

### 1\. (Optional) refresh the data

```bash
cd python
pip install -r requirements.txt
python download\_prezzi.py
```

This writes `data/prezzi\_chiusura.csv`. The CSV is already in the
repository, so this step is only needed to update the date range.

### 2\. Run the R analysis

Open RStudio in the project root and run `R/copule\_assignment.R`. The
script reads `data/prezzi\_chiusura.csv` and writes its outputs to the
`output/` folder.

Required R packages:

```r
install.packages(c("VarianceGamma", "mvtnorm", "writexl"))
```

## Selected results

|Statistic|Value|
|-|-:|
|Portfolio daily mean|0.00079|
|Portfolio daily sd|0.01124|
|VaR 95% (1-day)|-0.01774|
|VaR 99% (1-day)|-0.02846|

Marginal fit diagnostic (sd ratio model / empirical) — all assets in
`\[0.93, 0.99]`, confirming a coherent fit.

## References

* Embrechts, McNeil, Straumann, *Correlation and Dependence in Risk
Management: Properties and Pitfalls* (2002) — Sklar's theorem and
copula construction by inversion.
* Madan, Carr, Chang, *The Variance Gamma Process and Option Pricing*
(1998) — definition and properties of the Variance-Gamma distribution.
* Course lecture notes *Le Copule* (De Capitani, Univ. Milano-Bicocca).

## License

MIT.

\---

*Author: D'Isabella Simone — https://www.linkedin.com/in/simone-d-isabella-465595270/*

