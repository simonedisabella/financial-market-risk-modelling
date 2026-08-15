# Financial Market Risk Modelling

Group coursework in **R** on statistical inference, dependence modelling and portfolio tail risk. The original work has three main blocks: Newey-West inference and Monte Carlo coverage, Variance-Gamma marginals with a Gaussian copula, and a CreditRisk+-style loss simulation with Student-t dependence.

The public repository reorganises the scripts and omits raw market-data files. It does not turn the coursework into a trading project or claim that the credit block was calibrated to a real loan book.

## 1. Newey-West inference and AR-GARCH simulation

`R/01_inference_newey_west_ar_garch.R`

The first block implements Newey-West long-run variance estimation with Bartlett weights and uses it to build confidence intervals for the mean, standard deviation and left-tail VaR. It then fits an AR(1)-GARCH(1,1) model to Amazon returns and runs Monte Carlo experiments under Student-t and standardized chi-square innovations.

The simulation focuses on **coverage and interval width**: how dependence, sample size, innovation shape and GARCH persistence affect inference, and how much is lost by using an i.i.d. standard error when the simulated process is dependent.

![VaR coverage](figures/inference/coverage_var.png)

## 2. Variance-Gamma marginals and Gaussian copula

`R/02_variance_gamma_gaussian_copula.R`

Daily log returns for ten assets are fitted with Variance-Gamma marginals. Probability-integral transforms are mapped to Gaussian scores, the dependence matrix is estimated there, and the joint distribution is simulated with a Gaussian copula. Simulated **log returns are converted back to simple returns before portfolio aggregation**, so the equal-weight portfolio is formed on the correct return scale.

The `vgFit()` step is run on returns multiplied by 100 for numerical stability. The location, scale and skew parameters are mapped back to the original return scale before `pvg()` and `qvg()` are used.

![Simulated portfolio versus normal benchmark](figures/copula/portfolio_vs_normal.png)

## 3. CreditRisk+-style loss simulation with Student-t dependence

`R/03_creditrisk_t_copula.R`

This block starts from **20 synthetic credit buckets**. Their expected default counts and variances are deterministic coursework assumptions, chosen because no real credit-portfolio dataset was available. Gamma-Poisson mixing gives negative-binomial marginal default-count distributions; a Student-t copula is then used to impose cross-bucket dependence.

The simulation compares aggregate-loss distributions over several dependence levels, including a comonotonic benchmark, and reports VaR and Tail Conditional Expectation (TCE). In this synthetic setup, stronger dependence produces materially heavier portfolio loss tails; this is a model-specific result, not an empirical credit calibration.

![VaR versus dependence](figures/credit-risk/var_vs_dependence.png)

## Output archive

The repository contains **all 19 PNG outputs** from the coursework archive. They are indexed in [`figures/README.md`](figures/README.md). The three images above are only the compact selection shown on the landing page.

## Repository layout

```text
R/
├── 01_inference_newey_west_ar_garch.R
├── 02_variance_gamma_gaussian_copula.R
└── 03_creditrisk_t_copula.R

data/README.md
figures/
├── inference/
├── copula/
└── credit-risk/
```

Raw coursework inputs are not redistributed. The simulation sections are reproducible from code because the random-number generation is seeded; the empirical market-data sections require the original inputs described in [`data/README.md`](data/README.md).

## Technical notes

- Newey-West uses Bartlett weights `1 - j / (m + 1)` for maximum included lag `m`.
- The VaR inference exercise keeps the sign of the left-tail return quantile, matching the convention used in the coursework.
- A Gaussian copula has zero asymptotic tail dependence; the repository does not present it as a general heavy-tail dependence model simply because the marginals are Variance-Gamma.
- Credit-loss quantiles use the empirical inverse CDF (`quantile(..., type = 1)`) to respect the discrete loss support.
- TCE is the mean simulated loss conditional on the loss being at least the estimated VaR threshold.
- The CreditRisk+ bucket inputs are assumed/synthetic. The exercise studies the effect of the dependence structure; it is not an empirical credit calibration.

## Authorship

This repository contains **group coursework** and does not imply sole authorship.

Academic modelling on historical/synthetic data, not a production risk engine or live investment record.
