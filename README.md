# iconic

**I**nference with **C**ausal **O**bservational **N**egative-control **I**nstruments and **C**ontrols

An R package for benchmarking causal inference methods in multi-omic observational studies with unmeasured confounding, motivated by the SCENIC framework for studying maternal PFAS exposure and placental transcription.

---

## Background

Estimating causal effects of prenatal exposures on fetal outcomes is complicated by unmeasured confounding, factors like socioeconomic stress or diet that affect both maternal PFAS levels and placental gene expression. Standard regression adjusts for observed covariates but cannot remove this residual confounding.

The ICONIC framework leverages three sources of identification:

- **G** — a polygenic risk score for PFAS metabolism as a genetic instrument (Mendelian randomization)
- **W** — a panel of negative-control outcomes (transcripts not on the PFAS causal pathway) that share the same unmeasured confounders as Y
- **C** — observed covariates (fetal sex, gestational age, maternal ancestry)

This package provides a toy simulation where the ground truth is known exactly, allowing for a rigorous benchmarking of all four estimators.

---

## Installation

```r
remotes::install_github("karrixxa/iconic")
```

Dependencies: `AER`, `parallel` (both on CRAN).

---

## The Four Estimators

| Name    | Function       | Approach                                              |
|---------|----------------|-------------------------------------------------------|
| UNADJ   | *(internal)*   | Unadjusted OLS — bias reference                       |
| DIRECT  | `fit_direct()` | OLS with G, W, and covariates as controls             |
| COCA    | `fit_coca()`   | Negative-control ratio correction (delta method SE)   |
| IV/2SLS | `fit_iv2sls()` | Two-stage least squares using G as instrument         |
| PGC     | `fit_pgc()`    | Proxy G-component correction (3-step residualisation) |

---

## Causal Model

```
G ~ N(0,1)                               Genetic instrument (PFAS PRS)
U1, U2 ~ N(0,1)                          Unmeasured confounders
Z = scale(0.6*G + delta*U1 + eps)        Maternal PFAS (scaled)
M = alpha_M*Z + eps                       Mediator
Y_f = beta_M*M + beta_Z*Z + gamma_f*U1  Placental transcript f
W_f = omega*U1 + (1-omega)*U2 + eps      Negative-control transcript f

True total effect: tau = beta_Z + alpha_M * beta_M
```

---

## Key Findings (toy simulation)

- **IV/2SLS** is unbiased and correctly controls Type I error at 5% across all confounding strengths. Recommended for the real data analysis.
- **PGC** is unbiased on average but ~2× higher variance than IV/2SLS.
- **DIRECT** has structural bias that does not shrink with sample size.
- **COCA** becomes unstable at low proxy quality; Type I error approaches 100% under confounding.
- **UNADJ** is always severely biased; provided as a reference floor.

---

## Package Structure

```
iconic/
├── R/
│   ├── estimators.R    # fit_direct, fit_coca, fit_iv2sls, fit_pgc
│   ├── generate_data.R # generate_toy_data (internal DGP)
│   ├── run_methods.R   # run_methods, summarise_results (internal wrappers)
│   ├── simulation.R    # run_simulation, sweep_param, run_null_sim, sweep_null_by_conf
│   └── plots.R         # plotting helpers
├── tests/testthat/
├── DESCRIPTION
└── NAMESPACE
```

---

## License

MIT
