# scenic

**Simulation and Causal Estimation for Negative Control Inference**

`scenic` is an R package for simulating data under a causal model with unmeasured confounding, and for benchmarking four causal inference estimators:

| Estimator | Function | Description |
|-----------|----------|-------------|
| UNADJ     | *(internal)* | Unadjusted OLS baseline |
| DIRECT    | `fit_direct()` | OLS with instrument + Negative Controls as covariates |
| COCA      | `fit_coca()` | Negative-control outcome correction (ratio estimator) |
| IV2SLS    | `fit_iv2sls()` | Two-stage least squares with genetic instrument |
| PGC       | `fit_pgc()` | Proxy G-component correction (3-step) |

---

## Causal model

```
U  ~ N(0,1)                          Unmeasured confounder
U2 ~ N(0,1)                          Extra noise source
G  ~ N(0,1)                          Genetic instrument (G -> Z only)

Z  = 0.6*G + conf_str*U + noise       Exposure (scaled)
M  = alpha_M * Z + noise              Mediator
Y  = beta_M*M + beta_Z*Z + gamma*U + noise   Outcome
W  = w_signal*U + (1-w_signal)*U2 + noise    Negative control

True total effect of Z on Y: beta_Z + alpha_M * beta_M
```

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `beta_Z`  | 0.10    | Direct effect Z → Y |
| `alpha_M` | 0.50    | Z → Mediator path |
| `beta_M`  | 0.30    | Mediator → Y path |
| `conf_str`| 0.80    | Confounding strength (U → Z and U → Y) |
| `w_signal`| 0.70    | Quality of negative control W as U proxy |

---

## Package structure

```
scenic/
├── R/
│   ├── generate_data.R     # Internal: data-generating process
│   ├── estimators.R        # Exported: fit_direct, fit_coca, fit_iv2sls, fit_pgc
│   ├── run_methods.R       # Internal: loop over features, summarise
│   ├── simulation.R        # Exported: run_simulation, sweep_param, run_null_sim
│   └── plots.R             # Exported: plot_* helpers
├── scenic_cluster_analysis.R   # Top-level analysis script for cluster
├── DESCRIPTION
├── NAMESPACE
└── README.md
```
