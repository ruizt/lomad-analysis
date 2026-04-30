# numerical-validation

Monte Carlo and numerical checks of the theoretical results in the paper. Each script targets one result.

## Scripts

| Script | Paper result | What it checks |
|--------|-------------|----------------|
| `validate-clt.R` | Theorems 1 & 2 | Sampling distribution of $R_t$ matches the predicted CLT normal approximation. Checks standardized statistic against N(0,1) via QQ plot, Shapiro-Wilk, and interval coverage at multiple window sizes. |
| `validate-prop1-rho.R` | Proposition 1 ($\rho$) | Theoretical $\rho_t$ matches the mean empirical $R_t$ across replications. Checks convergence as $s$ increases. |
| `validate-prop1-V.R` | Proposition 1 ($V$) | Theoretical $V_t$ matches $s \cdot \text{Var}(R_t)$ across replications. |
| `validate-perturbation.R` | Perturbation bound | Under near-common trends ($|s_1 - s_2| < \varepsilon$), $\rho$ is a perturbation of the common-trend expression with remainder $O(\varepsilon\tau)$. Varies $\varepsilon$ from 0 to large. |
| `validate-gradient.R` | Gradient remark (Thm 2) | Analytic $\nabla g(\theta)$ matches numerical finite differences. Also verifies the centering argument: full 5×5 $V$ equals reduced 3×3 $V$. |
| `validate-sigma-matrix.R` | $\Sigma_{(3,4,5)}$ (Prop 1 proof) | MC verification that the Wick's-theorem entries match empirical covariances in centered coordinates. Tests white and ARMA noise. |
| `validate-rho-formula.R` | Appendix: $\rho_t$ under simulation model | General formula $\rho_t \approx (\tau_{\bar\nu}^2 - \tau_d^2)/(\tau_{\bar\nu}^2 + \tau_d^2 + \sigma^2)$ tracks empirical $R_t$ under smooth and rate coupling. |
| `validate-end-to-end.R` | Full pipeline | Runs `lomad_fit()` + `lomad_test()` on simulated data with known parameters. Compares estimated $\hat\rho$, $\hat V$, and CLT coverage against oracle values. Tests whether estimation error degrades the asymptotic approximation. |

## Usage

Run from the package root:

```r
devtools::load_all(".")
source("dev/numerical-validation/validate-clt.R")
# etc.
```

Scripts that generate plots save them as `.png` files in this directory.
