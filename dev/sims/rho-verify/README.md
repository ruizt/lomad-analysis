# rho-verify

Numerical verification of the general $\rho_t$ formula (Appendix, Eq. `rho-general`) against empirical local sample correlation $R_t$ from simulated data.

## What it checks

1. **General formula** $\rho_t \approx (\tau_{\bar\nu}^2 - \tau_d^2) / (\tau_{\bar\nu}^2 + \tau_d^2 + \sigma^2)$ tracks $R_t$ well across coupling methods (smooth, rate).
2. **Dropping the cross-covariance** $\text{Cov}_W(\bar\nu, d)$ introduces negligible error.
3. **The $w_t$-parametrized form** (Eq. `rho-w`) requires $w_t$ to vary slowly relative to $s$; when this fails, the general formula is much more accurate.

## Usage

```r
devtools::load_all(".")
source("dev/sims/rho-verify/verify-rho-formula.R")
```
