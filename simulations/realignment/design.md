# Realignment Study Design

## Purpose

Shows that local realignment — estimating the affine map on each window and
dividing it out — is not a workable preprocessing strategy: the estimated map
carries more error than it removes. Supplement demonstration.

## Scenarios

Both scenarios satisfy the local null everywhere (no true separation,
`d = 0`); every rejection reported anywhere in the study is a false alarm. The
base trend is shared exactly between scenarios; they differ only in the affine
map applied to the second series.

| Scenario | Map |
|----------|-----|
| `eq` | $a_t \equiv 0,\ b_t \equiv 1$ |
| `af` | drifting $a_t, b_t$ (cap 1.5% per window) |

Fixed parameters:

| Parameter | Value |
|-----------|-------|
| Series length *n* | 1250 |
| MA smoothing window *h* | 5 |
| Correlation window *s* | 50 |
| Trend | `sim_trends(n = 1250, d = 0, method = "rs", bw = 50, nb = 151, seed = 6001, affine_s = 50, affine_cap = 0.015)` |
| Noise | AR(1), $\phi = 0.5$, `lambda_target = 0.5`, `seed = 1001` |
| $\alpha$ | 0.05 |
| MSinference | `construct_grid(n)`, `estimate_lrv(q = 25, r_bar = 10, p = 1)`, `sim_runs = 1000` |

The `eq` scenario is derived from the `af` trend object by replacing the
second trend with the first.

## Methods and framings

| Cell | Scenario | Method | Framing |
|------|----------|--------|---------|
| 1 | `eq` | lomad | $s = 50$ |
| 2 | `eq` | MSinference | direct |
| 3 | `af` | lomad | $s = 50$ |
| 4 | `af` | lomad | $s = 80$ |
| 5 | `af` | lomad | $s = 100$ |
| 6 | `af` | MSinference | global standardization |
| 7 | `af` | MSinference | local realignment, map from the raw series |
| 8 | `af` | MSinference | local realignment, map from the $h$-smoothed series |

Realignment is a rolling regression of $y_2$ on $y_1$ over centered windows of
width *s*, divided out; unconstrained apart from a finiteness guard.

## Expected outputs

- `results/_raw/ms-{scenario}-{framing}[-seed{t}.{n}].rds` — one cached file per
  MSinference cell (~30-90 min each); gitignored, archived on Zenodo
- `results/realignment-results[-seed{t}.{n}].rds` — compiled study object, one
  per draw, both tracked
- `../_img/sfig-realignment.png` (showcase draw only),
  `../_tbl/stbl-realignment.csv`, `../_tbl/stbl-realignment-seed7307.2411.csv`

```bash
# showcase draw
Rscript simulations/realignment/run-realignment.R

# second draw
REALIGNMENT_TREND_SEED=7307 REALIGNMENT_NOISE_SEED=2411 \
  Rscript simulations/realignment/run-realignment.R

# figure and tables
Rscript simulations/realignment/realignment-results.R
```

Caches and compiled objects are namespaced by seed. Cached cells are never
recomputed; delete the file in `results/_raw/` to force a rerun.
