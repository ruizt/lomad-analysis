# Realignment Study Design

## Purpose

What does an equality-based multiscale trend comparison (MSinference;
Khismatullina and Vogt) do under local affine similarity, and can it be
rescued by removing the affine map first? One showcase series pair, two
scenarios, four framings. Companion demonstration for the supplement; the
argument it supports is that $\delta_t$ is the right estimand, not that the
comparator is deficient on its own terms.

## Scenarios

Both scenarios satisfy the local null everywhere (no true separation, `d = 0`);
every rejection reported anywhere in the study is a false alarm. The base
trend is shared exactly between scenarios; they differ only in the affine map
applied to the second series.

| Scenario | Map | Interpretation |
|----------|-----|----------------|
| `eq` | $a_t \equiv 0,\ b_t \equiv 1$ | trend equality; the comparator's home turf |
| `af` | drifting $a_t, b_t$ (cap 1.5% per window) | local affine similarity; the paper's null |

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
second trend with the first, which shares the base trend exactly.

## Methods and framings

| Cell | Scenario | Method | Framing |
|------|----------|--------|---------|
| 1 | `eq` | lomad | $s = 50$ |
| 2 | `eq` | MSinference | global standardization (harmless here; direct application) |
| 3 | `af` | lomad | $s = 50$ |
| 4 | `af` | lomad | $s = 80$ |
| 5 | `af` | lomad | $s = 100$ |
| 6 | `af` | MSinference | global standardization: remove one overall affine map |
| 7 | `af` | MSinference | local realignment (raw): rolling regression of $y_2$ on $y_1$, divide out |
| 8 | `af` | MSinference | local realignment (smoothed): map estimated from the $h$-smoothed series, applied to $y_2$ |

lomad is run across window sizes because $s = 50$ at $\bar\lambda = 0.5$,
$\phi = 0.5$ is the hardest cell of the power study; calibration there alone
invites the objection that the test would not have rejected regardless of the
data. At $s = 100$ the within-window affine drift is about twice the
generating cap, so those cells test lomad against a mildly violated version of
its own null. There is no positive control here: that lomad detects real
separation is established by the power study over 18 cells and 1500 replicates
per structure, which a single additional cell could only weaken.

Cells 5 and 6 differ only in what the rolling regression sees, and they fail
differently. Cell 5 is the naive plug-in: errors-in-variables attenuation
shrinks $\hat b_t$ toward zero by the window factor
$\lambda_{w}/(1+\lambda_{w})$, so the estimated map is systematically wrong,
not merely noisy (realized mean $\hat b \approx 0.61$ against a true
$b_t \in [0.97, 1.15]$; the factor varies by window because $\tau_w^2$ does).
Cell 6 is the version a careful practitioner would try; smoothing raises the
SNR entering the regression and reduces the attenuation (mean
$\hat b \approx 0.72$) but induces serial correlation that cuts the effective
sample size per window, so the map estimate is *more* variable, not less
(sd 0.38 vs 0.31; 96 sign flips vs 41). Instability in both cells concentrates
where the trend is locally flat — $\text{cor}(\log\tau_w^2, |\hat b_t - b_t|)
= -0.84$ — the same regime the paper identifies as uninformative for any
similarity inference. In both cells, dividing by $\hat b_t$ also rescales the
noise by $1/\hat b_t^2$, so the realigned series violates the comparator's
stationary long-run variance assumption by construction.

Realignment uses centered windows of width *s* (more favorable than the
trailing windows lomad uses) and is unconstrained apart from a finiteness
guard; $\hat b_t$ diagnostics (range, sd, sign flips) are recorded per cell.

An optional fifth MSinference cell (`eq` + raw realignment: realignment breaks
the test even when no map needs removing) can be added to `MS_CELLS` in
`run-realignment.R` if a referee asks.

## Endpoints

Per cell: the global test decision, and the count of flagged intervals or
rejected windows. Reported as a single table; the paired figure shows both
simulated series.

`reject` is the only endpoint comparable across methods. For MSinference it is
the global multiscale decision $\hat\Psi > q_{1-\alpha}$; for lomad it is
whether any window survives the BY step-up. The `flagged`/`total` counts are
reported for texture but their denominators are not comparable: lomad's are
test windows (1197 at $s = 50$), MSinference's are (location, bandwidth) grid
points (15250), and one interval flagged by MSinference is not one window
rejected by lomad.

### What $\hat\Psi$ is

The compiled `.rds` carries `stat` and `quant` for the MSinference cells; both
are dropped from the display table, since they are `NA` for lomad and inviting
cross-method comparison of them would be misleading. For the record,
MSinference reports

$$
\hat\Psi = \max_{(u,b)\in G}\left\{\left|\hat\psi(u,b)\right| - \lambda(b)\right\},
\qquad
\lambda(b) = \sqrt{2\log\{1/(2b)\}}
$$

the maximum over the grid $G$ of location–bandwidth pairs of the standardized
kernel-weighted difference statistic $\hat\psi$, each penalized by the additive
scale correction $\lambda(b)$ that puts all bandwidths on a common footing.
Here $G$ holds 15250 pairs over 61 bandwidths, giving $\lambda(b) \in
[1.18, 2.88]$. The critical value $q_{1-\alpha}$ is the $(1-\alpha)$ quantile
of the same maximum under the Gaussian multiplier null. Because the correction
is subtracted, $\hat\Psi$ is routinely negative under the null — the `eq` cell
returns $-0.975$ — which simply means no grid point exceeded its own
correction.

Note that `multiscale_test()` reports `stat` as a `max()` over a pairwise
matrix padded with structural zeros, which floors the value at 0 for `n_ts = 2`;
`run-realignment.R` extracts the true $(1,2)$ pairwise statistic instead. Test
decisions are unaffected, since the critical value is positive.

## Expected outputs

- `results/_raw/ms-{scenario}-{framing}[-seed{t}.{n}].rds` — one cached file per
  MSinference cell (slow: ~30-90 min each at `sim_runs = 1000`); gitignored,
  archived with the other raw results on Zenodo
- `results/realignment-results[-seed{t}.{n}].rds` — compiled study object, one
  per draw, both tracked
- `../_img/sfig-realignment.png` (showcase draw only),
  `../_tbl/stbl-realignment.csv`, `../_tbl/stbl-realignment-seed7307.2411.csv` —
  via `realignment-results.R`

## Workflow

```bash
# Development / fast pass: lomad cells only, MSinference cells left pending
REALIGNMENT_DRY=1 Rscript simulations/realignment/run-realignment.R

# Full run: executes any MSinference cell without a cache file (~2 h total)
Rscript simulations/realignment/run-realignment.R

# Smoke test of the MSinference harness (coarse quantile, ~6 min/cell)
SIM_RUNS=50 Rscript simulations/realignment/run-realignment.R

# Figure and table from the compiled object
Rscript simulations/realignment/realignment-results.R
```

Cached cells are never recomputed; delete the corresponding file in
`results/_raw/` to force a rerun. A cell run at a reduced `SIM_RUNS` is
stamped as such and must be deleted before the full run (the runner warns).

### Robustness to the draw

The study rests on one series pair, so it is also run on an unrelated draw to
confirm the showcase pair is not a freak one:

```bash
REALIGNMENT_TREND_SEED=7307 REALIGNMENT_NOISE_SEED=2411 \
  Rscript simulations/realignment/run-realignment.R
```

Both the MSinference caches and the compiled object are namespaced by seed, so
an alternate draw never reads or overwrites the default one. The alternate
draw is reported as a table only; a second copy of the figure would not say
anything the numbers do not.

The alternate draw lands at the 14th percentile of realized slope drift
against the showcase pair's 54th (measured over 200 draws at
`affine_cap = 0.015`), which makes it the more informative of the two: global
standardization is *adequate* there, and local realignment still fails, which
separates the failure of the repair from the size of the drift being repaired.

Note that `sim_trends()` consumes RNG differently depending on argument
*type*: `affine_s = 50` and `affine_s = 50L` yield different trends from the
same seed. The integer literals in `run-realignment.R` are therefore
load-bearing for reproducing the cached cells.
