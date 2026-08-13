# Power Study Design

## Objective

Characterise how reliably the test detects a *local* departure from affine
similarity, as a function of how large that departure actually is, across three
trend structures and the nuisance factors that govern estimation difficulty
(signal-to-noise, autocorrelation, window size).

## The null

Two series satisfy the null when they agree up to a *locally* affine map,

$$\nu_{2t} = a_t + b_t\,\nu_{1t},$$

with $a_t$ and $b_t$ free to drift, so long as they are near constant within any
one window. Series related this way are alike locally while their levels and
scales may differ completely across the series as a whole. Correlation is
invariant to shifting and rescaling either series, so this — not trend equality
— is what the test is calibrated against, and it is what the simulation has to
generate.

## Data-generating process

```
X_{kt} = ν_{kt} + Z_{kt},   Z_{kt} ~ AR(1, φ, σ²),   k = 1, 2
```

Trends are built in two stages. A coupling structure first separates a pair of
Fourier-basis trends, with $d$ scaling how far apart they are pushed. The affine
layer is then applied to the second, $\nu_2 \leftarrow a_t + b_t\nu_2$, with the
coefficients drawn as smoothed random walks and rescaled so that neither moves
by more than `affine_cap` over any window of length $s$.

That cap is what makes the pair *locally* affine similar: within a window the
map is effectively constant, while across the series the coefficients accumulate
freely, so windows far apart see genuinely different maps. At the settings below
$b_t$ swings about 17% from end to end.

## Trend structures

| Label | Generator                   | Key parameter                        |
|-------|-----------------------------|--------------------------------------|
| `fr`  | `sim_trends(method = "fr")` | rate *r* = 0.01, `bump = "gaussian"` |
| `rs`  | `sim_trends(method = "rs")` | bandwidth *bw* = 50                  |
| `rm`  | `sim_trends(method = "rm")` | bandwidth *bw* = 50                  |

`fr` concentrates its separation into brief evenly spaced events; `rs` and `rm`
spread it over episodes at irregular times, `rm` allowing the trends to cross.

**Pulse shape for `fr`.** Events use a gaussian pulse rather than a shape-2
gamma. The gamma leaves zero with non-zero slope, so the coupling weight has a
corner at each event onset, and difference-based noise estimation cannot cancel
a corner — Hall and Van Keilegom (2003, eqn 2.4) require a bounded derivative.
What survives differencing enters the residual autocovariance as a positive
additive bias, amplified near the unit root since the long-run noise variance
goes as (1+φ)/(1−φ). The gaussian pulse is matched on width and smooth at onset.
`fr` remains the weakest structure at high autocorrelation, which is the
contrast it is in the design to provide.

## Parameters

| Parameter              | Value                         |
|------------------------|-------------------------------|
| Correlation window *s* | 50, 100, 150                  |
| Series length *T*      | 25*s* — i.e. 1250, 2500, 3750 |
| AR(1) coefficient *φ*  | 0.3, 0.5, 0.7                 |
| Target SNR *λ*         | 0.5, 1.5                      |
| Smoothing window *h*   | 5                             |
| Fourier basis *nb*     | 151                           |
| Affine cap / bandwidth | 1.5% per window / 0.5*T*      |
| Replicates *S*         | 500                           |
| Significance level *α* | 0.05                          |

162 cells in total.

**The window is the design factor; *T* follows it.** Accumulated affine drift
depends on how many windows a series contains, so fixing *T*/*s* = 25 holds the
drift constant while *s* varies, and window size is not confounded with how much
affine variation the method faced. *nb* = 151 holds the shortest basis period at
*s*/3 for every *s*, fixing trend smoothness relative to the window too.

Long series are what let a per-window cap matter. The cap has to be small enough
that the map is constant *within* a window — a larger one would violate the null
outright — so length is the only way a small local drift accumulates into
genuinely different transformations at distant points. Power at matched
separation does not itself depend on *T*.

**φ stops at 0.7.** Beyond that the moving-average smoother cannot track the
trend, leaving trend signal in the residuals; the AR estimator absorbs it as
inflated $\hat\phi$ and the test turns conservative. At φ = 0.8, $\hat\phi$ pegs
at its 0.99 clamp in 8% of `rs` replicates and 17–67% of `fr` ones, and `fr`
loses detection entirely. At 0.7 there is no pegging for `rs` and `fr` retains
power, and 0.3 / 0.5 / 0.7 is evenly spaced.

## The *d* grid

| structure | *d* values    |
|-----------|---------------|
| `rs`      | 0, 0.60, 1.55 |
| `rm`      | 0, 0.36, 0.93 |
| `fr`      | 0, 0.27, 0.70 |

*d* is a generator knob, not an effect size, and never appears in a figure.
Separation is linear in it, but the structures distribute separation differently
in time, so a shared *d* gives them roughly 2× different *local* separation. The
values above are one base grid times a constant per structure (1.00 / 0.60 /
0.45), chosen so all three span a common range of realized δ_t with medians near
0, 0.20 and 0.45. See `_notes/delta-calibration.md`.

Three values suffice. A single *d* already spans most of the range — q10 0.036
to q90 0.747 at `rs`, *d* = 0.60 — so the grid shifts the distribution rather
than creating the coverage.

## Estimands

The estimand is **local**: the probability of rejecting window *t* given its
true separation, P(reject | δ_t). Global "did anything fire" power is not
reported — with 1200 to 3600 windows per series it saturates immediately, and
the sup-null it corresponds to is not the hypothesis under test.

Ground truth is

$$\delta_t = \sqrt{1 - r_t^2}, \qquad
  r_t = \operatorname{Corr}_{W_t}\!\left(\nu_1^{(h)}, \nu_2^{(h)}\right),$$

the RMS distance between the trends once the local affine map is removed by
least squares, measured on the noise-free trends smoothed exactly as the
observed series are. The generating coefficients are deliberately not used to
remove the map: fixing the slope at its window average charges the base
separation twice once *d* > 0, and δ_t then exceeds 1.

For each replicate, `run_rep()` returns two components.

### Summary (one row per replicate)

- **detected** — did any window get flagged (recorded, not reported)
- **n_win** — number of testable windows
- **delta_med** — median realized δ_t
- **lambda1**, **lambda2** — realized per-series SNR on Proposition 1's definition
- **b_range** — realized affine drift, so the layer's magnitude is recoverable

### Windows (the unit of analysis)

One row per window per replicate:

- **t** — time index
- **delta** — ground-truth δ_t
- **lambda** — min(λ₁, λ₂), the local SNR
- **p_raw** — raw p-value, so the study can be rethresholded at another α
- **rejected** — BY decision at α = 0.05

Every window is stored. The records are highly redundant — windows overlap by
*s* − 1 points — but subsampling them saves only disk and would discard the
adjacency information needed to say anything about contiguous flagged regions.
The sweep comes to about 4 GB.

## Expected outputs

### Per cell

- `{cell}.rds` — metadata + per-replicate summary data frame
- `{cell}-windows.rds` — the per-window records above

### Aggregated (after `collect-results.R`)

- `results/simulations-power-results.rds` — every replicate, every cell
- `results/simulations-power-curves.rds` — rejection rate by δ_t bin: the local power curves
- `results/simulations-power-auc.rds` — concordance between rejection and δ_t

Both aggregates pool over *d*, whose only job is to populate the δ_t range.
`collect-results.R` reduces each file as it reads it and keeps only per-bin
tallies: the sweep holds ~194 million window rows, which cannot be held in
memory at once.

Concordance is P(a randomly chosen rejected window is more separated than a
randomly chosen non-rejected one), computed as a Mann-Whitney statistic from
1000-bin tallies.

### Figures

Built by `simulations/simulation-results.R`, not by this study, and written to
`simulations/_img/`.

------------------------------------------------------------------------

## File layout

```
power/
├── design.md              ← you are here
├── simulation-template.R  ← local illustration mirroring tide/sim.R; writes nothing
├── collect-results.R      ← assembles fetched per-cell files into the compiled artifacts
├── results/
│   ├── simulations-power-results.rds       ← tracked
│   ├── simulations-power-curves.rds        ← tracked
│   ├── simulations-power-auc.rds           ← tracked
│   ├── _simulations-power.zip              ← archive of _raw/, Zenodo only
│   └── _raw/              ← per-cell .rds + -windows.rds fetched from Tide, Zenodo only
└── tide/                  ← Kubernetes scaffolding
    ├── sim.R              ← simulation script (mounted into container)
    ├── submit.sh          ← runs the full pipeline (PVC → jobs → wait → fetch)
    ├── submit_sweep.sh    ← submits one Job per (struct, d, s, phi, snr)
    ├── pvc.yaml           ← shared storage (create once)
    ├── job.yaml           ← Job spec template, filled by envsubst (single source of truth)
    ├── accessor.yaml      ← lightweight pod for file retrieval
    ├── fetch.sh           ← copies results from PVC to local machine
    └── test_one_job.sh    ← submits one small job to smoke-test sim.R
```

`_`-prefixed entries are gitignored and archived on Zenodo; the compiled `.rds`
files are small enough to track, so a fresh clone can rebuild every figure
without downloading anything.

`tide/` holds everything that talks to the cluster. `collect-results.R` sits
outside it because it only reads local files. `job.yaml` is the one copy of the
Job spec: `submit_sweep.sh` and `test_one_job.sh` fill its placeholders with
`envsubst`, so the spec cannot drift between the sweep and the smoke test.

------------------------------------------------------------------------

## Workflow

### Local development

Source `simulation-template.R`, which mirrors `tide/sim.R` at *S* = 20 and
writes nothing. To inspect a single replicate:

``` r
source("simulations/power/simulation-template.R")
run_rep(d = 0.60, struct = "rs", s_win = 100L, phi = 0.5, snr = 1.5, seed = 12345)
```

`tide/sim.R` is the source of truth. Change it first, then mirror the change
into `simulation-template.R`.

### On the cluster

``` bash
kubectl apply -n cal-poly-ruiz -f simulations/power/tide/pvc.yaml
bash simulations/power/tide/test_one_job.sh     # smoke test one cell
bash simulations/power/tide/submit_sweep.sh     # 162 jobs
bash simulations/power/tide/fetch.sh            # PVC -> results/_raw
Rscript simulations/power/collect-results.R
```

The image carries `lomad` and its dependencies; `sim.R` is mounted from a
ConfigMap that `submit_sweep.sh` refreshes on every run, so editing `sim.R`
does not require rebuilding the image.
