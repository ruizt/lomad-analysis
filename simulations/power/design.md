# Power Study Design

## Objective

Characterize how reliably the test detects a local departure from affine
similarity, as a function of how large that departure actually is, across three
trend structures and the nuisance factors that govern estimation difficulty
(signal-to-noise, autocorrelation, window size).

## Local affine similarity

Two series satisfy local affine similarity if $\nu_{2u} = a_t + b_t\nu_{1u}$ on 
a window $W_t$. The coefficients $a_t$ and $b_t$ free to drift with $t$, so long 
as they are near constant within any one window. 

Departures from local affine similarity are measured by mean square separation 
after affine realignment on a window $W_t$ of width $s$:

$$
\delta_t^2 = \frac{1}{\tau_{1t}^2} \min_{a_t, b_t > 0} 
\left\{ 
  \frac{1}{s}\sum_{u\in W_t}\left[ \nu_{1u} -  (a_t + b_t\nu_{2u})\right]^2 
\right\}
$$

We are interested in power relative to the effect size $\delta_t$.

## Data-generating process

$$
X_{kt} = \nu_{kt} + Z_{kt},   
\qquad
Z_{kt} \sim AR(1, \phi, \sigma^2),
\qquad k = 1, 2
$$

Trends are built in three stages:

1. Generate a pair of Fourier-basis trends with total $L^2$ separation of $d$.
2. Mix trends according to a mixing weight.
3. Apply a time-varying affine transformation $\nu_2 \leftarrow a_t + b_t\nu_2$, 
with $a_t, b_t$ drawn as smoothed random walks and rescaled so that neither moves
by more than `affine_cap` over any window of length $s$.

The cap is calibrated to preserve local affine similarity under $H_0$, but allows
for global drift of the affine map. At the settings below $b_t$ swings about 17% 
from end to end.

## Trend structures

Separation is generated globally but detected locally, and the two are not
interchangeable. A window of width $s$ sees at most $s$ terms of the global sum
defining $d$, so

$$\delta \leq d/\sqrt{s},$$

with equality approached only when the separation falls inside a single window.
The same total dissimilarity is therefore easy or hard to find depending on
whether it is concentrated into brief episodes or spread thinly across the
series. **How separation is distributed in time, not just how much of it there
is, is what this study varies.**

The distribution is controlled by the mixing weight of stage 2, which mixes the
base trends about their midpoint: $w_t \approx 1$ leaves them coupled, and
falling $w_t$ drives them apart. Three structures specify $w_t$ differently:

| Label | Generator                   | Key parameter                          |
|-------|-----------------------------|----------------------------------------|
| `fr`  | `sim_trends(method = "fr")` | rate *r* = 0.4/*s*, `bump = "gaussian"` |
| `rs`  | `sim_trends(method = "rs")` | bandwidth *bw* = 50                    |
| `rm`  | `sim_trends(method = "rm")` | bandwidth *bw* = 50                    |

`fr` concentrates its separation into brief evenly spaced events; `rs` and `rm`
spread it over episodes at irregular times, with `rm` allowing the trends to cross.

`fr`'s rate is set from the window rather than fixed: the pulse width it implies
scales as $1/(r^2T)$, so a constant rate shrinks the events as the design grows
until they fall below the smoothing bandwidth and are erased before the test
sees them. At $r = 0.4/s$ there are 10 events of width $s/4$ at every $s$.

## Simulation design

We simulate data according to factorial combinations of the following:

| Parameter              | Value                         |
|------------------------|-------------------------------|
| Correlation window *s* | 50, 100, 150                  |
| Series length *T*      | 25*s* — i.e. 1250, 2500, 3750 |
| AR(1) coefficient *φ*  | 0.3, 0.5, 0.7                 |
| Target SNR *λ*         | 0.5, 1.5                      |

This results in 162 cells in total. Throughout, we fix the following:

| Parameter               | Value                         |
|-------------------------|-------------------------------|
| Smoothing bandwidth *h* | 5                             |
| Fourier basis *nb*      | 151                           |
| Affine cap / bandwidth  | 1.5% per window / 0.5*T*      |
| Replicates *S*          | 500                           |
| Significance level *α*  | 0.05                          |

To generate data with varying effect sizes, we draw the initial Fourier trends
(before mixing and affine transformation) at the following levels of separation:

| structure | *d* values    |
|-----------|---------------|
| `rs`      | 0, 0.60, 1.55 |
| `rm`      | 0, 0.36, 0.93 |
| `fr`      | 0, 0.5, 1.7   |

*d* is a generator knob, not an effect size; the structures distribute separation
differently in time, so a shared *d* gives them roughly 2× different *local*
separation. Each structure instead takes the values that put its realized δ_t on
a common footing, with medians near 0, 0.20 and 0.50. `rs` and `rm` differ by a
constant factor (1.00 / 0.60); `fr` does not, because its δ_t saturates near
0.43 — separation confined to isolated events leaves most windows coupled, so
the median cannot climb further.

Three levels suffice: a single *d* already spans q25 0.10 to q75 0.42, so the
grid shifts the distribution rather than creating the coverage.

## Estimands

The estimand is **local**: the probability of rejecting window *t* given its
true separation, P(reject | δ_t).

Ground truth is given by $\delta_t = \sqrt{1 - (r_t\vee 0)^2}$ where 
$r_t = \operatorname{Corr}_{W_t}\!\left(\nu_1, \nu_2\right)$ is the empirical
correlation between true trends on $W_t$. See the methods paper for details. 

## Implementation

Replicates are executed by the wrapper `run_rep()`, which returns two components.

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

Every window is stored. The sweep comes to about 4 GB.

Expected outputs are stored per cell:

- `{cell}.rds` — metadata + per-replicate summary data frame
- `{cell}-windows.rds` — the per-window records above

The script `collect-results.R` then aggregates these and stores the files:

- `results/simulations-power-results.rds` — every replicate, every cell
- `results/simulations-power-curves.rds` — rejection rate by δ_t bin: the local power curves
- `results/simulations-power-auc.rds` — concordance between rejection and δ_t

Aggregation pools over *d* and bins by $\delta_t$.

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

### Local testing

Source `simulation-template.R`, which mirrors `tide/sim.R` at *S* = 20 and
writes nothing. To inspect a single replicate:

``` r
source("simulations/power/simulation-template.R")
run_rep(d = 0.60, struct = "rs", s_win = 100L, phi = 0.5, snr = 1.5, seed = 12345)
```

`tide/sim.R` is the source of truth. Change it first, then mirror the change
into `simulation-template.R`.

### Deployment

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
