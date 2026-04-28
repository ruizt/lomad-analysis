# Agent Instructions: Empirical Validation of Approximation Quality

## Objective

Assess the empirical accuracy of the specialization approximations for the local population correlation `rho` and the asymptotic variance term `V` under the signal-plus-noise model.

The goal is not to prove new theory. The goal is to determine when the approximations are accurate enough to support the methodology and to identify the regimes in which they degrade.

## Main questions

1.  How well does the common-trend approximation for `rho` match the empirical mean of the rolling correlation statistic?
2.  How well does the approximation for `V/s` match the empirical variance of the rolling correlation statistic?
3.  How do approximation errors change with:
    -   trend smoothness / local curvature,
    -   smoothing bandwidth `h`,
    -   local window size `s`,
    -   signal-to-noise ratio,
    -   degree of trend mismatch `delta`.

## Quantities to compare

For each simulation condition, estimate: - empirical mean of `R_t` across replicates, - empirical variance of `R_t` across replicates, - theoretical approximation for `rho`, - theoretical approximation for `V/s`.

Then compute: - absolute error in the mean approximation, - relative error in the mean approximation when sensible, - absolute error in the variance approximation, - relative error in the variance approximation when sensible.

## Simulation design

Use a focused factorial design rather than a large exhaustive grid.

### Core factors

-   series length `T`: choose 2 to 3 representative values,
-   smoothing bandwidth `h`: small / medium / large,
-   local window size `s`: small / medium / large,
-   signal-to-noise ratio: low / medium / high,
-   trend structure:
    -   shared trend,
    -   `delta`-similar trend with small mismatch,
    -   `delta`-similar trend with moderate mismatch,
-   trend smoothness:
    -   slowly varying,
    -   moderately varying,
    -   more locally curved.

### Noise model

Start with the same Gaussian noise model used in the paper's specialization section. Keep the dependence model simple and interpretable. AR(1) noise is sufficient for the main validation unless a richer model is already central to the paper.

### Replicates

Use enough replicates so that Monte Carlo error is clearly smaller than the approximation errors of interest.

## Recommended workflow

1.  Simulate trends and noise under one condition.
2.  Construct the observed rolling correlation statistic `R_t`.
3.  Compute the model-based approximation for `rho`.
4.  Compute the model-based approximation for `V`.
5.  Estimate the empirical mean and variance of `R_t` over replicates.
6.  Store errors and summary diagnostics.
7.  Repeat over all conditions.

## Minimal outputs

Produce the following:

### Tables

-   table of mean-approximation error by condition,
-   table of variance-approximation error by condition.

### Figures

-   empirical mean of `R_t` vs approximated `rho`,
-   empirical variance of `R_t` vs approximated `V/s`,
-   error vs `h`,
-   error vs `s`,
-   error vs signal-to-noise ratio,
-   error vs trend mismatch `delta`.

## Interpretation targets

Summarize results in terms of the following questions: - In what regimes is the approximation for `rho` reliable? - Is the approximation for `rho` consistently more accurate than the approximation for `V`? - When does the variance approximation begin to break down materially? - Are the most important determinants of approximation quality `h`, `s`, local curvature, or `delta`? - Is there a practical operating region in which the approximations are accurate enough for inference?

## Practical conclusion the study should support

The final deliverable should support one of the following statements: - the approximations are accurate in the regimes relevant to the paper, - the mean approximation is reliable but the variance approximation is only qualitative, - the approximations require empirical calibration outside a limited operating region.

## Constraints

-   Keep the study tightly scoped.
-   Prioritize interpretability over exhaustive coverage.
-   Do not introduce new methodology.
-   Treat this as a validation study for the existing specialization, not a second full simulation paper.

## Deliverable

Return: 1. a concise written summary of findings, 2. tables and plots, 3. a short recommendation on whether the manuscript should address approximation quality empirically in the main text, appendix, or not at all.
