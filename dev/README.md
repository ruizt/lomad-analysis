# dev/

This folder is for development work — simulation studies, exploratory
analysis, and prototyping. It is **not** part of the `lomad` package build.

## Workflow for contributors

At the top of any notebook or script, load the package functions with:

```r
devtools::load_all()
```

This makes all functions in `R/` available without having to install the
package. After making changes to `R/`, re-run `devtools::load_all()` to
pick them up.

## Structure

```
dev/
├── notebooks/    # Quarto notebooks for simulation studies and analysis
├── scripts/      # Standalone R scripts for prototyping and exploration
└── simulation/   # Power study infrastructure (see simulation-study.md)
```

See `simulation-study.md` for a full description of the simulation study design and the expected contents of `dev/simulation/`.

## Keeping files out of version control

Two options, both covered by `.gitignore`:

- **Prefix with `_`** — for individual files alongside tracked ones (e.g. `_my_scratch.R`)
- **`scratch/` folder** — drop anything in `dev/scratch/` and it won't be tracked

## Adding new functions

If a helper function you write in a script or notebook turns out to be
generally useful, move it to `R/` and add a roxygen2 header. Run
`devtools::document()` to regenerate the documentation.
