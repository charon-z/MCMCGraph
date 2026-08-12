# cran-comments

## Test environments
- local macOS (Apple silicon), R 4.6.0
- GitHub Actions R-CMD-check (see `.github/workflows/R-CMD-check.yaml`)

## R CMD check results

`R CMD check --as-cran` returns **0 errors, 0 warnings**, with the following NOTEs:

1. **New submission.** This is the first release of the package.

2. **Assignment to the global environment** in `R/zzz_sad_dist.R`
   (`dSADmvnorm`, `rSADmvnorm`).
   This is intentional and required. BPFC registers a user-defined
   distribution (the structured-antedependence multivariate normal) with the
   `nimble` MCMC engine. During model compilation, `nimble` resolves the
   density and random-generation functions *by name* from the global
   environment, so the package exposes these two functions there at load time
   via a single, idempotent, guarded helper (`.mcmcgraph_expose_sad_to_global`).
   No user objects are overwritten (the helper only assigns when the name is
   absent unless `force = TRUE`), and the package functions themselves carry
   the canonical definitions. This mirrors the registration pattern used by
   other `nimble`-based packages.

## Notes for reviewers
- `nimble` is in `Depends` (not `Imports`) because the engine must be attached
  to the search path for model compilation to resolve its internal options
  (e.g. `getNimbleOption`); this matches the convention of the `nimble`
  ecosystem (e.g. `nimbleEcology`, `nimbleHMC`).
- Examples that fit a model are wrapped in `\donttest{}` because each fit
  compiles C++ through `nimble` and exceeds the 5-second example budget.
