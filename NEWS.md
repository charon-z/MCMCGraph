# BPFC 0.1.1

- Added `plot_mixing_trace()` for publication-ready mixing-proportion MCMC
  diagnostics with explicit axis labels, mathematical component notation,
  running means and posterior means.
- Added a resumable computational-performance benchmark covering BPFC-MCMC,
  BPFC-EM, Gaussian mixtures, K-means and agglomerative clustering, with
  end-to-end timing and peak-memory recording.
- Fixed benchmark summarization so methods without an MCMC iteration count are
  retained in the aggregated performance table.
- Declared testthat edition 3.

# BPFC 0.1.0

First public release accompanying the *Briefings in Bioinformatics*
manuscript.

## Features

* `run_mcmc_binary()` fits Bayesian finite-mixture functional clustering of
  paired (K = 2) longitudinal trajectories at a fixed number of clusters `J`.
* Paired cluster means use Legendre orthogonal polynomial bases (order 4).
* Within-individual residual correlation uses a block-diagonal SAD(1) covariance
  with a single shared time-correlation parameter `phi` and time-specific
  innovation variances.
* `fit_many_J()` / `eval_bic()` / `plot_bic()` provide BIC-based model selection
  over a grid of `J`.
* `plot_trace_density()` provides MCMC trace and density diagnostics.
* `as_binary_data()` accepts both wide (`n x 2d`) and long input layouts.
* Label switching is resolved with an ECR-style relabeling step; the fit returns
  MAP labels, posterior cluster-membership probabilities, and per-feature
  assignment uncertainty.

## Reproducibility

* `inst/reproduce/` contains the full simulation suite (stability, the
  joint/concatenation/independent ablation with two mirror scenarios, secondary
  ablations, and convergence diagnostics) plus the LINCS dose-response workflow.
* All simulation inputs are regenerated deterministically from stored scenario
  parameters and fixed seeds.
