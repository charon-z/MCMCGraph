# config.R ----------------------------------------------------------------
# Shared configuration for the MCMCGraph reproduce suite. All driver scripts
# (01-04) source this first. Paths are resolved relative to the project root,
# which is taken from the environment variable MCG_ROOT or, failing that, the
# current working directory.

MCG_ROOT     <- Sys.getenv("MCG_ROOT", unset = normalizePath(getwd()))
MCG_DATA_DIR <- file.path(MCG_ROOT, "benchmark", "data")     # sim_truth_K*.rds
MCG_OUT_DIR  <- file.path(MCG_ROOT, "results", "reproduce")  # all outputs land here
MCG_LIB_DIR  <- file.path(MCG_ROOT, "MCMCGraph", "inst", "reproduce")

dir.create(MCG_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (sub in c("sim1A", "sim1B", "sim1C", "sim1D", "tables", "figures"))
  dir.create(file.path(MCG_OUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)

# Experiment scale. Override any of these via environment variables so the same
# scripts run a quick R = 5 validation or the full R = 50 study.
MCG_R     <- as.integer(Sys.getenv("MCG_R",     unset = "50"))   # replicates
MCG_NITER <- as.integer(Sys.getenv("MCG_NITER", unset = "8000")) # MCMC iterations
MCG_EM_INIT <- as.integer(Sys.getenv("MCG_EM_INIT", unset = "10"))
MCG_SEED_BASE <- 20260530L

MCG_SCENARIOS <- c(3L, 5L, 8L)   # K_true for the three paper scenarios

# Per-replicate seed: distinct and reproducible across scenarios/replicates.
mcg_seed <- function(scenario, rep) MCG_SEED_BASE + scenario * 100000L + rep

# Which replicates this worker should run (for parallel splitting). Default: all.
mcg_reps <- function() {
  r <- Sys.getenv("MCG_REPS", unset = "")
  if (nzchar(r)) return(eval(parse(text = r)))
  seq_len(MCG_R)
}

mcg_scenarios <- function() {
  s <- Sys.getenv("MCG_SCENARIOS", unset = "")
  if (nzchar(s)) return(eval(parse(text = s)))
  MCG_SCENARIOS
}

message(sprintf("[config] ROOT=%s  R=%d  NITER=%d  EM_INIT=%d",
                MCG_ROOT, MCG_R, MCG_NITER, MCG_EM_INIT))
