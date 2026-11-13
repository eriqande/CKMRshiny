# Run from the repository root. This needs microscoping, but not Shiny.
source("R/laugh-simulation.R")

# Small, stationary population: exact ages 0-3, then a 4+ group.
# Both sexes share survival and reproductive schedules. These geometric numbers
# give annual survival 0.5 and integer age/sex targets after package expansion.
design <- list(
  nata = c(512, 256, 128, 64, 64),
  wata = c(0, 1, 2, 4, 8),
  pmatata = c(0, 0.5, 1, 1, 1),
  samp_nya = rbind(
    c(10, 6, 6, 6, 0), # Year 1: sample numbers at ages 0, 1, 2, 3, 4+
    c(10, 6, 6, 6, 0), # Year 2
    c(10, 6, 6, 6, 0)  # Year 3
  ),
  MAX_JUVE_AGE = 2,
  MIN_CATCH_CURVE_AGE = 0,
  beqy_POP = FALSE,
  etriv = 0.01,
  extraplus = FALSE
)

# The same list can go to the analytic function and the simulator.
expected <- do.call(microscoping::ckmr_laugh_test, design)
simulated <- do.call(simulate_laugh_test,
  c(design, list(num_reps = 100, seed = 1401)))

# Or run both and assemble a comparison, including Monte Carlo standard errors.
checked <- check_laugh_test(design, num_reps = 500, seed = 1401)
checked$comparison

# Inspect what was actually simulated, including rounding and terminal survival.
checked$simulation$preparation$ages
checked$simulation$preparation$samples
checked$simulation$preparation$diagnostics
head(checked$simulation$first$pedigree)
head(checked$simulation$first$samples)

# z_difference measures discrepancy relative to SIMULATION uncertainty, not
# estimator precision. No automatic pass/fail claim is made: persistent deviations
# deserve investigation, including rounding, terminal ages, and census timing.
# Samples are nonlethal observations; counts apply the package's POP/GGP timing
# restrictions. This does not test additional mortality from lethal sampling.
