# CKMRshiny

A teaching app for exploring parent-offspring close-kin mark-recapture. The
first scenario uses independent random mating, equal expected reproductive
output, and independent simple random samples of females, males, and offspring.

## Run

From this repository in R:

```r
install.packages(c("shiny", "ggplot2", "tibble")) # once, if needed
shiny::runApp()
```

Run a trial, then choose a saved trial and replicate. The population and sampling
figures show that replicate; pairs and estimation figures summarize the trial.
Editing the preparation controls only affects the next trial. Every replicate
generates its own population, reproduction, and samples. All parents, including
those with zero offspring or zero sampled offspring, appear in the relevant
histograms.

The moment estimate for each sex is `sampled parents * sampled kids / pairs`.
The Poisson pseudo-likelihood assumes that the pair count has mean
`sampled parents * sampled kids / abundance`. Its continuous maximum gives the
same estimate. The exact pair-count model in this equal-weight scenario is
binomial with `sampled kids` trials and probability `sampled parents / abundance`;
the Poisson is a working approximation. Zero matches give an unbounded estimate, explicitly counted in the
figures and table. Finite-estimate medians, means, and CVs exclude those replicates.
CV is the sample standard deviation divided by the finite mean, displayed as a
percentage; it is unavailable with fewer than two finite estimates. These are
conditional summaries and should be interpreted alongside the unbounded count.
For Poisson pseudo-MLE rows, the **Finite 95% CI inclusion** column reports the number of finite-replicate Wald
intervals containing the true abundance, divided by the number of finite
replicates, with the percentage alongside. This is coverage conditional on a
finite estimate; unbounded replicates are excluded from its denominator. With
no finite estimates, it displays `0 / 0 (unavailable)`. Moment-estimator rows
show a dash in this column.

The estimation area also shows approximate 95% Wald confidence intervals,
`Nhat +/- 1.96 * Nhat / sqrt(matches)`. Expected Fisher information under the
Poisson working model is `C / N^3`, where `C = sampled parents * sampled kids`.
The plug-in asymptotic variance is consequently `Nhat^2 / matches`. Intervals are
shown once per sex because both current estimators give identical results. They
are sorted from bottom to top by lower limit, then upper limit, then replicate
ID; the selected replicate has a thick line and a point. Red vertical lines mark
true abundance. Negative lower limits are preserved, and the normal approximation
can be poor with few matches. Zero-match replicates have no Wald interval and are
counted explicitly, including a note when the selected replicate has none.

The app retains up to 30 trials in the current session. Limits distinguish
temporary simulation memory from retained sample storage:

- **One replicate:** at most 10 million population individuals (`Nf + Nm + Nk`).
- **One trial:** at most 20 million sampled individuals (`num_reps * (Sf + Sm + Sk)`).
- **One session:** at most 40 million sampled individuals across saved trials.
- **Replication:** at most 500 replicates per trial.

Population size summed across replicates is not a storage limit. These counts
are practical memory safeguards rather than exact byte budgets; summaries,
pair tables, and estimates also take space. Nothing is silently discarded from
an existing trial. Trial archives can be downloaded as RDS and the
scenario-specific register as CSV. Closing the app
ends in-memory storage; download before closing. RDS files can be read and
reanalyzed in R; uploading/restoring them in the app is not implemented yet.

By default each replicate retains **complete samples plus population summaries**:
every sampled adult's attributes, every sampled offspring's pedigree, observed
kin pairs, estimates, and frequency tables of offspring counts for the full
population (including zero-offspring parents). Unsampled individual records are
discarded after each replicate is summarized. Figures use these compact results
directly without regenerating populations. RDS archives contain all retained
data and the seeds needed to regenerate full populations on demand.

## Reusable functions

`R/pop-simulation.R` separates settings validation, population generation,
reproduction, sampling, pair identification, estimation, and trial orchestration.
`R/pop-plots.R` provides summaries and plots. `R/pop-app.R` contains the UI and
session logic. Each function has roxygen documentation. This is not yet an R
package; no package installation or documentation generation is required.

The original exploratory `R/pop-fundamentals.R` is preserved. Shiny autoloading
of `R/` is disabled by `R/_disable_autoload.R`, and `app.R` explicitly sources
only the production files.

```r
source("R/pop-simulation.R")
source("R/pop-plots.R")
settings <- pop_settings(Nf = 1000, Nm = 1500, Nk = 3000,
                         Sf = 100, Sm = 100, Sk = 300)
trial <- pop_run_trial(settings, num_reps = 100, seed = 1001)
pop_plot_pairs(trial, rep_id = 1)
pop_trial_table(list(trial))

# Reproduce a compact replicate exactly from its saved stage seeds:
saved_rep <- trial$reps[[1]]
replicate <- pop_regenerate_rep(saved_rep)

# Regenerate the same population and reproduction, then take a new sample:
resampled <- pop_regenerate_rep(saved_rep, Sf = 200, Sm = 200, Sk = 300,
                               sampling_seed = 9001)
resampled$estimates

# Re-estimate using only the saved samples, without regenerating anything:
new_pairs <- pop_identify_pairs(saved_rep$samples)
pop_estimate(new_pairs, Sf = 100, Sm = 100, Sk = 300)

# Optionally retain the entire regenerated population and pedigree in R:
full_rep <- pop_regenerate_rep(saved_rep, retain_population = TRUE)
full_rep$population
full_rep$pedigree
```

Replicate `r` uses `first seed + r - 1` to derive three distinct stage seeds.
Population, reproduction, and sampling use offsets 0, 104729, and 209759,
wrapping within the positive R integer range. All three seeds are saved
explicitly. Changing a sampling seed never changes the population or its
reproduction. Simulation uses explicit Mersenne-Twister, Inversion, and Rejection
RNG settings and restores the caller's random state. Trials record schema
version 2, the R version, and RNG settings. Regeneration requires compatible
simulation code and RNG behavior; earlier trial schemas are not supported by
the regeneration helper. Regeneration/resampling is currently an R API, not an
app control. Sample data and summaries remain usable without regeneration.

## Checks

```sh
Rscript tests/check-pop.R
```

Checks cover validation, RNG isolation, unequal sex sizes, census sampling,
zero matches, repeated matches to a parent, theoretical match means, the
likelihood maximum, archive round trips, compact/full equivalence, stage-seed
regeneration, resampling without changing reproduction, larger populations,
sample/population budget boundaries, plotting edge cases, and Shiny trial
immutability. Shiny session checks are skipped when Shiny is unavailable.
