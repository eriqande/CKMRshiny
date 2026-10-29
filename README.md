# CKMRshiny

A teaching app for exploring parent-offspring close-kin mark-recapture. The
first scenario uses independent random mating, equal expected reproductive
output, and independent simple random samples of females, males, and offspring.
The second tab, **Randomly unequal reproduction**, adds independently assigned
gamma reproductive sampling weights, while keeping survey sampling uniform.
The third tab, **Size-based reproduction**, connects age and individual growth
to length-based reproductive sampling weights, still with uniform adult surveys.

## Run

From this repository in R:

```r
install.packages(c("shiny", "ggplot2", "tibble")) # once, if needed
shiny::runApp()
```

Run a trial, then choose a saved trial and replicate. The population and sampling
figures show that replicate; pair-count and estimation figures summarize the trial.
The size-based pair-trait comparison also follows the selected replicate.
Editing the preparation controls only affects the next trial. Every replicate
generates its own population, reproduction, and samples. All parents, including
those with zero offspring or zero sampled offspring, appear in the relevant
histograms.

Population/Reproduction and Sampling each have an independent button to switch
the histogram y-axis between linear and log10 scales. This changes the view only.
Log10 bars start at 0.5 parents so counts of one remain visible; empty bins are
omitted. Parents with zero offspring still appear in the occupied bin at x = 0.

Under the section 01 Population / Reproduction figure in all tabs, a table reports realized **Effective N**
and **Effective N / N** separately for females and males. These are calculated
from **all offspring in the selected replicate's population, not just the sample**.
If adult i has k_i offspring and K is the total cohort size, the same-parent
probability for two distinct offspring is `h = sum(k_i * (k_i - 1)) / (K * (K - 1))`.
Effective N is `1 / h`; its ratio uses the full census adult count, including
zero-offspring adults. This finite-cohort effective number can exceed census N.
It is infinite when no same-parent offspring pairs exist, and undefined when
there are fewer than two offspring. The calculation uses retained population
frequency summaries, requiring neither regeneration nor additional individuals.

The moment estimate for each sex is `sampled parents * sampled kids / pairs`.
The Poisson pseudo-likelihood assumes that the pair count has mean
`sampled parents * sampled kids / abundance`. Its continuous maximum gives the
same estimate. The exact pair-count model in the equal-reproduction scenario is
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

Each tab has its own trial selection, parameters, replicate controls, register,
and RDS/CSV exports. Simulations are not shared across tabs. The app retains up
to 30 trials across all tabs in the current session. Limits distinguish
temporary simulation memory from retained sample storage:

- **One replicate:** at most 10 million population individuals (`Nf + Nm + Nk`).
- **One trial:** at most 20 million sampled individuals (`num_reps * (Sf + Sm + Sk)`).
- **One session:** at most 40 million sampled individuals across saved trials in all tabs.
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

## Randomly unequal reproduction

Section **00** previews Gamma(shape, scale = 1) densities in the female and male
colors, with separate shape inputs and a linear/log10 density y-axis button.
Previews update immediately before any trial exists. The sidebar checkbox
**Use random reproductive sampling weights** exists only on this tab and is
checked by default. Unchecking it runs equal reproduction in this tab. Saved
trials retain the checkbox state and both shapes; editing previews never rewrites
an existing trial. The supported shape range is 0.05 to 1000, defaulting to 1.

Reproductive sampling weights are used to choose parents when generating
offspring. They do not represent body mass, and they do not affect which adults
are sampled for the observed data. Only relative sampling weights within sex
matter. With scale fixed at one, the raw mean is the shape and the theoretical
CV is `1 / sqrt(shape)`. Each replicate draws a new independent set of adult
reproductive sampling weights and holds them fixed during reproduction. The
normalized reproductive sampling weights have a Dirichlet distribution; the
marginal offspring-count vector for each sex is Dirichlet-multinomial.

Density previews show raw gamma densities to the 99.5th percentile. For shapes
below one, the singularity at zero is omitted by starting at 0.001 of the upper
plot range; the density is not renormalized. Exact zero densities are omitted on
log10 axes. Sampled adults retain their reproductive sampling weights, and each
weighted replicate retains compact population summaries of those sampling
weights (means, CVs, ranges, and quantiles) alongside the offspring summaries.

Sections **01-04** use the same figure and estimator functions as tab one.
Uniform adult sampling preserves the expected pair count, but reproductive
inequality can inflate its variance. The existing Poisson Wald intervals do not
account for this inflation. A Dirichlet-multinomial estimator is not yet added.

## Size-based reproduction

The introduction, **Bigger adults have more offspring**, distinguishes reproductive
sampling weights from adult survey selection. The latter remains uniform.
Female and male settings are independent, with default length exponents 3 and 1.

Section **00** provides von Bertalanffy growth and length-output curves, plus
numerical ERRO distribution previews (12,000 fixed draws per sex, independent
of trial RNG). Adult age is a shifted exponential with mortality rate Z,
conditioned on the minimum/maximum age interval; ages are truncated, not capped.
Individual Linf and K are correlated lognormal traits parameterized by arithmetic
means and CVs; the correlation control applies on the log scale. Length is
`Linf * (1 - exp(-K * (age - t0)))`, with `t0 < minimum adult age`.
Reproductive propensity is `length^exponent * epsilon`, where epsilon is
independent mean-one lognormal variation with the specified conditional CV.
Dividing by the realized mean propensity within each sex gives mean-one ERRO.
Only these relative reproductive sampling weights affect parental assignments.
A zero exponent and zero reproductive CV give equal reproduction.

Section **01** shows a random display subset of at most 500 adults per sex
on growth and length-ERRO plots, followed by ERRO histograms for **every** adult,
realized family-size histograms, and realized effective parent numbers.
Section **02** repeats the two relationship plots for the **actual survey samples**
with matching axes. Above 10,000 adults combined, every adult contributes to
an 80-by-80 display grid; point area represents the number in each cell.

Section **03** adds normalized population, sampled-adult, and identified-pair
histograms, toggled between ERRO and length. **Each identified pair contributes
one observation**; a parent with five matched offspring contributes five times.
The dark dashed curve gives the reproductive model expectation for pair parents,
averaging over reproduction and uniform surveys. For mean-one ERRO r with
population distribution f(r), this is `r * f(r)`. The expected probability of a
bin is its summed population ERRO divided by total ERRO, computed before
unsampled adults are discarded. The length expectation uses the same reproductive
weighting. It is not a conditional prediction based on one realized survey or
pedigree. Observed histograms are normalized separately; no outline is drawn
for a sex with zero pairs. Full-population bin edges are shared by all sources.
Each sex has a right-hand **Identified pairs per bin** axis with whole-number
ticks, using `count = density * identified pairs * bin width`. This axis applies
to the solid identified-pair outline; all distributions retain their density
interpretation on the left. The count axis is omitted when a sex has no matches.
Panels stack on narrow screens and remain side by side on wider displays.

Saved trials retain all sex-specific model settings, complete sampled adult
traits, up to 500 display adults per sex, and compact full-population histogram
counts and ERRO sums. Full populations regenerate exactly from stage seeds.
The display subset uses its own RNG stream so it cannot change reproduction
or survey sampling. Growth traits require more transient memory per adult than
the first two scenarios; individual-count limits are not exact byte limits.
Existing point estimators and Poisson intervals remain unchanged.

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
pop_effective_numbers(trial$reps[[1]])

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
Rscript tests/check-random-reproduction.R
Rscript tests/check-size-reproduction.R
```

Checks cover validation, RNG isolation, unequal sex sizes, census sampling,
zero matches, repeated matches to a parent, theoretical match means, the
likelihood maximum, archive round trips, compact/full equivalence, stage-seed
regeneration, resampling without changing reproduction, larger populations,
sample/population budget boundaries, plotting edge cases, and Shiny trial
immutability. Shiny session checks are skipped when Shiny is unavailable.
The second check script also covers gamma draws, sampling-weight rescaling,
uniform survey sampling, overdispersed pair counts, weighted regeneration,
preview boundaries, tab isolation, module namespaces, and shared storage budgets.

The size-based checks cover truncated ages, correlated growth, reproductive CV,
equal-reproduction limits, exact regeneration/resampling, repeated pair parents,
normalized histogram expectations, zero-match plots, and immutable trial settings.
