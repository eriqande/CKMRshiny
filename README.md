# CKMRshiny

For a static, step-by-step companion to Mark Bravington's Acme Archipelago
analysis, see [the tidy R Markdown examples](tidy_rmd_examples/README.md).
The notebook starts from his original data, keeps preparation visible, and
checks agreement with the original likelihood when his packages are available.

One launcher offers **Parent-offspring (POP)**, **Siblings (HSP)**, and
**Sampling design / Laugh test**. The top buttons switch directly between
them; **Switch app** returns to the chooser. The apps remain in the same session,
so controls, separate saved-trial registries, and design-editor values survive
switching. A browser restart starts a new session.

The parent-offspring app's
first scenario uses independent random mating, equal expected reproductive
output, and independent simple random samples of females, males, and offspring.
The second tab, **Randomly unequal reproduction**, adds independently assigned
gamma reproductive sampling weights, while keeping survey sampling uniform.
The third tab, **Size-based reproduction**, connects age and individual growth
to length-based reproductive sampling weights, still with uniform adult surveys.
The fourth tab, **Breaking Stuff!**, uses the same growth/reproduction model but
adds separate female and male length windows for adult sampling. Its estimator
still assumes uniform sampling and targets total abundance, deliberately leaving
the resulting bias uncorrected.

## Run

From this repository in R:

```r
install.packages(c("shiny", "ggplot2", "tibble")) # once, if needed
shiny::runApp()
```

## Sampling design and expected kin counts

The third app prepares inputs for and runs `microscoping::ckmr_laugh_test()`.
Calculations require an installed version with the `samp_nya` matrix interface
(tested with 1.1.2); the input editor and exports work without it. Restart Shiny
after updating the package. Explanations beside each input describe its meaning
and the package's treatment, including the equal-sex and steady-state assumptions.

Four aligned age curves describe population numbers (both sexes), biological
weight, proportion mature, and relative sampling allocation scores. Draw a stroke
across ages with a mouse or touch pointer, or open a curve's numeric table to
edit values or paste the entire vector. Numeric tables have a bounded horizontal
scroll window, visible scroll track, and sticky Age/Value labels. Focus the
window and use arrow keys, Home, or End as a keyboard alternative to scrolling.
Presets, multiplication controls, and display-scale controls help refine a shape. Maturity is bounded by
zero and one; other values must be finite and nonnegative. Strokes send a single
completed edit to Shiny, and intervening ages are filled as the pointer moves.

The final age category is a plus group: `12+` means every age from 12 onward.
Its population entry is the total number of those animals, and its sampling
entry allocates samples to the whole group. **Set age grid** starts fresh example
curves on a new grid, retaining the current total abundance. The button is
deliberate; merely typing a new age does not reset a drawing. **Undo last curve
edit** restores previous curves and grids, retaining up to 25 changes. It does
not undo annual sample counts or comparison settings. **Set total abundance**
rescales the population shape; editing ages individually may change that total.

**Sampling years, numbers, and age composition** groups annual sample counts
with a year selector, the sampling editor, and allocation previews. Append a
number to the comma-separated counts to add a year. Curves belong to positions
in that list (Year 1, Year 2, etc.); zero counts retain gap years. The selected
year's colored line is prominent, with other active years muted behind it.
Drawing, numeric edits, pasting, presets, and scaling apply only to that year.
**Copy this curve to all active years** is an undoable operation and leaves
hidden curves unchanged.

Removing final entries hides those years without deleting their curves from
session memory. Restoring those year positions recovers their previous values;
a genuinely new year inherits the currently selected curve. Changing sample
counts alone never changes the curves. Age-grid changes and Reset all curves
explicitly reinitialize active and cached curves; Undo restores the snapshot.

Sampling scores are normalized separately for each year to obtain `samp_pata`,
the age composition of the sample, **not age-specific probabilities of catching
population individuals**. The prepared matrix has rows in consecutive-year order
and columns in age order, ending at the plus group, with explicit dimension names.
Each row sums to one, even when its year's sample count is zero. Multiplying a
whole row by a positive constant leaves that year's proportions unchanged.
Only active years enter the matrix and exports. The selected year's normalized
proportions appear below the editor; a heatmap and table show fractional expected
allocations across all active years. These are allocations rather than simulated
or rounded sample counts.

The juvenile comparison cutoff is distinct from maturity. Advanced explanations
cover the catch-curve fitting start and whether adults sampled in an offspring's
birth year can count as possible parents. The final view exposes all numeric
schedules, normalized sampling proportions, and biological weight times maturity.
CSV downloads contain one row per active year and age, including raw scores,
normalized proportions, annual effort, and expected allocations. The R-input
download contains exact function arguments: `samp_nya` is an unnamed year-by-age
matrix of sample **numbers**, obtained by multiplying each year's proportions
by its annual effort. Its example package call is commented out.

**Calculate expected kin pairs** runs the installed package on a snapshot of the
current inputs. It saves expected juvenile-offspring POPs, adult-offspring POPs,
juvenile HSPs and GGPs, and juvenile/adult sample allocations. Select earlier
calculations to compare them; editing inputs never changes a saved result.
The app indicates when the current inputs differ from the displayed calculation.
Results can be exported as CSV or as an RDS containing the input snapshot,
exact function arguments, package version, and preparation details. Values stay
in the Shiny session when switching between the three apps; export them to keep
them beyond that session. These are rough expected counts, not precision estimates.

The results also show the package's smoothed population and expanded plus group,
without changing the curves in the editor. A preflight uses its `smooth_n` helper
and its default plus-group expansion to bound the large work arrays before
calculation (five million cells in one principal array). Designs exceeding that
interactive budget can still be exported for separate R calculations. This
preflight is checked against 1.1.2; kin counts always come from the installed
`ckmr_laugh_test()`, with `extraplus = FALSE`.

`R/laugh-inputs.R` contains documented preparation/validation functions;
`R/laugh-calculation.R` converts the inputs and calls the installed package;
`R/laugh-app.R` holds the interface and session history. The pointer editor and
its styles are in `www/laugh-editor.js` and `www/laugh-editor.css`. Illustrative
defaults use a declining age distribution; they are not assessment data.

## Standalone simulation check of the laugh test

`R/laugh-simulation.R` is a small verification tool, independent of Shiny. It
borrows only microscoping's private demographic helpers (`smooth_n` and
`get_psurv`); it independently generates pedigrees and counts relationships.
Its demographic preparation currently follows microscoping 1.1.2. It does not
copy the package's kin-probability formulas or modify the installed package.

Start with [examples/check-laugh-test.R](examples/check-laugh-test.R). A single
named list describes both calculations:

```r
source("R/laugh-simulation.R")
design <- list(
  nata = c(512, 256, 128, 64, 64),
  wata = c(0, 1, 1, 1, 1),
  pmatata = c(0, 1, 1, 1, 1),
  samp_nya = rbind(c(10, 6, 6, 6, 0),
                  c(10, 6, 6, 6, 0),
                  c(10, 6, 6, 6, 0)),
  MAX_JUVE_AGE = 2
)
checked <- check_laugh_test(design, num_reps = 500, seed = 1401)
checked$comparison
checked$simulation$preparation$ages
checked$simulation$preparation$diagnostics
head(checked$simulation$first$pedigree)
head(checked$simulation$first$samples)
```

You can also call `do.call(microscoping::ckmr_laugh_test, design)` and
`do.call(simulate_laugh_test, c(design, list(num_reps = 100, seed = 1401)))`
separately. The legacy `samp_pata`/`samp_ny` input is supported too. Returned
tables are ordinary data frames. The seed is local: simulations restore the
caller's random state, including when they fail.

The simulator uses equal sex schedules and maintains fixed whole-number
age/sex abundances by uniform random thinning, with constant recruitment.
Each offspring receives independently drawn parents, proportional to
`wata * pmatata`. These are expected outputs, without persistent individual
maturity states, mate fidelity, or additional reproductive heterogeneity.
Age-zero output must be zero. Founders are replaced through three expanded
lifespans of burn-in before sampling, so sampled grandparent relationships
can be recovered. Fractional age/sex sample allocations are stochastically
rounded independently, preserving their expectations. Samples within a
year/sex/age category are drawn without replacement; recaptures across years
are allowed.
The retained pedigree includes birth and death years, allowing reconstruction
of the simulated censuses; a missing death year means alive at the final census.

The samples are nonlethal observations. Kin counts explicitly apply the
package's eligibility rules: a POP parent must be sampled after the offspring's
birth year (or in that year with `beqy_POP = TRUE`); a GGP grandparent must be
sampled strictly after the intermediate parent's birth year. This checks those
restricted comparison counts, not population changes caused by lethal sampling.
Shared-parent matches count both parents of full siblings; grandparent matches
count each pedigree path. HSP/GGP summaries include only different birth
cohorts and two juvenile sample events.

This is an approximate finite-population benchmark. Inspect the returned age
table and diagnostics: fractional population abundances are rounded per sex,
and the simulator kills individuals after the final expanded age. The package's
HSP survival calculation can extrapolate beyond that age. Consequently, Monte
Carlo error is not the only possible source of differences. The comparison
reports simulation means, Monte Carlo standard errors, and standardized
differences, without an automatic pass/fail label. A zero Monte Carlo SE gives
`NA` for the standardized difference, rather than claiming agreement. These
standard errors describe simulation uncertainty, not CKMR estimator precision.

Oversized examples stop with an explanatory error before population allocation.
Defaults allow at most 20,000 standing individuals, 200,000 pedigree records per
replicate, 500,000 sample-event pairs per replicate, and twenty million records
or pair checks across the batch. The comparison also checks the analytic work
grid before invoking microscoping. This tool intentionally supports only
`extraplus = FALSE`; it is intended for small implementation checks.

## Sibling experiments

The sibling app follows the same sections 00-04, colors, trial/replicate controls,
histogram scale toggles, estimate distributions, and sorted interval displays.
It has four tabs:

1. **Equal reproduction:** sibling comparisons within one randomly mated cohort.
2. **Lucky litters:** gamma reproductive sampling weights represent unequal
   contributions to the surviving offspring cohort, compared within one episode.
3. **Cross-cohort siblings:** the same adults all survive to a second episode,
   but their gamma sampling weights are independently redrawn. Only offspring
   from different episodes are compared.
4. **Persistent size effects:** the same adults grow for one year before episode
   two, retaining their individual Linf and K. A mean-one lognormal multiplier
   varies annual growth rate (default CV 0.25). Positive increments cannot shrink
   fish or overshoot their own Linf. Fecundities are recalculated from new lengths;
   additional reproductive output CV remains zero. Uncheck growth to compare with
   fixed adult traits; set growth-rate CV to zero for deterministic von B growth.
   All adults survive, even if they pass the initial maximum-age cutoff.

Adults are not survey-sampled in these experiments. Each offspring has an
independently selected mother and father; kin is identified perfectly and the
simulation distinguishes maternal from paternal sharing. A full sibling pair
contributes one maternal and one paternal match. The app explicitly notes that
brood-level single-sire mating would require different treatment.

Within-cohort comparisons number `C = choose(Sk, 2)`; across cohorts they number
`C = Sk1 * Sk2`. Matches are counted from family sizes, without enumerating all
comparisons: `choose(k_i, 2)` within a cohort or `k_i1 * k_i2` across cohorts.
Both naive estimators use `Nhat = C / matches`. The Poisson Wald intervals are
working approximations; family clustering can affect their coverage. Full-sibling
counts are also reported separately as a diagnostic.

The sharing probabilities conditional on adult sampling weights are
`sum(q_i1^2)` within a cohort and `sum(q_i1 * q_i2)` across cohorts. Independent
seasonal luck restores the latter's mean to `1/N`; persistent size differences
remain. Population effective numbers are calculated from **realized population
offspring counts**, including cross-cohort `1 / sum(k_i1*k_i2/(Nk1*Nk2))`.
These realized reciprocals can exceed census abundance.

Cross-cohort plots compare the same adults' ERROs or realized family sizes in
both episodes, including lengths in the size scenario. The population trait
figures can show either episode. Size-based diagnostic histograms use episode-one
traits and show population adults, parents represented per sampled offspring,
and parents represented per sibling match. Their dashed expectations use the
sample-size-weighted mixture of episode ERROs and the product ERRO1 * ERRO2,
respectively (the latter reduces to ERRO squared without growth). Relative ERRO
can decrease even when a fish grows, because it is normalized within sex and episode.
Parent traits and zero-family
histograms are diagnostic simulation truth, not additional survey observations.

Complete offspring samples, family counts, summaries, settings, and seeds are
retained. Unsampled populations are discarded; `hsp_regenerate_rep()` reconstructs
them. Both apps share 30 trial slots and a 40-million-sample session budget. A
sibling replicate's population budget counts the adults once plus offspring from
all its episodes, up to 10 million individuals. Retained data and temporary memory
also include summaries and display subsets; individual-count limits are not byte limits.

Sibling functions live in `R/hsp-simulation.R`, `R/hsp-plots.R`, and
`R/hsp-app.R`, with the chooser in `R/ckmr-app.R`. They reuse the POP adult
generation, reproduction, RNG isolation, growth controls, estimate summaries,
and estimation/interval plots. For example:

```r
source("R/pop-simulation.R")
source("R/hsp-simulation.R")
settings <- hsp_settings(scenario = "cross_lucky", female_shape = 0.2)
trial <- hsp_run_trial(settings, num_reps = 100, seed = 1001)
hsp_trial_table(list(trial))
full_rep <- hsp_regenerate_rep(trial$reps[[1]], retain_population = TRUE)
```

## Parent-offspring experiments

Run a trial, then choose a saved trial and replicate. The population and sampling
figures show that replicate; pair-count and estimation figures summarize the trial.
The size-based pair-trait comparison also follows the selected replicate.
Editing the preparation controls only affects the next trial. Every replicate
generates its own population, reproduction, and samples. All parents, including
those with zero offspring or zero sampled offspring, appear in the relevant
histograms.

In **Breaking Stuff!**, first run a uniform-sampling baseline. Move the two-ended
length sliders beneath the section 02 ERRO plots to preview eligibility. Drag
either endpoint or the whole interval. Excluded regions and display adults turn
gray; actual sample overlays are shown only when the windows match the saved
trial. Pair counts and estimation results remain unchanged while previewing.
Click **Do Size-Selective Sampling** to apply the windows to every replicate and
save a new trial linked to its source. This regenerates the same populations and
reproduction from their stage seeds, preserves the exact sampled offspring, and
samples adults uniformly without replacement inside the windows. Changing one
sex's window does not change the other sex's sample. Reapplying identical windows
reproduces identical samples; a full-range window restores the uniform baseline.
No trial is saved if any replicate has fewer eligible adults than requested.
Source IDs and applied windows appear in the register and exports. Derived
trials consume the same shared trial/sample budgets as new population trials.

The expected pair-parent curve in section 03 uses ERRO sums restricted to eligible
adults for selected trials, while the population histogram still includes every
adult. This is a diagnostic expectation, not a correction to the abundance fit.

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
the regeneration helper. Length-window resampling is exposed in Breaking Stuff!;
other changes to sampling seeds and sizes remain available through the R API.
Sample data and summaries remain usable without regeneration.

## Checks

```sh
Rscript tests/check-pop.R
Rscript tests/check-random-reproduction.R
Rscript tests/check-size-reproduction.R
Rscript tests/check-selectivity.R
Rscript tests/check-hsp.R
Rscript tests/check-laugh-editor.R
Rscript tests/check-laugh-calculation.R
Rscript tests/check-laugh-simulation.R
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
The selectivity checks cover fixed offspring samples, independent sex-specific
sampling, full-window restoration, exact regeneration, eligible ERRO expectations,
failure without partial saves, induced estimation bias, and four-tab budgets.

The sibling checks compare family-count results with explicit pair enumeration,
check full-sibling double contributions and large count arithmetic, independent
seasonal luck, fixed size effects, census effective numbers, compact regeneration,
zero-match plots, scenario isolation, switching without losing trials, unique UI
IDs, and shared POP/HSP storage budgets. Growth checks cover deterministic von B
increments, nonshrinking and bounded stochastic increments, no mortality, updated
two-episode enrichment weights, and replay of legacy fixed-size archives.

The design-editor checks cover plus-group totals, row normalization, matrix
orientation and dimension names, year-specific expected allocations and gap years,
input bounds, stale-edit rejection, active-year editing, hidden-year restoration,
new-year inheritance, copy-to-all/undo, undo across age grids, clear/reset workflows,
and state preservation across the three-app launcher. Browser checks also cover
painting, the year dropdown, and scrolling all the way to an older plus group.
