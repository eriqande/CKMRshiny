#' Validate the settings for a simple parent-offspring experiment
#'
#' All sizes must be positive whole numbers. Parents and offspring are sampled
#' without replacement from three disjoint groups. Each replicate generates a
#' new population and pedigree; the population sizes are fixed across replicates.
#' @param Nf,Nm Numbers of adult females and males.
#' @param Nk Number of offspring in the population.
#' @param Sf,Sm,Sk Numbers of females, males, and offspring sampled.
#' @param use_weights Assign independent gamma reproductive sampling weights.
#' @param female_shape,male_shape Gamma shapes; scale is fixed at one. Supported
#' shapes range from 0.05 to 1000 to avoid numerical extremes in draws/previews.
#' @param size_model Optional list with Female and Male settings from
#' pop_size_parameters(). Overrides gamma reproduction when supplied.
#' @return A named list of validated sizes and reproduction settings.
#' @export
pop_settings <- function(Nf = 1000, Nm = 1000, Nk = 2000,
                         Sf = 100, Sm = 100, Sk = 200,
                         use_weights = FALSE, female_shape = 1, male_shape = 1,
                         size_model = NULL) {
  values <- list(Nf = Nf, Nm = Nm, Nk = Nk, Sf = Sf, Sm = Sm, Sk = Sk)
  for (name in names(values)) {
    .pop_integer(values[[name]], name)
  }
  if (Sf > Nf || Sm > Nm || Sk > Nk) {
    stop("Each sample size must be no larger than its population size.", call. = FALSE)
  }
  if (!is.logical(use_weights) || length(use_weights) != 1L || is.na(use_weights)) {
    stop("use_weights must be TRUE or FALSE.", call. = FALSE)
  }
  .pop_shape(female_shape, "Female gamma shape")
  .pop_shape(male_shape, "Male gamma shape")
  result <- c(lapply(values, as.integer), list(use_weights = use_weights,
                                    female_shape = female_shape, male_shape = male_shape))
  if (!is.null(size_model)) {
    if (!is.list(size_model) || !setequal(names(size_model), c("Female", "Male")))
      stop("size_model needs Female and Male parameter lists.", call. = FALSE)
    result$size_model <- lapply(size_model[c("Female", "Male")], function(p)
      do.call(pop_size_parameters, p))
  }
  result
}

#' Validate a gamma shape in the supported educational range
#' @param x Numeric scalar shape.
#' @param name Label used in error messages.
#' @return Invisibly returns x or raises an informative error.
#' @keywords internal
.pop_shape <- function(x, name = "Gamma shape") {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0.05 || x > 1000) {
    stop(name, " must be between 0.05 and 1000.", call. = FALSE)
  }
  invisible(x)
}

#' Check a positive integer scalar
#' @param x Value to check.
#' @param name Label included in an error message.
#' @return Invisibly returns the original value, or raises an informative error.
#' @keywords internal
.pop_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be a positive whole number no larger than ",
         .Machine$integer.max, ".", call. = FALSE)
  }
  invisible(x)
}

#' Evaluate an expression with a reproducible, locally scoped random seed
#'
#' Explicit RNG kinds make saved seeds independent of the caller's RNG settings.
#' The caller's RNG kinds and random state (including an absent state) are restored
#' even if the expression fails. This helper does not affect parallel RNG streams.
#' @param seed Positive integer seed.
#' @param code Expression evaluated in the caller's environment.
#' @return The value of the expression.
#' @keywords internal
.pop_with_seed <- function(seed, code) {
  .pop_integer(seed, "seed")
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(code)
}

#' Generate adults with equal or gamma reproductive sampling weights
#'
#' With equal reproduction every reproductive sampling weight is one. Otherwise,
#' each adult independently receives Gamma(shape, scale = 1), fixed for that
#' replicate's reproduction. Only relative sampling weights within sex affect
#' parent assignment; survey sampling remains uniform. Parent IDs are unique within sex;
#' the pair (sex, parent_id) identifies an adult. Offspring form a separate group.
#' @param Nf,Nm Numbers of females and males.
#' @param use_weights Whether to assign random reproductive sampling weights.
#' @param female_shape,male_shape Gamma shapes, each between 0.05 and 1000.
#' @return A tibble with parent_id, sex, and reproductive_weight columns. The
#' latter stores raw reproductive sampling weights, not biological body mass.
#' @export
pop_generate_population <- function(Nf, Nm, use_weights = FALSE,
                                    female_shape = 1, male_shape = 1) {
  .pop_integer(Nf, "Nf")
  .pop_integer(Nm, "Nm")
  .pop_shape(female_shape, "Female gamma shape")
  .pop_shape(male_shape, "Male gamma shape")
  if (!is.logical(use_weights) || length(use_weights) != 1L || is.na(use_weights)) {
    stop("use_weights must be TRUE or FALSE.", call. = FALSE)
  }
  sampling_weights <- if (use_weights) {
    c(stats::rgamma(Nf, shape = female_shape, scale = 1),
      stats::rgamma(Nm, shape = male_shape, scale = 1))
  } else rep(1, Nf + Nm)
  tibble::tibble(parent_id = c(seq_len(Nf), seq_len(Nm)),
                 sex = rep(c("Female", "Male"), c(Nf, Nm)),
                 reproductive_weight = sampling_weights)
}

#' Generate one cohort of offspring by independent parental assignments
#'
#' Each offspring receives one mother and one father, independently, with
#' probabilities proportional to the adults' reproductive sampling weights. There are no
#' pair bonds, mortality, overlapping cohorts, or errors in kin identification.
#' @param population Adult tibble from pop_generate_population().
#' @param Nk Number of offspring to produce.
#' @return A tibble with kid_id, mother_id, and father_id columns.
#' @export
pop_reproduce <- function(population, Nk) {
  .pop_integer(Nk, "Nk")
  draw <- function(sex) {
    # Sample row indices: sample(x) has surprising behavior for a single ID.
    parents <- population[population$sex == sex, ]
    if (!nrow(parents) || any(!is.finite(parents$reproductive_weight)) ||
        any(parents$reproductive_weight < 0) || sum(parents$reproductive_weight) <= 0) {
      stop("Each sex needs finite, nonnegative reproductive sampling weights and a positive total.",
           call. = FALSE)
    }
    parents$parent_id[sample.int(nrow(parents), Nk, replace = TRUE,
                                  prob = parents$reproductive_weight)]
  }
  tibble::tibble(kid_id = seq_len(Nk), mother_id = draw("Female"),
                 father_id = draw("Male"))
}

#' Draw simple random samples of parents and offspring
#'
#' All three samples are independent, without replacement. Reusing a saved
#' population and pedigree with a new locally scoped seed allows resampling
#' without regenerating reproduction.
#' @param population Adult population tibble.
#' @param pedigree Offspring tibble from pop_reproduce().
#' @param Sf,Sm,Sk Sample sizes for females, males, and offspring.
#' @return A list containing females and males as complete adult tibbles, and
#' kids as a complete sampled pedigree tibble. Unsampled adults are not retained.
#' @export
pop_sample <- function(population, pedigree, Sf, Sm, Sk) {
  females <- population$parent_id[population$sex == "Female"]
  males <- population$parent_id[population$sex == "Male"]
  pop_settings(length(females), length(males), nrow(pedigree), Sf, Sm, Sk)
  female_rows <- which(population$sex == "Female")
  male_rows <- which(population$sex == "Male")
  list(females = population[female_rows[sample.int(length(females), Sf)], ],
       males = population[male_rows[sample.int(length(males), Sm)], ],
       kids = pedigree[sample.int(nrow(pedigree), Sk), ])
}

#' Identify all parent-offspring pairs in the samples
#'
#' Kin is observed perfectly. A parent matched to three offspring contributes
#' three rows. An offspring may contribute both a maternal and a paternal pair.
#' @param samples Complete samples from pop_sample(). Only sampled data are used.
#' @return A tibble with sex, parent_id, and kid_id, including zero-row results.
#' @export
pop_identify_pairs <- function(samples) {
  kids <- samples$kids
  mothers <- kids$mother_id %in% samples$females$parent_id
  fathers <- kids$father_id %in% samples$males$parent_id
  tibble::tibble(sex = c(rep("Female", sum(mothers)), rep("Male", sum(fathers))),
                 parent_id = c(kids$mother_id[mothers], kids$father_id[fathers]),
                 kid_id = c(kids$kid_id[mothers], kids$kid_id[fathers]))
}

#' Estimate female and male abundance from parent-offspring matches
#'
#' The moment estimator is S_parent * S_kid / matches. The Poisson
#' pseudo-likelihood treats the match count as Poisson with mean
#' S_parent * S_kid / N, with N >= S_parent and continuous N. Its maximum is at
#' the same estimate; it is calculated analytically rather than optimized
#' numerically. In this equal-weight scenario the exact count distribution is
#' binomial(S_kid, S_parent / N); the Poisson is a working approximation. No
#' uncertainty intervals are computed here. Zero matches give an unbounded
#' estimate for both methods.
#' @param pairs Observed pair tibble from pop_identify_pairs().
#' @param Sf,Sm,Sk Sample sizes.
#' @return A tibble with sex, method, matches, estimate, and finite columns.
#' @export
pop_estimate <- function(pairs, Sf, Sm, Sk) {
  for (name in c("Sf", "Sm", "Sk")) .pop_integer(get(name), name)
  matches <- vapply(c("Female", "Male"), function(s) sum(pairs$sex == s), integer(1))
  if (any(matches > Sk)) stop("A sex cannot have more pairs than sampled offspring.")
  estimates <- c(Sf, Sm) * as.double(Sk) / matches
  tibble::tibble(sex = rep(c("Female", "Male"), 2),
                 method = rep(c("Moments", "Poisson pseudo-MLE"), each = 2),
                 matches = rep(unname(matches), 2), estimate = rep(unname(estimates), 2),
                 finite = is.finite(estimate))
}

#' Summarize realized reproductive output, including zero-offspring parents
#'
#' This frequency table is sufficient for the population reproduction figures.
#' It contains one row per observed offspring count and sex, rather than one row
#' per adult. It can also summarize sampled parents against sampled offspring;
#' offspring assigned to parents outside the supplied adults are ignored.
#' @param population Adult tibble, with IDs unique within sex.
#' @param pedigree Offspring pedigree tibble.
#' @return A tibble with sex, offspring, and n_parents (frequency) columns.
#' @export
pop_summarize_reproduction <- function(population, pedigree) {
  tibble::as_tibble(do.call(rbind, lapply(c("Female", "Male"), function(sex) {
    ids <- population$parent_id[population$sex == sex]
    parent_ids <- if (sex == "Female") pedigree$mother_id else pedigree$father_id
    counts <- tabulate(match(parent_ids, ids), nbins = length(ids))
    frequencies <- table(counts)
    data.frame(sex = sex, offspring = as.integer(names(frequencies)),
               n_parents = as.integer(frequencies))
  })))
}

#' Calculate realized effective parent numbers from full-population family sizes
#'
#' For each sex, k_i is the number of offspring assigned to adult i, including
#' adults with zero offspring, and K = sum(k_i). Two distinct offspring drawn
#' uniformly without replacement share a parent with probability
#' h = sum(k_i * (k_i - 1)) / (K * (K - 1)). Effective N is 1/h and the ratio
#' divides this by the census adult count. Calculations use the retained
#' full-population frequency table, never survey samples or expected sampling
#' weights. Counts are promoted to doubles before multiplication to avoid R
#' integer overflow. No individual records are expanded or regenerated.
#'
#' When K >= 2 but no sibling pairs exist, h = 0 and Effective N is Inf. Effective
#' N can exceed census N under this finite-cohort definition and is not capped.
#' When K < 2, the probability and effective numbers are undefined (NA).
#' @param replicate A saved replicate with reproduction_summary containing sex,
#' offspring, and n_parents for the entire population.
#' @return A tibble with sex, census_n, n_offspring, same_parent_probability,
#' effective_n, and effective_ratio, separately for females and males.
#' @export
pop_effective_numbers <- function(replicate) {
  counts <- replicate$reproduction_summary
  tibble::as_tibble(do.call(rbind, lapply(c("Female", "Male"), function(sex) {
    rows <- counts$sex == sex
    k <- as.double(counts$offspring[rows])
    frequency <- as.double(counts$n_parents[rows])
    census <- sum(frequency)
    total <- sum(k * frequency)
    h <- if (total >= 2) sum(k * (k - 1) * frequency) / (total * (total - 1)) else NA_real_
    effective <- if (is.na(h)) NA_real_ else if (h == 0) Inf else 1 / h
    data.frame(sex = sex, census_n = census, n_offspring = total,
               same_parent_probability = h, effective_n = effective,
               effective_ratio = if (census > 0) effective / census else NA_real_)
  })))
}

#' Derive stable, distinct seeds for the stages of a replicate
#'
#' Stage offsets are 0, 104729, and 209759, wrapping into 1..integer.max.
#' Derivation uses double arithmetic to avoid integer overflow. The resulting
#' seeds are saved explicitly so regeneration does not depend on re-deriving
#' them, and changing sampling cannot affect reproduction draws.
#' @param seed Positive integer replicate seed.
#' @return A named list of population, reproduction, and sampling seeds.
#' @keywords internal
.pop_stage_seeds <- function(seed) {
  .pop_integer(seed, "seed")
  offsets <- c(population = 0, reproduction = 104729, sampling = 209759)
  as.list(stats::setNames(as.integer((as.double(seed) - 1 + offsets) %%
                                      .Machine$integer.max + 1), names(offsets)))
}

#' Generate a replicate from explicit stage seeds and compact its retained data
#'
#' Full population and pedigree objects exist only within this call unless
#' retain_population is TRUE. No closures or attributes capture discarded rows.
#' Samples retain every sampled adult attribute and the sampled kids' complete
#' pedigree. Summaries count all adults, including parents with zero offspring.
#' @param settings Validated population and sampling settings.
#' @param seed Replicate identifier seed.
#' @param seeds Named stage seeds saved with the result.
#' @param retain_population Whether to additionally retain all population rows.
#' @return A schema-version-2 replicate object.
#' @keywords internal
.pop_build_rep <- function(settings, seed, seeds, retain_population) {
  if (!is.logical(retain_population) || length(retain_population) != 1L ||
      is.na(retain_population)) stop("retain_population must be TRUE or FALSE.")
  for (name in c("population", "reproduction", "sampling")) {
    .pop_integer(seeds[[name]], paste(name, "seed"))
  }
  population <- .pop_with_seed(seeds$population,
    if (!is.null(settings$size_model)) pop_generate_size_population(settings$Nf, settings$Nm, settings$size_model)
    else pop_generate_population(settings$Nf, settings$Nm, settings$use_weights,
                             settings$female_shape, settings$male_shape))
  pedigree <- .pop_with_seed(seeds$reproduction, pop_reproduce(population, settings$Nk))
  samples <- .pop_with_seed(seeds$sampling,
                            pop_sample(population, pedigree, settings$Sf, settings$Sm, settings$Sk))
  pairs <- pop_identify_pairs(samples)
  result <- list(schema_version = 2L, seed = as.integer(seed), seeds = seeds,
                 settings = settings, storage = if (retain_population) "full" else "compact",
                 reproduction_summary = pop_summarize_reproduction(population, pedigree),
                 sampling_weight_summary = if (settings$use_weights || !is.null(settings$size_model)) pop_summarize_sampling_weights(population) else NULL,
                 samples = samples, pairs = pairs,
                 estimates = pop_estimate(pairs, settings$Sf, settings$Sm, settings$Sk))
  if (!is.null(settings$size_model)) {
    result$size_summary <- pop_summarize_size(population)
    # A separate display stream cannot change traits, reproduction, or sampling.
    display_seed <- as.integer((as.double(seed) + 314158) %% .Machine$integer.max + 1)
    result$display_adults <- .pop_with_seed(display_seed, do.call(rbind,
      lapply(c("Female", "Male"), function(s) {
        adults <- population[population$sex == s, ]
        adults[sample.int(nrow(adults), min(500L, nrow(adults))), ]
      })))
  }
  if (retain_population) {
    result$population <- population
    result$pedigree <- pedigree
  }
  result
}

#' Simulate a replicate, retaining complete samples and population summaries
#'
#' Population generation, reproduction, and sampling use separate saved seeds.
#' Unsampled individuals are discarded after summarizing by default. The full
#' population can be regenerated on demand with pop_regenerate_rep().
#' @param settings Validated settings from pop_settings().
#' @param seed Replicate seed used to derive three stage seeds. RNG state is preserved.
#' @param retain_population Also retain the full adult population and pedigree.
#' Defaults to FALSE. This affects storage, not simulation results.
#' @return A list with settings, stage seeds, reproduction_summary, complete
#' samples, pairs, and estimates; optionally population and pedigree as well.
#' @export
pop_simulate_rep <- function(settings, seed, retain_population = FALSE) {
  settings <- do.call(pop_settings, settings)
  .pop_build_rep(settings, seed, .pop_stage_seeds(seed), retain_population)
}

#' Regenerate or resample a saved replicate from its stage seeds
#'
#' With defaults this exactly reproduces a compact replicate. Changing only
#' sampling_seed or sample sizes preserves its population and reproduction.
#' Set retain_population = TRUE to inspect or work with all individuals in R.
#' Regeneration assumes compatible simulation code and R/RNG behavior; trials
#' record schema, R version, and RNG kinds. Earlier schema objects are rejected.
#' @param replicate Saved schema-version-2 replicate.
#' @param Sf,Sm,Sk Sample sizes; default to the saved sample sizes.
#' @param sampling_seed Optional replacement seed for sampling only.
#' @param retain_population Whether to retain the regenerated full population.
#' @return A newly generated replicate, leaving the saved object unchanged.
#' @export
pop_regenerate_rep <- function(replicate, Sf = replicate$settings$Sf,
                               Sm = replicate$settings$Sm, Sk = replicate$settings$Sk,
                               sampling_seed = NULL, retain_population = FALSE) {
  if (!identical(replicate$schema_version, 2L)) {
    stop("Regeneration requires a schema-version-2 replicate.", call. = FALSE)
  }
  settings <- replicate$settings
  settings[c("Sf", "Sm", "Sk")] <- list(Sf, Sm, Sk)
  settings <- do.call(pop_settings, settings)
  seeds <- replicate$seeds
  if (!is.null(sampling_seed)) seeds$sampling <- sampling_seed
  .pop_build_rep(settings, replicate$seed, seeds, retain_population)
}

#' Run an immutable trial with independently generated populations
#'
#' Replicate r uses seed + r - 1. Complete samples and population frequency
#' summaries are retained; full populations can be regenerated from seeds.
#' Limits on memory/storage are
#' applied by the app, not this reusable function. Seeds reproduce draws under
#' the recorded RNG convention; record R/package versions for long-term archives.
#' @param settings Settings from pop_settings().
#' @param num_reps Number of independent population replicates.
#' @param seed Starting seed; the last replicate seed must fit in an R integer.
#' @param trial_id Positive integer identifier.
#' @param progress Optional function called with (completed, total) after each rep.
#' @return A list with metadata, settings, compact reps, and a replicate-summary tibble.
#' @export
pop_run_trial <- function(settings, num_reps = 100, seed = 1001, trial_id = 1,
                          progress = NULL) {
  settings <- do.call(pop_settings, settings)
  .pop_integer(num_reps, "num_reps")
  .pop_integer(seed, "seed")
  .pop_integer(trial_id, "trial_id")
  if (as.double(seed) + num_reps - 1 > .Machine$integer.max) {
    stop("The last replicate seed exceeds the integer limit.", call. = FALSE)
  }
  reps <- lapply(seq_len(num_reps), function(r) {
    ans <- pop_simulate_rep(settings, as.integer(as.double(seed) + r - 1))
    if (!is.null(progress)) progress(r, num_reps)
    ans
  })
  summary <- do.call(rbind, lapply(seq_along(reps), function(r) {
    x <- reps[[r]]$estimates
    x$rep_id <- r
    x$seed <- reps[[r]]$seed
    x$truth <- ifelse(x$sex == "Female", settings$Nf, settings$Nm)
    x
  }))
  list(schema_version = 2L, trial_id = as.integer(trial_id),
       scenario = if (!is.null(settings$size_model)) "Size-based reproduction" else
         if (settings$use_weights) "Randomly unequal reproduction" else "Equal reproduction",
       num_reps = as.integer(num_reps), seed = as.integer(seed), settings = settings,
       created_at = Sys.time(), R_version = R.version.string,
       rng = c("Mersenne-Twister", "Inversion", "Rejection"),
       reps = reps, summary = tibble::as_tibble(summary))
}

#' Compute approximate 95 percent Wald intervals for abundance estimates
#'
#' With C = sampled parents times sampled offspring and m = observed matches,
#' the Poisson working model has expected information C/N^3. Substituting the
#' estimate C/m gives variance Nhat^2/m and endpoints Nhat +/- 1.96 * SE.
#' Negative lower limits are deliberately preserved. Zero-match estimates have
#' no Wald interval, represented by NA variance, SE, and endpoints. This is an
#' asymptotic approximation, especially unreliable for small match counts.
#' @param estimates Tibble with estimate and matches columns, as returned by
#' pop_estimate() or the summary component of pop_run_trial().
#' @return The input tibble with variance, se, lower, upper, and interval_available
#' columns appended. This does not mutate the input or require population data.
#' @export
pop_wald_intervals <- function(estimates) {
  available <- is.finite(estimates$estimate) & estimates$matches > 0
  se <- rep(NA_real_, nrow(estimates))
  se[available] <- estimates$estimate[available] / sqrt(estimates$matches[available])
  estimates$variance <- se^2
  estimates$se <- se
  estimates$lower <- estimates$estimate - 1.96 * se
  estimates$upper <- estimates$estimate + 1.96 * se
  estimates$interval_available <- available
  estimates
}

#' Summarize abundance estimates separately by sex and method
#'
#' Median, mean, and CV exclude nonfinite estimates. CV is the sample standard
#' deviation (denominator n - 1) divided by the finite mean, returned as a ratio.
#' With no finite estimates, all three summaries are NA; with one finite estimate
#' the CV is NA. The CV is also undefined for a zero finite mean. Counts retain
#' all replicates so finite-only summaries can be interpreted conditionally.
#' CI inclusion is the fraction of finite replicates whose 95 percent Wald
#' interval contains the true abundance, including equality at either endpoint.
#' Its denominator excludes unbounded estimates; with no finite replicates the
#' fraction is NA. This is conditional coverage, not coverage across all reps.
#' @param estimates A replicate-summary tibble with sex, method, truth, and
#' estimate and matches columns, such as the summary element of pop_run_trial().
#' @return A tibble with grouping columns, finite_median, finite_mean, finite_cv,
#' num_finite, num_unbounded, num_reps, num_ci_included, and finite_ci_inclusion.
#' Abundances are not rounded.
#' @export
pop_summarize_estimates <- function(estimates) {
  keys <- unique(estimates[c("sex", "method", "truth")])
  if (!nrow(keys)) return(tibble::tibble())
  intervals <- pop_wald_intervals(estimates)
  included <- intervals$interval_available & intervals$lower <= intervals$truth &
    intervals$upper >= intervals$truth
  tibble::as_tibble(do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    rows <- estimates$sex == key$sex & estimates$method == key$method &
      estimates$truth == key$truth
    values <- estimates$estimate[rows]
    finite <- values[is.finite(values)]
    n <- length(finite)
    average <- if (n) mean(finite) else NA_real_
    cbind(key, data.frame(
      finite_median = if (n) stats::median(finite) else NA_real_,
      finite_mean = average,
      finite_cv = if (n >= 2 && average != 0) stats::sd(finite) / average else NA_real_,
      num_finite = n, num_unbounded = sum(!is.finite(values)), num_reps = length(values),
      num_ci_included = sum(included[rows], na.rm = TRUE),
      finite_ci_inclusion = if (n) sum(included[rows], na.rm = TRUE) / n else NA_real_
    ))
  })))
}

#' Make the scenario-one trial register
#'
#' The register is a compact view; trial objects retain samples and summaries.
#' Other scenarios can supply their own register views over their trial objects.
#' @param trials List of saved trial objects.
#' @param include_weights Include the random-reproduction checkbox and gamma
#' shapes in the register, as used by the second tab. Defaults to FALSE.
#' @return A tibble with settings, seed range, and match/estimation summaries.
#' @export
pop_trial_table <- function(trials, include_weights = FALSE) {
  if (!length(trials)) return(tibble::tibble())
  tibble::as_tibble(do.call(rbind, lapply(trials, function(trial) {
    x <- trial$summary[trial$summary$method == "Moments", ]
    female <- x$sex == "Female"
    row <- data.frame(trial_id = trial$trial_id, num_reps = trial$num_reps,
               as.data.frame(trial$settings[c("Nf", "Nm", "Nk", "Sf", "Sm", "Sk")]),
               samples_retained = sum(unlist(trial$settings[c("Sf", "Sm", "Sk")])) * trial$num_reps,
               seed_first = trial$seed,
               seed_last = as.integer(as.double(trial$seed) + trial$num_reps - 1),
               mean_MOP = mean(x$matches[female]), mean_FOP = mean(x$matches[!female]),
               unbounded_female = sum(!x$finite[female]),
               unbounded_male = sum(!x$finite[!female]))
    if (include_weights) {
      row$use_sampling_weights <- trial$settings$use_weights
      row$female_gamma_shape <- trial$settings$female_shape
      row$male_gamma_shape <- trial$settings$male_shape
      row <- row[c("trial_id", "use_sampling_weights",
                   setdiff(names(row), c("trial_id", "use_sampling_weights")))]
    }
    row
  })))
}

#' Retain compact population summaries of reproductive sampling weights
#' @param population Complete adult population from pop_generate_population().
#' @return A tibble of sex-specific counts, raw means, realized CVs, ranges,
#' and 1st/50th/99th percentiles. No full adult vectors are captured.
#' @export
pop_summarize_sampling_weights <- function(population) {
  tibble::as_tibble(do.call(rbind, lapply(c("Female", "Male"), function(sex) {
    w <- population$reproductive_weight[population$sex == sex]
    q <- stats::quantile(w, c(0.01, 0.5, 0.99), names = FALSE)
    data.frame(sex = sex, n_adults = length(w), mean = mean(w),
               cv = if (length(w) > 1) stats::sd(w) / mean(w) else NA_real_,
               minimum = min(w), q01 = q[1], median = q[2], q99 = q[3], maximum = max(w))
  })))
}

#' Specify a sex-specific growth, survival, and reproduction model
#'
#' Ages follow a shifted exponential with rate Z, conditioned on age <= max_age.
#' Linf and K are correlated lognormal traits with specified arithmetic means
#' and CVs; rho is their correlation on the log scale. Length follows
#' Linf * (1 - exp(-K * (age - t0))). Reproductive propensity is proportional
#' to length^exponent times independent mean-one lognormal noise with CV output_cv.
#' ERRO is this propensity divided by its realized within-sex population mean.
#' @param Linf,K Arithmetic means of individual growth parameters.
#' @param t0 Age at zero length, strictly below min_age.
#' @param Z Continuous mortality rate for the adult age distribution.
#' @param min_age,max_age Youngest and oldest adult ages (truncated, not capped).
#' @param Linf_cv,K_cv Coefficients of variation of growth parameters.
#' @param rho Correlation of log(Linf) and log(K), between -1 and 1.
#' @param exponent Power of length determining expected reproductive output.
#' @param output_cv Conditional reproductive propensity CV at a given length.
#' @return A validated parameter list; lengths and ages use arbitrary consistent units.
#' @export
pop_size_parameters <- function(Linf = 100, K = 0.2, t0 = 0, Z = 0.2,
                                min_age = 1, max_age = 30, Linf_cv = 0.1,
                                K_cv = 0.15, rho = -0.5, exponent = 3,
                                output_cv = 0.3) {
  p <- as.list(environment())
  bounds <- list(Linf = c(1, 500), K = c(0.01, 3), t0 = c(-20, 100),
                 Z = c(0.01, 5), min_age = c(0, 100), max_age = c(0.01, 200),
                 Linf_cv = c(0, 3), K_cv = c(0, 3), rho = c(-1, 1),
                 exponent = c(0, 8), output_cv = c(0, 3))
  for (name in names(p)) {
    x <- p[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < bounds[[name]][1] || x > bounds[[name]][2])
      stop(name, " must be between ", bounds[[name]][1], " and ", bounds[[name]][2], ".", call. = FALSE)
  }
  if (max_age <= min_age || t0 >= min_age)
    stop("Require t0 < minimum adult age < maximum adult age.", call. = FALSE)
  p
}

#' Draw adult traits and positive, mean-one relative reproductive propensities
#'
#' Uses inverse-CDF sampling for the truncated exponential (no pile-up at the
#' maximum age). Log propensities are centered before exponentiation for numerical
#' stability. No survey or offspring draws occur here.
#' @param n Number of adults.
#' @param parameters Validated sex-specific growth parameters.
#' @return A data frame with age, individual Linf/K, length, and reproductive_weight.
#' @export
pop_size_traits <- function(n, parameters = pop_size_parameters()) {
  .pop_integer(n, "n")
  p <- do.call(pop_size_parameters, parameters)
  age <- p$min_age - log1p(-stats::runif(n) * (-expm1(-p$Z * (p$max_age - p$min_age)))) / p$Z
  z1 <- stats::rnorm(n)
  z2 <- p$rho * z1 + sqrt(max(0, 1 - p$rho^2)) * stats::rnorm(n)
  sL <- sqrt(log1p(p$Linf_cv^2))
  sK <- sqrt(log1p(p$K_cv^2))
  Linf <- p$Linf * exp(sL * z1 - sL^2 / 2)
  K <- p$K * exp(sK * z2 - sK^2 / 2)
  length <- Linf * (-expm1(-K * (age - p$t0)))
  sr <- sqrt(log1p(p$output_cv^2))
  log_output <- p$exponent * log(length) + sr * stats::rnorm(n) - sr^2 / 2
  relative <- exp(log_output - max(log_output))
  data.frame(age = age, individual_Linf = Linf, individual_K = K, length = length,
             expected_reproductive_weight = exp(p$exponent * log(length) - max(log_output)) / mean(relative),
             reproductive_weight = relative / mean(relative))
}

#' Generate both sexes with length-based reproductive sampling weights
#' @param Nf,Nm Numbers of female and male adults.
#' @param size_model Named Female and Male parameter lists.
#' @return Adult tibble retaining individual growth, length, and ERRO attributes.
#' @export
pop_generate_size_population <- function(Nf, Nm, size_model) {
  tibble::as_tibble(do.call(rbind, lapply(c("Female", "Male"), function(s) {
    n <- if (s == "Female") Nf else Nm
    cbind(data.frame(parent_id = seq_len(n), sex = s), pop_size_traits(n, size_model[[s]]))
  })))
}

#' Retain full-population trait histograms and offspring-parent expectations
#'
#' Each sex and variable uses 60 equal-width bins covering every adult. Counts
#' summarize f(x); summed ERRO in each bin summarizes the distribution of parents
#' of randomly chosen offspring. For ERRO this is the bin integral of r*f(r),
#' using realized mean-one ERRO, with no kernel approximation or tail truncation.
#' @param population Complete size-based adult population.
#' @return A compact tibble with bin endpoints, adult counts, ERRO sums, and sex.
#' @export
pop_summarize_size <- function(population) {
  tibble::as_tibble(do.call(rbind, lapply(c("Female", "Male"), function(s) {
    adults <- population[population$sex == s, ]
    do.call(rbind, lapply(c("length", "reproductive_weight"), function(variable) {
      x <- adults[[variable]]
      bounds <- range(c(0, x))
      breaks <- seq(bounds[1], max(bounds[2], 1e-8) * (1 + 1e-10), length.out = 61)
      bin <- findInterval(x, breaks, all.inside = TRUE)
      weighted <- rowsum(adults$reproductive_weight, bin, reorder = FALSE)
      totals <- numeric(60)
      totals[as.integer(rownames(weighted))] <- weighted[, 1]
      data.frame(sex = s, variable = variable, lower = head(breaks, -1),
                 upper = tail(breaks, -1), count = tabulate(bin, nbins = 60),
                 erro_sum = totals)
    }))
  })))
}

#' Flatten size-based settings into a scenario-specific trial register
#' @param trials Saved size-based trials.
#' @return The shared register plus every sex-specific growth/reproduction setting.
#' @export
pop_size_trial_table <- function(trials) {
  result <- pop_trial_table(trials)
  if (!length(trials)) return(result)
  for (sex in c("Female", "Male")) for (name in names(pop_size_parameters())) {
    result[[paste(tolower(sex), name, sep = "_")]] <- vapply(trials,
      function(t) t$settings$size_model[[sex]][[name]], numeric(1))
  }
  result
}
