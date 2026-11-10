# Sibling experiments reuse the POP adult, reproduction, RNG, and summary code.

#' Validate a within-cohort or cross-cohort sibling experiment
#' @param Nf,Nm Numbers of adult females and males, constant between episodes.
#' @param Nk1,Nk2 Offspring population sizes in episodes one and two.
#' @param Sk1,Sk2 Uniform offspring sample sizes without replacement.
#' @param scenario One of equal, lucky, cross_lucky, or size.
#' @param use_weights Enable gamma reproductive sampling weights in lucky scenarios.
#' @param female_shape,male_shape Gamma shapes, with scale fixed at one.
#' @param size_model Female and Male parameter lists from pop_size_parameters().
#' In the size scenario output_cv is forced to zero: size effects persist and
#' there is no additional lucky-litter effect independent of length.
#' @param growth Advance size-based adults by one year between episodes.
#' @param growth_cv CV of a mean-one lognormal multiplier on each adult's annual
#' growth rate. Zero gives its deterministic von Bertalanffy increment.
#' @return Validated settings; unused second-episode sizes are stored as zero.
#' @export
hsp_settings <- function(Nf = 1000, Nm = 1000, Nk1 = 2000, Nk2 = 2000,
                         Sk1 = 200, Sk2 = 200, scenario = c("equal", "lucky", "cross_lucky", "size"),
                         use_weights = TRUE, female_shape = 1, male_shape = 1, size_model = NULL,
                         growth = TRUE, growth_cv = .25) {
  scenario <- match.arg(scenario)
  cross <- scenario %in% c("cross_lucky", "size")
  base <- pop_settings(Nf, Nm, Nk1, Sf = 1, Sm = 1, Sk = Sk1,
    use_weights = use_weights, female_shape = female_shape, male_shape = male_shape)
  if (!cross && Sk1 < 2) stop("Within-cohort comparisons need at least two sampled offspring.")
  if (cross) {
    .pop_integer(Nk2, "Episode 2 offspring")
    .pop_integer(Sk2, "Episode 2 sample")
    if (Sk2 > Nk2) stop("Episode 2 sample cannot exceed its offspring population.")
  }
  if (scenario == "size") {
    if (!is.logical(growth) || length(growth) != 1L || is.na(growth)) stop("growth must be TRUE or FALSE.")
    if (!is.numeric(growth_cv) || length(growth_cv) != 1L || !is.finite(growth_cv) || growth_cv < 0 || growth_cv > 3)
      stop("Growth-rate CV must be between zero and 3.")
    if (is.null(size_model)) size_model <- list(Female = pop_size_parameters(exponent = 3),
                                               Male = pop_size_parameters(exponent = 1))
    size_model <- pop_settings(Nf, Nm, Nk1, Sf = 1, Sm = 1, Sk = Sk1, size_model = size_model)$size_model
    size_model <- lapply(size_model, function(p) { p$output_cv <- 0; p })
  } else size_model <- NULL
  list(Nf = base$Nf, Nm = base$Nm, Nk1 = base$Nk, Nk2 = if (cross) as.integer(Nk2) else 0L,
    Sk1 = base$Sk, Sk2 = if (cross) as.integer(Sk2) else 0L, scenario = scenario, cross = cross,
    use_weights = scenario %in% c("lucky", "cross_lucky") && use_weights,
    female_shape = female_shape, male_shape = male_shape, size_model = size_model,
    growth = scenario == "size" && growth, growth_cv = if (scenario == "size") growth_cv else 0)
}

#' Grow size-based adults for one year without shrinking or crossing their Linf
#'
#' Given remaining growth D = max(0, Linf - length), the increment is
#' D * (1 - exp(-individual_K * multiplier)). The positive mean-one lognormal
#' multiplier varies independently among adults. This is variation in growth
#' rate, not a specified CV of the resulting increment. Individual Linf/K remain
#' fixed. Ages advance by one even beyond the initial age-distribution cutoff;
#' no mortality is applied. Relative reproductive sampling weights are recalculated
#' from new lengths and normalized separately within sex.
#' @param adults Population from pop_generate_size_population().
#' @param size_model Female and Male length-exponent parameter lists.
#' @param growth_cv CV of the annual growth-rate multiplier, between zero and 3.
#' @return Adult tibble with unchanged IDs and growth traits, advanced ages and
#' nondecreasing lengths, and updated mean-one length-based sampling weights.
#' @export
hsp_grow_adults <- function(adults, size_model, growth_cv = .25) {
  stopifnot(length(growth_cv) == 1L, is.finite(growth_cv), growth_cv >= 0, growth_cv <= 3)
  sigma <- sqrt(log1p(growth_cv^2))
  multiplier <- if (growth_cv == 0) rep(1, nrow(adults)) else
    exp(sigma * stats::rnorm(nrow(adults)) - sigma^2 / 2)
  remaining <- pmax(0, adults$individual_Linf - adults$length)
  adults$length <- adults$length + remaining * (-expm1(-adults$individual_K * multiplier))
  adults$age <- adults$age + 1
  for (sex in c("Female", "Male")) {
    rows <- which(adults$sex == sex)
    log_output <- size_model[[sex]]$exponent * log(adults$length[rows])
    weight <- exp(log_output - max(log_output))
    adults$reproductive_weight[rows] <- weight / mean(weight)
    adults$expected_reproductive_weight[rows] <- adults$reproductive_weight[rows]
  }
  adults
}

#' Count sibling comparisons without constructing all individual pairs
#' @param first,second Offspring pedigree tibbles, with mother_id and father_id.
#' Supply second only for cross-cohort comparisons. IDs may repeat between cohorts.
#' @param Nf,Nm Adult population sizes; IDs within sex run from one to census N.
#' @return List with family counts, sex-specific match counts, full-sibling pair
#' count, and number of distinct offspring comparisons. A full sibling contributes
#' once to each sex's match count. All arithmetic products use doubles.
#' @export
hsp_identify_matches <- function(first, Nf, Nm, second = NULL) {
  cross <- !is.null(second)
  families <- lapply(c("Female", "Male"), function(sex) {
    key <- if (sex == "Female") "mother_id" else "father_id"
    n <- if (sex == "Female") Nf else Nm
    k1 <- as.double(tabulate(first[[key]], nbins = n))
    k2 <- if (cross) as.double(tabulate(second[[key]], nbins = n)) else rep(0, n)
    matches <- if (cross) k1 * k2 else k1 * (k1 - 1) / 2
    data.frame(sex = sex, parent_id = seq_len(n), k1, k2, matches)
  })
  # Full siblings share the parental couple; this is only a diagnostic count.
  couples <- function(x) table(paste(x$mother_id, x$father_id, sep = ":"))
  a <- couples(first)
  if (cross) {
    b <- couples(second)
    common <- intersect(names(a), names(b))
    full <- sum(as.double(a[common]) * as.double(b[common]))
  } else full <- sum(as.double(a) * (as.double(a) - 1) / 2)
  list(families = tibble::as_tibble(do.call(rbind, families)),
       matches = vapply(families, function(x) sum(x$matches), numeric(1)),
       full_siblings = full,
       comparisons = if (cross) as.double(nrow(first)) * nrow(second) else choose(nrow(first), 2))
}

#' Estimate adult abundance from maternal and paternal sibling matches
#' @param matches Two counts, female then male; full siblings count in both.
#' @param comparisons Number of unordered within-cohort or cross-cohort comparisons.
#' @return Same estimate schema as pop_estimate(), allowing shared plots and
#' Wald summaries. Both estimators use C/m; zero matches give infinity.
#' @export
hsp_estimate <- function(matches, comparisons) {
  stopifnot(length(matches) == 2, all(is.finite(matches)), all(matches >= 0),
    length(comparisons) == 1, is.finite(comparisons), comparisons > 0, all(matches <= comparisons))
  estimates <- comparisons / matches
  tibble::tibble(sex = rep(c("Female", "Male"), 2),
    method = rep(c("Moments", "Poisson pseudo-MLE"), each = 2), matches = rep(matches, 2),
    estimate = rep(estimates, 2), finite = is.finite(rep(estimates, 2)))
}

#' Make a compact histogram of family counts including zero-offspring adults
#' @param counts Numeric offspring counts for every adult of one sex.
#' @param sex Female or Male.
#' @return Reproduction frequency data compatible with the POP plotting language.
#' @keywords internal
.hsp_frequency <- function(counts, sex) {
  frequencies <- table(counts)
  data.frame(sex = sex, offspring = as.double(names(frequencies)), n_parents = as.double(frequencies))
}

#' Simulate one sibling replicate with fixed adult identities across episodes
#' @param settings Validated hsp_settings().
#' @param seed Positive replicate seed; population, episode luck, reproduction,
#' sampling, and display draws use saved, distinct deterministic streams.
#' @param retain_population Retain both complete populations and pedigrees for checks.
#' @return Compact replicate with complete offspring samples, family summaries,
#' effective numbers, match counts, display adults, and estimates. Census adults
#' and unsampled offspring are discarded by default. Seeds permit regeneration.
#' @export
hsp_simulate_rep <- function(settings, seed, retain_population = FALSE) {
  # Remove the derived flag before revalidating public settings.
  settings$cross <- NULL
  s <- do.call(hsp_settings, settings)
  .pop_integer(seed, "seed")
  offsets <- c(population = 0, luck_2 = 524287, reproduction_1 = 104729,
               reproduction_2 = 314159, sampling_1 = 209759, sampling_2 = 734003, display = 999983,
               growth = 15485863)
  seeds <- as.list(as.integer((as.double(seed) - 1 + offsets) %% .Machine$integer.max + 1))
  names(seeds) <- names(offsets)
  adults1 <- .pop_with_seed(seeds$population,
    if (s$scenario == "size") pop_generate_size_population(s$Nf, s$Nm, s$size_model)
    else pop_generate_population(s$Nf, s$Nm, s$use_weights, s$female_shape, s$male_shape))
  # Identical IDs, independent seasonal luck. Size-based adults can grow, while
  # keeping their individual Linf/K and the link between length and fecundity.
  adults2 <- if (s$scenario == "cross_lucky") .pop_with_seed(seeds$luck_2,
    pop_generate_population(s$Nf, s$Nm, s$use_weights, s$female_shape, s$male_shape)) else
      if (s$scenario == "size" && s$growth) .pop_with_seed(seeds$growth,
        hsp_grow_adults(adults1, s$size_model, s$growth_cv)) else adults1
  pedigree1 <- .pop_with_seed(seeds$reproduction_1, pop_reproduce(adults1, s$Nk1))
  pedigree2 <- if (s$cross) .pop_with_seed(seeds$reproduction_2, pop_reproduce(adults2, s$Nk2)) else NULL
  sample1 <- .pop_with_seed(seeds$sampling_1, pedigree1[sample.int(s$Nk1, s$Sk1), ])
  sample2 <- if (s$cross) .pop_with_seed(seeds$sampling_2, pedigree2[sample.int(s$Nk2, s$Sk2), ]) else NULL
  found <- hsp_identify_matches(sample1, s$Nf, s$Nm, sample2)
  population_summaries <- sample_summaries <- effective <- displays <- list()
  for (i in 1:2) {
    sex <- c("Female", "Male")[i]
    a <- adults1[adults1$sex == sex, ]
    b <- adults2[adults2$sex == sex, ]
    n <- nrow(a)
    key <- if (sex == "Female") "mother_id" else "father_id"
    k1 <- as.double(tabulate(pedigree1[[key]], nbins = n))
    k2 <- if (s$cross) as.double(tabulate(pedigree2[[key]], nbins = n)) else rep(0, n)
    f <- found$families[found$families$sex == sex, ]
    q1 <- a$reproductive_weight / sum(a$reproductive_weight)
    q2 <- b$reproductive_weight / sum(b$reproductive_weight)
    ps <- list(); ss <- list()
    for (episode in seq_len(if (s$cross) 2L else 1L)) {
      ps[[episode]] <- cbind(.hsp_frequency(if (episode == 1) k1 else k2, sex), episode = episode)
      ss[[episode]] <- cbind(.hsp_frequency(if (episode == 1) f$k1 else f$k2, sex), episode = episode)
    }
    population_summaries[[i]] <- do.call(rbind, ps)
    sample_summaries[[i]] <- do.call(rbind, ss)
    realized <- if (s$cross) sum(k1 * k2) / (as.double(s$Nk1) * s$Nk2) else
      if (s$Nk1 >= 2) sum(k1 * (k1 - 1)) / (as.double(s$Nk1) * (s$Nk1 - 1)) else NA_real_
    propensity <- if (s$cross) sum(q1 * q2) else sum(q1^2)
    effective[[i]] <- data.frame(sex = sex, census_n = n, same_parent_probability = realized,
      effective_n = 1 / realized, effective_ratio = 1 / (realized * n),
      weight_probability = propensity, expected_sample_matches = found$comparisons * propensity)
    # Bounded diagnostic subset, independent of all simulation/sample draws.
    idx <- .pop_with_seed(as.integer((as.double(seeds$display) + i - 1) %% .Machine$integer.max + 1),
                          sample.int(n, min(500L, n)))
    shown <- a[idx, ]
    shown$weight_1 <- n * q1[idx]; shown$weight_2 <- n * q2[idx]
    if (s$scenario == "size") {
      shown$length_2 <- b$length[idx]; shown$age_2 <- b$age[idx]
    }
    shown$offspring_1 <- k1[idx]; shown$offspring_2 <- k2[idx]
    displays[[i]] <- shown
  }
  result <- list(schema_version = 2L, seed = as.integer(seed), seeds = seeds, settings = s,
    samples = list(cohort_1 = sample1, cohort_2 = sample2),
    # Families with no sampled offspring need no individual row in the archive.
    families = found$families[found$families$k1 + found$families$k2 > 0, ],
    matches = stats::setNames(found$matches, c("Female", "Male")),
    full_siblings = found$full_siblings, comparisons = found$comparisons,
    reproduction_summary = tibble::as_tibble(do.call(rbind, population_summaries)),
    sampling_summary = tibble::as_tibble(do.call(rbind, sample_summaries)),
    effective = tibble::as_tibble(do.call(rbind, effective)),
    display_adults = tibble::as_tibble(do.call(rbind, displays)),
    estimates = hsp_estimate(found$matches, found$comparisons))
  if (s$scenario == "size") {
    traits <- pop_summarize_size(adults1)
    traits$sample_sum <- traits$match_sum <- traits$sibling_weight_sum <- traits$sample_weight_sum <- 0
    for (sex in c("Female", "Male")) for (variable in c("length", "reproductive_weight")) {
      a <- adults1[adults1$sex == sex, ]
      b <- adults2[adults2$sex == sex, ]
      f <- found$families[found$families$sex == sex, ]
      rows <- which(traits$sex == sex & traits$variable == variable)
      bins <- traits[rows, ]
      bin <- findInterval(a[[variable]], c(bins$lower, tail(bins$upper, 1)), all.inside = TRUE)
      for (column in c("sample_sum", "match_sum", "sibling_weight_sum", "sample_weight_sum")) {
        values <- switch(column, sample_sum = f$k1 + f$k2, match_sum = f$matches,
                          sibling_weight_sum = a$reproductive_weight * b$reproductive_weight,
                          sample_weight_sum = s$Sk1 * a$reproductive_weight + s$Sk2 * b$reproductive_weight)
        sums <- rowsum(values, bin, reorder = FALSE)
        traits[[column]][rows[as.integer(rownames(sums))]] <- sums[, 1]
      }
    }
    result$size_summary <- traits
    result$size_summary_2 <- pop_summarize_size(adults2)
  }
  if (retain_population) result <- c(result, list(population_1 = adults1, population_2 = adults2,
                                                 pedigree_1 = pedigree1, pedigree_2 = pedigree2))
  result
}

#' Reproduce a sibling replicate from its saved settings and seed
#' @param replicate Saved sibling replicate with schema_version 1 or 2.
#' @param retain_population Include complete adult and offspring populations.
#' @return Regenerated replicate; requires compatible code and R RNG behavior.
#' @export
hsp_regenerate_rep <- function(replicate, retain_population = FALSE) {
  if (!replicate$schema_version %in% c(1L, 2L)) stop("Unsupported sibling archive schema.")
  settings <- replicate$settings
  # Before growth was introduced, size effects were fixed. Replay those samples
  # without growth; the returned archive uses the current schema and summaries.
  if (replicate$schema_version == 1L) { settings$growth <- FALSE; settings$growth_cv <- 0 }
  hsp_simulate_rep(settings, replicate$seed, retain_population)
}

#' Run independent sibling replicates and retain a trial-level estimate table
#' @param settings Validated sibling settings.
#' @param num_reps Number of independently simulated populations.
#' @param seed First replicate seed; subsequent seeds increase by one.
#' @param trial_id Trial identifier in this scenario's registry.
#' @param progress Optional function(completed, total).
#' @return Trial with complete samples, compact summaries, settings and seeds.
#' @export
hsp_run_trial <- function(settings, num_reps = 100, seed = 1001, trial_id = 1, progress = NULL) {
  settings$cross <- NULL
  settings <- do.call(hsp_settings, settings)
  .pop_integer(num_reps, "Replicates"); .pop_integer(seed, "Seed"); .pop_integer(trial_id, "Trial ID")
  if (as.double(seed) + num_reps - 1 > .Machine$integer.max) stop("Last replicate seed exceeds the integer limit.")
  reps <- lapply(seq_len(num_reps), function(i) {
    result <- hsp_simulate_rep(settings, as.integer(as.double(seed) + i - 1))
    if (!is.null(progress)) progress(i, num_reps)
    result
  })
  summary <- do.call(rbind, lapply(seq_along(reps), function(i) {
    x <- reps[[i]]$estimates
    x$rep_id <- i; x$seed <- reps[[i]]$seed
    x$truth <- ifelse(x$sex == "Female", settings$Nf, settings$Nm)
    x
  }))
  list(schema_version = 2L, trial_id = as.integer(trial_id), scenario = settings$scenario,
    num_reps = as.integer(num_reps), seed = as.integer(seed), settings = settings,
    created_at = Sys.time(), R_version = R.version.string,
    rng = c("Mersenne-Twister", "Inversion", "Rejection"),
    reps = reps, summary = tibble::as_tibble(summary))
}

#' Tabulate saved sibling trials with scenario-specific settings
#' @param trials List of sibling trial objects.
#' @return Tibble with seeds, cohort sizes, sampling settings, and match means;
#' lucky and size scenarios also include their reproductive parameters.
#' @export
hsp_trial_table <- function(trials) {
  if (!length(trials)) return(tibble::tibble())
  tibble::as_tibble(do.call(rbind, lapply(trials, function(t) {
    s <- t$settings; x <- t$summary[t$summary$method == "Moments", ]
    row <- data.frame(trial_id = t$trial_id, num_reps = t$num_reps,
      Nf = s$Nf, Nm = s$Nm, Nk1 = s$Nk1, Nk2 = s$Nk2, Sk1 = s$Sk1, Sk2 = s$Sk2,
      seed_first = t$seed, seed_last = as.double(t$seed) + t$num_reps - 1,
      mean_maternal_matches = mean(x$matches[x$sex == "Female"]),
      mean_paternal_matches = mean(x$matches[x$sex == "Male"]),
      samples_retained = (s$Sk1 + s$Sk2) * t$num_reps)
    if (!s$cross) {
      row <- row[setdiff(names(row), c("Nk2", "Sk2"))]
      names(row)[names(row) == "Nk1"] <- "Nk"
      names(row)[names(row) == "Sk1"] <- "Sk"
    }
    if (s$scenario %in% c("lucky", "cross_lucky")) {
      row$use_sampling_weights <- s$use_weights
      row <- row[c("trial_id", "use_sampling_weights", setdiff(names(row), c("trial_id", "use_sampling_weights")))]
      row$female_gamma_shape <- s$female_shape; row$male_gamma_shape <- s$male_shape
    }
    if (s$scenario == "size") {
      row$one_year_growth <- isTRUE(s$growth)
      row$growth_rate_cv <- if (is.null(s$growth_cv)) 0 else s$growth_cv
      for (sex in c("Female", "Male")) for (name in names(s$size_model[[sex]]))
        row[[paste(tolower(sex), name, sep = "_")]] <- s$size_model[[sex]][[name]]
    }
    row
  })))
}
