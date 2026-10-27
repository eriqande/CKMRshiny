# Run from the repository root: Rscript tests/check-pop.R
# No testing framework is needed; Shiny server checks run if Shiny is installed.
source("R/pop-simulation.R")
source("R/pop-plots.R")

expect_error <- function(code) {
  result <- tryCatch({ force(code); FALSE }, error = function(e) TRUE)
  stopifnot(result)
}

# Reject malformed sizes, oversampling, and seeds that would overflow.
expect_error(pop_settings(Nf = 0))
expect_error(pop_settings(Sk = 1.5))
expect_error(pop_settings(Sf = NA_real_))
expect_error(pop_settings(Sf = 1001))
expect_error(pop_run_trial(pop_settings(), 2, .Machine$integer.max))

# Reproducibility includes restoration of the caller's RNG state and kinds.
set.seed(987)
before <- .Random.seed
kind <- RNGkind()
settings <- pop_settings(Nf = 17, Nm = 31, Nk = 80, Sf = 8, Sm = 11, Sk = 40)
a <- pop_simulate_rep(settings, 42)
b <- pop_simulate_rep(settings, 42)
stopifnot(identical(a, b), identical(before, .Random.seed), identical(kind, RNGkind()))
expect_error(.pop_with_seed(42, stop("deliberate failure")))
stopifnot(identical(before, .Random.seed))
rm(.Random.seed, envir = .GlobalEnv)
invisible(pop_simulate_rep(settings, 42))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

# Compact storage has no full population/pedigree; regeneration exactly restores
# summaries and samples, including at the integer seed boundary tested below.
stopifnot(is.null(a$population), is.null(a$pedigree), a$storage == "compact",
          identical(a, pop_regenerate_rep(a)), length(unique(unlist(a$seeds))) == 3)
full <- pop_regenerate_rep(a, retain_population = TRUE)
stopifnot(identical(a$reproduction_summary, full$reproduction_summary),
          identical(a$samples, full$samples), identical(a$pairs, full$pairs))
resampled <- pop_regenerate_rep(a, Sf = 10, sampling_seed = 9001,
                               retain_population = TRUE)
stopifnot(identical(full$population, resampled$population),
          identical(full$pedigree, resampled$pedigree),
          identical(full$reproduction_summary, resampled$reproduction_summary),
          nrow(resampled$samples$females) == 10,
          !identical(full$samples, resampled$samples),
          identical(a, pop_simulate_rep(settings, 42)))
expect_error(pop_regenerate_rep(list(schema_version = 1L)))
expect_error(pop_regenerate_rep(a, sampling_seed = 0))
expect_error(pop_regenerate_rep(a, Sf = 99))

# Unequal sex sizes expose parent-pool mixups; all samples are unique and in range.
stopifnot(all(full$pedigree$mother_id %in% seq_len(17)),
          all(full$pedigree$father_id %in% seq_len(31)),
          length(unique(a$samples$females$parent_id)) == 8,
          length(unique(a$samples$males$parent_id)) == 11,
          length(unique(a$samples$kids$kid_id)) == 40,
          identical(a$pairs, pop_identify_pairs(a$samples)))
counts <- pop_parent_counts(a)
observed <- pop_parent_counts(a, TRUE)
for (sex in c("Female", "Male")) {
  stopifnot(sum(counts$offspring[counts$sex == sex] * counts$n_parents[counts$sex == sex]) == 80,
            sum(counts$n_parents[counts$sex == sex]) == if (sex == "Female") 17 else 31,
            sum(observed$offspring[observed$sex == sex] * observed$n_parents[observed$sex == sex]) == sum(a$pairs$sex == sex),
            sum(observed$n_parents[observed$sex == sex]) == if (sex == "Female") 8 else 11)
}

# A census gives every pair, including the singleton sample.int edge case.
census <- pop_simulate_rep(pop_settings(1, 3, 20, 1, 3, 20), 5)
stopifnot(sum(census$pairs$sex == "Female") == 20,
          sum(census$pairs$sex == "Male") == 20,
          identical(census$estimates$estimate, c(1, 3, 1, 3)))
last_seed <- pop_run_trial(pop_settings(1, 1, 1, 1, 1, 1), 1, .Machine$integer.max)
stopifnot(last_seed$reps[[1]]$seed == .Machine$integer.max,
          pop_trial_table(list(last_seed))$seed_last == .Machine$integer.max,
          identical(last_seed$reps[[1]], pop_regenerate_rep(last_seed$reps[[1]])))

# A known sample with no observed parents must preserve zero counts and Inf.
pedigree <- tibble::tibble(kid_id = 1:3, mother_id = c(1L, 1L, 2L), father_id = c(2L, 2L, 1L))
known_sample <- list(females = tibble::tibble(parent_id = 3L),
                     males = tibble::tibble(parent_id = 3L), kids = pedigree)
empty <- pop_identify_pairs(known_sample)
stopifnot(nrow(empty) == 0, all(is.infinite(pop_estimate(empty, 1, 1, 3)$estimate)))
known_sample$females$parent_id <- 1L
known_sample$males$parent_id <- 2L
pairs <- pop_identify_pairs(known_sample)
stopifnot(sum(pairs$sex == "Female") == 2, sum(pairs$sex == "Male") == 2)

# Trials keep one summary per sex/method/rep, including zero-pair replicates.
sparse <- pop_run_trial(pop_settings(50, 70, 10, 1, 1, 1), 30, 301)
stopifnot(nrow(sparse$summary) == 120, length(sparse$reps) == 30,
          all(vapply(sparse$reps, function(x) is.list(x$samples), logical(1))),
          identical(sparse$reps[[3]], pop_simulate_rep(sparse$settings, 303)),
          all(sparse$summary$estimate[sparse$summary$method == "Moments"] ==
              sparse$summary$estimate[sparse$summary$method == "Poisson pseudo-MLE"]))
register <- pop_trial_table(list(sparse))
stopifnot(register$num_reps == 30, register$seed_first == 301,
          register$seed_last == 330, register$unbounded_female > 0,
          register$samples_retained == 90)
archive <- tempfile(fileext = ".rds")
saveRDS(sparse, archive)
stopifnot(identical(readRDS(archive), sparse))
unlink(archive)

# Larger populations stay compact, and population totals across reps no longer
# impose a storage cap. Compare retained sizes against the same full replicate.
large_settings <- pop_settings(50000, 50000, 200000, 20, 20, 40)
large <- pop_run_trial(large_settings, 10, 8001)
large_full <- pop_regenerate_rep(large$reps[[1]], retain_population = TRUE)
stopifnot(as.numeric(object.size(large$reps[[1]])) < as.numeric(object.size(large_full)) / 20,
          all(vapply(large$reps, function(x) is.null(x$population) && is.null(x$pedigree), logical(1))))
cat("Compact large-population trial:", format(object.size(large), units = "auto"),
    "retained for 3 million population individuals across reps.\n")

# Finite-only summaries exclude Inf and handle zero/one finite observation.
wald <- pop_wald_intervals(tibble::tibble(estimate = c(1000, 1000, Inf),
                                         matches = c(4L, 1L, 0L)))
stopifnot(wald$variance[1] == 250000, wald$se[1] == 500,
          abs(wald$lower[1] - 20) < 1e-10, abs(wald$upper[1] - 1980) < 1e-10,
          wald$lower[2] < 0, abs(wald$lower[2] + 960) < 1e-10,
          is.na(wald$lower[3]), is.na(wald$upper[3]),
          !wald$interval_available[3])
# Tied intervals keep replicate-ID order, and ranks belong to sorted intervals.
interval_fixture <- list(summary = tibble::tibble(
  sex = rep(c("Female", "Male"), each = 4), method = "Moments",
  estimate = rep(c(1000, 1000, 1000, Inf), 2),
  matches = rep(c(4L, 1L, 4L, 0L), 2), rep_id = rep(c(8L, 5L, 2L, 1L), 2),
  truth = 1000), settings = list(Nf = 1000, Nm = 1000))
interval_data <- pop_interval_data(interval_fixture)
for (sex in c("Female", "Male")) {
  finite <- interval_data[interval_data$sex == sex & interval_data$interval_available, ]
  stopifnot(identical(finite$rep_id, c(5L, 2L, 8L)),
            identical(finite$interval_rank, 1:3))
}

# Finite-only summaries exclude Inf and handle zero/one finite observation.
estimate_check <- pop_summarize_estimates(tibble::tibble(
  sex = c(rep("Female", 3), rep("Male", 3), rep("Female", 2), rep("Male", 2)),
  method = c(rep("Moments", 6), rep("Zero", 2), rep("Constant", 2)),
  truth = 1000, estimate = c(Inf, 500, 1000, Inf, Inf, 200, Inf, Inf, 100, 100),
  matches = c(0L, 4L, 4L, 0L, 0L, 1L, 0L, 0L, 4L, 4L)
))
stopifnot(estimate_check$finite_mean[1] == 750,
          estimate_check$finite_median[1] == 750,
          abs(estimate_check$finite_cv[1] - sqrt(2) / 3) < 1e-12,
          estimate_check$num_finite[1] == 2, estimate_check$num_unbounded[1] == 1,
          estimate_check$finite_mean[2] == 200, is.na(estimate_check$finite_cv[2]),
          is.na(estimate_check$finite_mean[3]), is.na(estimate_check$finite_median[3]),
          is.na(estimate_check$finite_cv[3]), estimate_check$num_unbounded[3] == 2,
          estimate_check$finite_cv[4] == 0,
          estimate_check$num_ci_included[1] == 1,
          estimate_check$finite_ci_inclusion[1] == 0.5,
          estimate_check$finite_ci_inclusion[2] == 0,
          estimate_check$num_ci_included[3] == 0,
          is.na(estimate_check$finite_ci_inclusion[3]))
# Coverage includes either endpoint, without excluding an interval's boundary.
boundaries <- pop_summarize_estimates(tibble::tibble(
  sex = "Female", method = "Moments", estimate = 1000, matches = 4L,
  truth = c(wald$lower[1], wald$upper[1])))
stopifnot(all(boundaries$num_ci_included == 1), all(boundaries$finite_ci_inclusion == 1))

# Direct likelihood evaluation checks the claimed analytic maximum.
m <- 7
S <- 12
K <- 20
Nhat <- S * K / m
loglik <- function(N) stats::dpois(m, lambda = S * K / N, log = TRUE)
stopifnot(loglik(Nhat) > loglik(Nhat * 0.95), loglik(Nhat) > loglik(Nhat * 1.05))

# Independent simulations agree with the theoretical means for unequal sexes.
means_trial <- pop_run_trial(pop_settings(20, 40, 60, 10, 10, 20), 1000, 2001)
means <- subset(means_trial$summary, method == "Moments")
stopifnot(abs(mean(means$matches[means$sex == "Female"]) - 10) < 0.4,
          abs(mean(means$matches[means$sex == "Male"]) - 5) < 0.4)

# Every plotting path must work, including no finite estimates and a census.
# Weighted histogram heights must represent individuals, not summary-table rows.
population_bars <- ggplot2::ggplot_build(pop_plot_reproduction(a))$data[[1]]
sample_bars <- ggplot2::ggplot_build(pop_plot_reproduction(a, TRUE))$data[[1]]
stopifnot(sum(population_bars$count) == 17 + 31,
          sum(sample_bars$count) == 8 + 11)
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file, width = 10, height = 6)
print(pop_plot_reproduction(a))
print(pop_plot_reproduction(a, TRUE))
print(pop_plot_pairs(sparse))
print(pop_plot_estimates(sparse))
all_zero <- sparse
all_zero$summary$estimate <- Inf
all_zero$summary$finite <- FALSE
print(pop_plot_estimates(all_zero))
print(pop_plot_intervals(sparse))
print(pop_plot_intervals(all_zero))
print(pop_plot_intervals(interval_fixture, rep_id = 2))
print(pop_plot_intervals(interval_fixture, rep_id = 1))
print(pop_plot_pairs(pop_run_trial(pop_settings(1, 3, 20, 1, 3, 20), 1, 1)))
grDevices::dev.off()
unlink(plot_file)

if (requireNamespace("shiny", quietly = TRUE)) {
  source("R/pop-app.R")
  # Budget boundaries count retained samples, not populations across replicates.
  boundary <- pop_settings(50000, 50000, 100000, 50000, 50000, 100000)
  stopifnot(.pop_check_budget(boundary, 100) == 2e7,
            .pop_check_budget(boundary, 100, stored_samples = 2e7) == 2e7,
            .pop_check_budget(pop_settings(4999999, 1, 5000000, 1, 1, 1), 500) == 1500)
  expect_error(.pop_check_budget(boundary, 101))
  expect_error(.pop_check_budget(boundary, 100, stored_samples = 2e7 + 1))
  expect_error(.pop_check_budget(pop_settings(5000000, 1, 5000000, 1, 1, 1), 1))
  expect_error(.pop_check_budget(boundary, 1, num_trials = 30))
  html <- as.character(pop_app_ui())
  stopifnot(grepl("parent-offspring pair fundamentals", html, fixed = TRUE),
            grepl("An increasingly interesting series of simple parent-offspring scenarios", html, fixed = TRUE),
            grepl("Start with the simplest possible population", html, fixed = TRUE))
  shiny::testServer(pop_app_server, {
    session$setInputs(Nf = 17, Nm = 31, Nk = 80, Sf = 8, Sm = 11, Sk = 40,
                      num_reps = 3, seed = 42, run_trial = 1)
    stopifnot(length(trials()) == 1)
    stopifnot(stored_individuals() == 3 * (8 + 11 + 40),
              is.null(trials()[[1]]$reps[[1]]$population))
    session$setInputs(trial_id = "1", rep_id = 2)
    original <- current_trial()
    stopifnot(current_rep()$seed == 43)
    # Input edits don't rewrite the saved trial.
    session$setInputs(Nf = 20, seed = 100)
    stopifnot(identical(current_trial(), original))
    session$setInputs(run_trial = 2)
    stopifnot(length(trials()) == 2, trials()[["2"]]$settings$Nf == 20,
              identical(trials()[["1"]], original))
    session$setInputs(trial_id = "2", rep_id = 3)
    stopifnot(current_rep()$seed == 102)
    # Invalid requests fail without appending a partial trial.
    session$setInputs(Sf = 99, run_trial = 3)
    stopifnot(length(trials()) == 2)
    session$setInputs(Sf = 8, num_reps = 501, run_trial = 4)
    stopifnot(length(trials()) == 2)
    session$setInputs(num_reps = 1, Nf = 10000000, run_trial = 5)
    stopifnot(length(trials()) == 2)

    # Regression: a very sparse male estimate distribution must not display the
    # true abundance as its finite median. Use an explicit summary fixture with
    # 98 zero-pair reps, one one-pair rep, and one two-pair rep for each method.
    # Sm * Sk = 1000 gives finite estimates 1000 and 500, whose median is 750.
    edge <- original
    edge$num_reps <- 100L
    edge$summary <- do.call(rbind, lapply(seq_len(100), function(r) {
      male_matches <- if (r <= 98) 0L else if (r == 99) 1L else 2L
      pairs <- tibble::tibble(sex = c(rep("Female", 10), rep("Male", male_matches)))
      result <- pop_estimate(pairs, Sf = 10, Sm = 10, Sk = 100)
      result$truth <- ifelse(result$sex == "Female", 100L, 50000L)
      result$rep_id <- r
      result$seed <- r
      result
    }))
    trials(list("1" = edge))
    session$setInputs(trial_id = "1", rep_id = 1)
    displayed <- output$estimation_table
    stopifnot(grepl('Male </td> <td> Moments </td> <td align="right"> 50000 </td> <td> 750 </td> <td> 750 </td> <td> 47.1% </td> <td> - </td> <td> 98 / 100', displayed, fixed = TRUE),
              grepl('Male </td> <td> Poisson pseudo-MLE </td> <td align="right"> 50000 </td> <td> 750 </td> <td> 750 </td> <td> 47.1% </td> <td> 0 / 2 (0.0%) </td> <td> 98 / 100', displayed, fixed = TRUE))
  })
  cat("Shiny session checks passed.\n")
} else {
  cat("Shiny not installed: server checks skipped.\n")
}
cat("Simulation, estimation, archive, and plotting checks passed.\n")
