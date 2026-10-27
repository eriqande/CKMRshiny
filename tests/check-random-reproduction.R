# Run from the repository root: Rscript tests/check-random-reproduction.R
source("R/pop-simulation.R")
source("R/pop-plots.R")
expect_error <- function(code) {
  stopifnot(tryCatch({ force(code); FALSE }, error = function(e) TRUE))
}
expect_error(pop_settings(use_weights = NA))
expect_error(pop_settings(use_weights = 1))
expect_error(pop_settings(female_shape = 0))
expect_error(pop_settings(male_shape = Inf))
expect_error(pop_settings(female_shape = 0.001))
expect_error(pop_settings(male_shape = 1001))

settings <- pop_settings(40, 70, 300, 10, 20, 60,
                         use_weights = TRUE, female_shape = 0.2, male_shape = 2)
set.seed(180)
old_seed <- .Random.seed
replicate <- pop_simulate_rep(settings, 99)
stopifnot(identical(old_seed, .Random.seed),
          identical(replicate, pop_regenerate_rep(replicate)),
          is.null(replicate$population), is.null(replicate$pedigree),
          nrow(replicate$sampling_weight_summary) == 2,
          identical(replicate$sampling_weight_summary$n_adults, c(40L, 70L)))
full <- pop_regenerate_rep(replicate, retain_population = TRUE)
resampled <- pop_regenerate_rep(replicate, Sf = 20, sampling_seed = 100,
                               retain_population = TRUE)
stopifnot(all(is.finite(full$population$reproductive_weight)),
          all(full$population$reproductive_weight >= 0),
          length(unique(full$population$reproductive_weight)) > 2,
          identical(full$population, resampled$population),
          identical(full$pedigree, resampled$pedigree),
          identical(full$reproduction_summary, resampled$reproduction_summary),
          identical(full$sampling_weight_summary, resampled$sampling_weight_summary))
effective <- pop_effective_numbers(replicate)
stopifnot(identical(effective, pop_effective_numbers(resampled)))
for (i in seq_len(nrow(effective))) {
  ids <- if (effective$sex[i] == "Female") full$pedigree$mother_id else full$pedigree$father_id
  counts <- tabulate(ids, nbins = effective$census_n[i])
  h <- sum(choose(counts, 2)) / choose(nrow(full$pedigree), 2)
  stopifnot(abs(h - effective$same_parent_probability[i]) < 1e-12)
}

# Changing reproductive sampling weights cannot change which IDs are surveyed:
# the survey has its own seed and never uses those sampling weights as probs.
equal_settings <- settings
equal_settings$use_weights <- FALSE
equal <- pop_simulate_rep(equal_settings, 99)
stopifnot(identical(replicate$samples$females$parent_id, equal$samples$females$parent_id),
          identical(replicate$samples$males$parent_id, equal$samples$males$parent_id),
          identical(replicate$samples$kids$kid_id, equal$samples$kids$kid_id),
          all(equal$samples$females$reproductive_weight == 1),
          all(equal$samples$males$reproductive_weight == 1))
baseline <- pop_simulate_rep(pop_settings(40, 70, 300, 10, 20, 60), 99)
stopifnot(identical(baseline$samples, equal$samples), identical(baseline$pairs, equal$pairs))

# Within-sex multiplication of reproductive sampling weights changes no
# parental probabilities; compare identical draws with a power-of-two factor.
scaled <- full$population
scaled$reproductive_weight[scaled$sex == "Female"] <-
  scaled$reproductive_weight[scaled$sex == "Female"] * 8
scaled$reproductive_weight[scaled$sex == "Male"] <-
  scaled$reproductive_weight[scaled$sex == "Male"] * 0.25
stopifnot(identical(.pop_with_seed(400, pop_reproduce(full$population, 500)),
                    .pop_with_seed(400, pop_reproduce(scaled, 500))))

# Large draws check both independent sex-specific gamma parameterizations.
population <- .pop_with_seed(100, pop_generate_population(10000, 12000, TRUE, 0.5, 3))
summary <- pop_summarize_sampling_weights(population)
stopifnot(abs(summary$mean[1] - 0.5) < 0.025,
          abs(summary$mean[2] - 3) < 0.1,
          abs(summary$cv[1] - sqrt(2)) < 0.12,
          abs(summary$cv[2] - 1 / sqrt(3)) < 0.04)

# Uniform surveys preserve the expected match count despite strong reproductive
# inequality. The beta-binomial marginal yields variance 50*.2*.8*(50+5)/(1+5).
trial <- pop_run_trial(pop_settings(50, 80, 100, 10, 16, 50,
                                    use_weights = TRUE, female_shape = 0.1, male_shape = 1),
                       num_reps = 1000, seed = 4001)
female <- subset(trial$summary, method == "Moments" & sex == "Female")
male <- subset(trial$summary, method == "Moments" & sex == "Male")
expected_variance <- 50 * 0.2 * 0.8 * (50 + 5) / (1 + 5)
stopifnot(abs(mean(female$matches) - 10) < 1.5,
          abs(mean(male$matches) - 10) < 1,
          abs(stats::var(female$matches) / expected_variance - 1) < 0.25)
register <- pop_trial_table(list(trial), include_weights = TRUE)
stopifnot(register$use_sampling_weights, register$female_gamma_shape == 0.1,
          register$male_gamma_shape == 1,
          !"female_gamma_shape" %in% names(pop_trial_table(list(trial))))

# Theoretical density previews neither use RNG nor generate infinite curve data.
old_seed <- .Random.seed
for (shape in c(0.05, 0.2, 1, 2, 1000)) {
  d <- pop_gamma_density(shape)
  stopifnot(all(is.finite(d$density)), all(d$density >= 0),
            abs(max(d$sampling_weight) - qgamma(0.995, shape)) < 1e-8)
}
stopifnot(identical(old_seed, .Random.seed))
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file, width = 8, height = 5)
for (shape in c(0.05, 1, 1000)) {
  print(pop_plot_gamma(shape, "Female", "linear"))
  print(pop_plot_gamma(shape, "Male", "log10"))
}
print(pop_plot_reproduction(replicate))
print(pop_plot_reproduction(replicate, TRUE, "log10"))
print(pop_plot_pairs(trial))
print(pop_plot_intervals(trial))
grDevices::dev.off()
unlink(plot_file)

if (requireNamespace("shiny", quietly = TRUE)) {
  source("R/pop-app.R")
  html <- as.character(pop_app_ui())
  stopifnot(grepl('id="unequal-use_weights"', html, fixed = TRUE),
            !grepl('id="use_weights"', html, fixed = TRUE),
            grepl('id="unequal-female_shape"', html, fixed = TRUE),
            grepl("00", html, fixed = TRUE))
  ids <- regmatches(html, gregexpr(' id="[^"]+"', html))[[1]]
  stopifnot(!anyDuplicated(ids))

  shiny::testServer(.pop_unequal_server, {
    session$setInputs(Nf = 40, Nm = 70, Nk = 300, Sf = 10, Sm = 20, Sk = 60,
                      num_reps = 3, seed = 99, use_weights = TRUE,
                      female_shape = 0.2, male_shape = 2, run_trial = 1)
    handles <- session$returned
    stopifnot(length(handles$trials()) == 1)
    session$setInputs(trial_id = "1", rep_id = 2)
    original <- handles$current_trial()
    stopifnot(grepl("Female </td> <td> 40 </td>", output$effective_numbers_table, fixed = TRUE),
              grepl("Male </td> <td> 70 </td>", output$effective_numbers_table, fixed = TRUE))
    stopifnot(original$settings$use_weights, original$settings$female_shape == 0.2,
              grepl(session$ns("rep_id"), output$rep_control$html, fixed = TRUE))
    session$setInputs(female_shape = 0.05, male_shape = 4, gamma_y_toggle = 1)
    stopifnot(identical(handles$current_trial(), original),
              grepl("447.2%", output$female_gamma_summary, fixed = TRUE),
              grepl("50.0%", output$male_gamma_summary, fixed = TRUE))
    # Equal reproduction within tab two keeps the tab-two register and controls.
    session$setInputs(use_weights = FALSE, seed = 300, run_trial = 2)
    stopifnot(length(handles$trials()) == 2,
              !handles$trials()[[2]]$settings$use_weights,
              all(handles$trials()[[2]]$reps[[1]]$samples$females$reproductive_weight == 1),
              identical(handles$trials()[[1]], original))
    session$setInputs(female_shape = -1, run_trial = 3)
    stopifnot(length(handles$trials()) == 2)
  })

  # One real root session runs both tabs with distinct namespaces, saved results,
  # selectors, and shared sample/trial budgets.
  shiny::testServer(pop_app_server, {
    session$setInputs(Nf = 40, Nm = 70, Nk = 300, Sf = 10, Sm = 20, Sk = 60,
                      num_reps = 2, seed = 99, run_trial = 1)
    session$setInputs(`unequal-Nf` = 40, `unequal-Nm` = 70, `unequal-Nk` = 300,
                      `unequal-Sf` = 10, `unequal-Sm` = 20, `unequal-Sk` = 60,
                      `unequal-num_reps` = 3, `unequal-seed` = 99,
                      `unequal-use_weights` = TRUE, `unequal-female_shape` = 0.2,
                      `unequal-male_shape` = 2, `unequal-run_trial` = 1)
    stopifnot(length(equal$trials()) == 1, length(unequal$trials()) == 1,
              !equal$trials()[[1]]$settings$use_weights,
              unequal$trials()[[1]]$settings$use_weights,
              usage$samples_equal == 180, usage$samples_unequal == 270)
    first <- equal$trials()[[1]]
    session$setInputs(`unequal-female_shape` = 0.5, `unequal-run_trial` = 2)
    stopifnot(length(unequal$trials()) == 2, identical(equal$trials()[[1]], first))
    usage$samples_unequal <- .pop_app_limits()$session_samples
    session$setInputs(run_trial = 2)
    stopifnot(length(equal$trials()) == 1)
  })
  cat("Two-tab UI/session checks passed.\n")
}
cat("Random reproductive sampling-weight checks passed.\n")
