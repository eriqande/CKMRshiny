# Run from repository root: Rscript tests/check-size-reproduction.R
source("R/pop-simulation.R")
source("R/pop-plots.R")
source("R/pop-app.R")
expect_error <- function(code) stopifnot(inherits(tryCatch({force(code); NULL}, error = identity), "error"))
expect_error(pop_size_parameters(t0 = 1))
expect_error(pop_size_parameters(max_age = 1))
expect_error(pop_size_parameters(rho = 1.1))
expect_error(pop_size_parameters(output_cv = -1))
expect_error(pop_settings(size_model = list(Female = pop_size_parameters())))

model <- list(Female = pop_size_parameters(), Male = pop_size_parameters(exponent = 1))
settings <- pop_settings(Nf = 1000, Nm = 1500, Nk = 4000, Sf = 400, Sm = 500, Sk = 1000,
                         size_model = model)
set.seed(32)
state <- .Random.seed
r <- pop_simulate_rep(settings, 456)
stopifnot(identical(state, .Random.seed), identical(r, pop_regenerate_rep(r)),
          is.null(r$population), nrow(r$display_adults) == 1000,
          nrow(r$samples$females) == 400)
full <- pop_regenerate_rep(r, retain_population = TRUE)
new_sample <- pop_regenerate_rep(r, Sf = 100, sampling_seed = 717)
stopifnot(identical(r$size_summary, new_sample$size_summary),
          identical(r$display_adults, new_sample$display_adults),
          identical(r$reproduction_summary, new_sample$reproduction_summary))
for (s in c("Female", "Male")) {
  adults <- full$population[full$population$sex == s, ]
  stopifnot(abs(mean(adults$reproductive_weight) - 1) < 1e-12,
            all(adults$age >= 1 & adults$age < 30), all(adults$length > 0))
  expected_length <- adults$individual_Linf * (1 - exp(-adults$individual_K * adults$age))
  stopifnot(max(abs(adults$length - expected_length)) < 1e-12)
}
# Population normalization removes absolute fecundity scale; zero exponent and
# zero residual CV give exact equal reproduction irrespective of length or age.
equal <- .pop_with_seed(99, pop_size_traits(500, pop_size_parameters(exponent = 0, output_cv = 0)))
stopifnot(all(equal$reproductive_weight == 1))
# Check draws against the truncated age CDF and intended log-trait correlation.
traits <- .pop_with_seed(399, pop_size_traits(50000, model$Female))
p <- model$Female
u <- (1 - exp(-p$Z * (traits$age - p$min_age))) / (1 - exp(-p$Z * (p$max_age - p$min_age)))
stopifnot(abs(mean(u) - .5) < .005, abs(mean(u < .25) - .25) < .01,
          abs(cor(log(traits$individual_Linf), log(traits$individual_K)) - p$rho) < .02)
relative_noise <- traits$reproductive_weight / traits$expected_reproductive_weight
stopifnot(abs(mean(relative_noise) - 1) < .01,
          abs(sd(relative_noise) / mean(relative_noise) - p$output_cv) < .01)

# Repeated matched parents contribute once per pair, not once per unique adult.
pop <- data.frame(sex = rep(c("Female", "Male"), each = 2), parent_id = rep(1:2, 2),
                  length = rep(c(10, 20), 2), reproductive_weight = rep(c(.5, 1.5), 2))
known <- list(size_summary = pop_summarize_size(pop),
              samples = list(females = pop[1:2, ], males = pop[3:4, ]),
              pairs = data.frame(sex = rep("Female", 4), parent_id = c(1, 2, 2, 2)))
for (variable in c("length", "reproductive_weight")) {
  d <- pop_pair_trait_data(known, variable)
  f <- d[d$sex == "Female", ]
  match_bins <- f[f$source == "Identified pairs" & f$count > 0, ]
  stopifnot(identical(as.numeric(match_bins$count), c(1, 3)), all(match_bins$total == 4))
  for (source in unique(f$source)) {
    x <- f[f$source == source, ]
    stopifnot(abs(sum(x$density * (x$upper - x$lower)) - 1) < 1e-12)
  }
  expected <- f[f$source == "Expected pair parents" & f$count > 0, ]
  stopifnot(max(abs(expected$count / expected$total - c(.25, .75))) < 1e-12)
  stopifnot(all(d$density[d$sex == "Male" & d$source == "Identified pairs"] == 0))
}
# Full-population summaries, not display subsamples, determine the expected curve.
for (variable in c("length", "reproductive_weight")) {
  d <- pop_pair_trait_data(r, variable)
  for (s in c("Female", "Male")) {
    a <- full$population[full$population$sex == s, ]
    x <- d[d$sex == s & d$source == "Expected pair parents", ]
    bins <- cut(a[[variable]], c(x$lower, tail(x$upper, 1)), include.lowest = TRUE, right = FALSE)
    direct <- tapply(a$reproductive_weight, bins, sum)
    direct[is.na(direct)] <- 0
    stopifnot(max(abs(x$count - direct)) < 1e-10)
  }
}
preview <- pop_size_preview(model)
stopifnot(identical(preview, pop_size_preview(model)), identical(state, .Random.seed))
# The right axes recover discrete counts exactly, independently for each sex
# and variable, including repeated parents. Zero matches have no count axis.
for (variable in c("length", "reproductive_weight")) {
  plot <- pop_plot_pair_traits(r, variable)
  data <- pop_pair_trait_data(r, variable)
  for (s in c("Female", "Male")) {
    bins <- data[data$sex == s & data$source == "Identified pairs", ]
    axis <- plot$panels[[s]]$scales$get_scales("y")$secondary.axis
    stopifnot(max(abs(axis$trans(bins$density) - bins$count)) < 1e-10,
              all(axis$breaks(c(0, 7.2)) %% 1 == 0))
  }
  no_pairs <- pop_plot_pair_traits(known, variable)
  stopifnot(inherits(no_pairs$panels$Male$scales$get_scales("y")$secondary.axis, "waiver"))
}
file <- tempfile(fileext = ".pdf")
grDevices::pdf(file, width = 10, height = 5)
for (variable in c("length", "reproductive_weight")) {
  print(pop_plot_pair_traits(r, variable))
  print(pop_plot_pair_traits(known, variable))
}
for (kind in c("growth", "output")) {
  print(pop_plot_size_relationship(r$display_adults, model, kind, r$size_summary))
  print(pop_plot_size_relationship(preview, model, kind)) # Exercises large-sample binning.
}
print(pop_plot_erro(r$size_summary))
print(pop_plot_erro(r$size_summary, "log10"))
grDevices::dev.off()
unlink(file)
# Exercise the stacked layout on a narrow device too.
file <- tempfile(fileext = ".pdf")
grDevices::pdf(file, width = 4, height = 9)
print(pop_plot_pair_traits(known))
grDevices::dev.off()
unlink(file)

# Real module integration: immutable saved settings, independent registries,
# validation, shared budgets, trait toggles, and numeric flattened CSV columns.
inputs <- list()
for (s in c("Female", "Male")) for (name in names(model[[s]]))
  inputs[[paste(tolower(s), name, sep = "_")]] <- model[[s]][[name]]
shiny::testServer(.pop_size_scenario_server, {
  handles <- session$returned
  do.call(session$setInputs, c(inputs, list(Nf = 100, Nm = 120, Nk = 400,
    Sf = 50, Sm = 60, Sk = 200, num_reps = 3, seed = 711, run_trial = 1)))
  stopifnot(length(handles$trials()) == 1, handles$stored_individuals() == 930)
  session$setInputs(trial_id = "1", rep_id = 1)
  original <- handles$current_trial()
  stopifnot(original$scenario == "Size-based reproduction", nrow(pop_size_trial_table(handles$trials())) == 1,
            pop_size_trial_table(handles$trials())$female_exponent == 3)
  session$setInputs(female_exponent = 0, pair_trait_toggle = 1, erro_y_toggle = 1)
  stopifnot(identical(original, handles$current_trial()))
  # Exercise renderers after the saved trial has been selected.
  stopifnot(!is.null(output$pair_traits_plot), !is.null(output$sample_growth),
            !is.null(output$population_output), !is.null(output$population_erro))
  session$setInputs(female_t0 = 4, run_trial = 2)
  stopifnot(length(handles$trials()) == 1)
})
html <- as.character(pop_app_ui())
ids <- regmatches(html, gregexpr(' id="[^"]+"', html))[[1]]
stopifnot(!anyDuplicated(ids), grepl('id="size-pair_trait_toggle"', html, fixed = TRUE),
          !grepl('id="size-use_weights"', html, fixed = TRUE))
cat("Size-based reproduction, pair enrichment, regeneration, and UI checks passed.\n")

# Three independent registries consume the same session budget.
shiny::testServer(pop_app_server, {
  prefixed <- stats::setNames(inputs, paste0("size-", names(inputs)))
  do.call(session$setInputs, c(prefixed, list(`size-Nf` = 100, `size-Nm` = 120,
    `size-Nk` = 400, `size-Sf` = 50, `size-Sm` = 60, `size-Sk` = 200,
    `size-num_reps` = 2, `size-seed` = 711, `size-run_trial` = 1)))
  stopifnot(length(size$trials()) == 1, length(equal$trials()) == 0,
            length(unequal$trials()) == 0, usage$samples_size == 620,
            usage$trials_size == 1)
  session$setInputs(Nf = 100, Nm = 120, Nk = 400, Sf = 50, Sm = 60, Sk = 200,
                    num_reps = 2, seed = 711, run_trial = 1)
  stopifnot(length(equal$trials()) == 1, usage$samples_equal == 620)
  usage$samples_size <- .pop_app_limits()$session_samples
  session$setInputs(run_trial = 2)
  stopifnot(length(equal$trials()) == 1)
})
cat("Three-tab shared budget and namespace checks passed.\n")
