# Run from repository root: Rscript tests/check-selectivity.R
source("R/pop-simulation.R")
source("R/pop-plots.R")
source("R/pop-app.R")
expect_error <- function(code) stopifnot(inherits(tryCatch({force(code); NULL}, error = identity), "error"))
model <- list(Female = pop_size_parameters(), Male = pop_size_parameters(exponent = 1))
s <- pop_settings(Nf = 400, Nm = 500, Nk = 1600, Sf = 40, Sm = 50, Sk = 500, size_model = model)
set.seed(772)
state <- .Random.seed
baseline <- pop_run_trial(s, 4, 451)
original <- serialize(baseline, NULL)
full_windows <- .pop_trial_length_bounds(baseline)
full_range <- pop_resample_size_trial(baseline, full_windows, 2)
stopifnot(identical(state, .Random.seed), identical(original, serialize(baseline, NULL)),
          full_range$source_trial_id == 1L, nrow(full_range$summary) == nrow(baseline$summary))
for (i in 1:4) stopifnot(identical(baseline$reps[[i]]$samples, full_range$reps[[i]]$samples),
                        identical(baseline$reps[[i]]$pairs, full_range$reps[[i]]$pairs))
w <- list(Female = c(55, full_windows$Female[2]), Male = full_windows$Male)
selected <- pop_resample_size_trial(baseline, w, 3)
other <- pop_resample_size_trial(baseline, list(Female = c(65, full_windows$Female[2]), Male = full_windows$Male), 4)
for (i in 1:4) {
  r <- selected$reps[[i]]
  stopifnot(identical(r, pop_regenerate_rep(r)),
    identical(r$samples$kids, baseline$reps[[i]]$samples$kids),
    identical(r$samples$males, baseline$reps[[i]]$samples$males),
    identical(r$samples$males, other$reps[[i]]$samples$males),
    identical(r$samples$kids, other$reps[[i]]$samples$kids),
    identical(r$reproduction_summary, baseline$reps[[i]]$reproduction_summary),
    identical(r$display_adults, baseline$reps[[i]]$display_adults),
    identical(r$size_summary, baseline$reps[[i]]$size_summary),
    all(r$samples$females$length >= 55), !anyDuplicated(r$samples$females$parent_id),
    nrow(r$samples$females) == s$Sf, is.null(r$population), is.null(r$pedigree),
    identical(r$estimates, pop_estimate(pop_identify_pairs(r$samples), s$Sf, s$Sm, s$Sk)))
  # Widening a selected source back to full range restores the original baseline.
  restored <- pop_resample_size_rep(r, full_windows)
  stopifnot(identical(restored$samples, baseline$reps[[i]]$samples))
  population <- pop_regenerate_rep(baseline$reps[[i]], retain_population = TRUE)$population
  # Compare expected-pair curves with direct sums over eligible population adults.
  for (sex in c("Female", "Male")) for (variable in c("length", "reproductive_weight")) {
    a <- population[population$sex == sex & population$length >= w[[sex]][1] & population$length <= w[[sex]][2], ]
    stopifnot(nrow(a) == r$eligible_counts[[sex]])
    data <- pop_pair_trait_data(r, variable)
    bins <- data[data$sex == sex & data$source == "Expected pair parents", ]
    b <- cut(a[[variable]], c(bins$lower, tail(bins$upper, 1)), include.lowest = TRUE, right = FALSE)
    direct <- tapply(a$reproductive_weight, b, sum)
    direct[is.na(direct)] <- 0
    stopifnot(max(abs(bins$count - direct)) < 1e-10)
  }
}
expect_error(pop_length_windows(list(Female = c(2,1), Male = c(0,100))))
expect_error(pop_length_windows(list(Female = c(NA,100), Male = c(0,100))))
expect_error(pop_resample_size_trial(baseline, list(Female = c(1e8,1e9), Male = full_windows$Male), 5))
# Failure after a successful replicate must still return no partial trial.
cutoffs <- vapply(baseline$reps, function(r) {
  a <- pop_regenerate_rep(r, retain_population=TRUE)$population
  sort(a$length[a$sex == "Female"], decreasing=TRUE)[s$Sf]
}, numeric(1))
order <- c(which.max(cutoffs), setdiff(seq_along(cutoffs), which.max(cutoffs)))
reordered <- baseline
reordered$reps <- baseline$reps[order]
completed <- 0L
failure <- tryCatch(pop_resample_size_trial(reordered,
  list(Female=c(max(cutoffs), Inf), Male=c(0, Inf)), 5,
  progress=function(done,total) completed <<- done), error=identity)
stopifnot(inherits(failure,"error"), completed == 1, grepl("Replicate 2",conditionMessage(failure)))
stopifnot(identical(original, serialize(baseline, NULL)), identical(state, .Random.seed))
register <- pop_selectivity_trial_table(list(baseline, selected))
stopifnot(is.na(register$source_trial_id[1]), register$source_trial_id[2] == 1,
          register$female_length_min[2] == 55, register$sampling[1] == "Uniform baseline")
# Exercise gray-region previews and actual-sample overlays without changing RNG.
for (kind in c("growth", "output")) for (show in c(FALSE, TRUE)) {
  p <- pop_plot_selectivity(selected$reps[[1]], w, kind, show)
  stopifnot(inherits(p, "ggplot"), length(ggplot2::ggplot_build(p)$data) >= 5)
}
# An unchanged estimator must show the bias induced by sampling large females.
big_s <- pop_settings(Nf = 1000, Nm = 1000, Nk = 2000, Sf = 100, Sm = 100, Sk = 500, size_model = model)
comparison <- pop_run_trial(big_s, 60, 1601)
large <- pop_resample_size_trial(comparison, list(Female = c(75, Inf), Male = c(0, Inf)), 2)
f <- large$summary$sex == "Female" & large$summary$method == "Moments"
b <- comparison$summary$sex == "Female" & comparison$summary$method == "Moments"
stopifnot(mean(large$summary$matches[f]) > 1.5 * mean(comparison$summary$matches[b]),
          mean(large$summary$estimate[f]) < .8 * big_s$Nf)

inputs <- list()
for (sex in names(model)) for (name in names(model[[sex]]))
  inputs[[paste(tolower(sex), name, sep = "_")]] <- model[[sex]][[name]]
module_usage <- .pop_session_usage()
shiny::testServer(.pop_size_scenario_server, args = list(selective = TRUE, shared_usage = module_usage), {
  handles <- session$returned
  do.call(session$setInputs, c(inputs, list(Nf=120, Nm=150, Nk=400, Sf=20, Sm=30, Sk=120,
    num_reps=3, seed=711, run_trial=1)))
  session$setInputs(trial_id="1", rep_id=1)
  baseline <- handles$current_trial()
  bounds <- handles$selection$bounds()
  session$setInputs(female_length_window=bounds$Female, male_length_window=bounds$Male)
  stopifnot(handles$selection$applied(), length(handles$trials()) == 1)
  session$setInputs(female_length_window=c(45,bounds$Female[2]))
  stopifnot(!handles$selection$applied(), identical(baseline, handles$current_trial()),
            grepl("Preview only", output$selection_status, fixed=TRUE))
  session$setInputs(apply_selection=1)
  stopifnot(length(handles$trials()) == 2, handles$stored_individuals() == 1020)
  session$setInputs(trial_id="2")
  stopifnot(handles$selection$applied(), handles$current_trial()$source_trial_id == 1,
            identical(baseline, handles$trials()[["1"]]),
            !is.null(output$sample_output), !is.null(output$sample_growth), !is.null(output$pair_traits_plot),
            grepl("Eligible individuals", output$sampling_table, fixed=TRUE))
  # Even invalid edits to next-trial controls cannot alter resampling inputs.
  session$setInputs(female_t0=99, female_length_window=c(1e8,1e9), apply_selection=2)
  stopifnot(length(handles$trials()) == 2)
  shared_usage$samples_equal <- .pop_app_limits()$session_samples
  session$setInputs(female_length_window=bounds$Female, apply_selection=3)
  stopifnot(length(handles$trials()) == 2)
})
shiny::testServer(pop_app_server, {
  prefixed <- stats::setNames(inputs, paste0("breaking-", names(inputs)))
  do.call(session$setInputs, c(prefixed, list(`breaking-Nf`=100, `breaking-Nm`=120,
    `breaking-Nk`=400, `breaking-Sf`=20, `breaking-Sm`=30, `breaking-Sk`=120,
    `breaking-num_reps`=2, `breaking-seed`=711, `breaking-run_trial`=1)))
  stopifnot(length(breaking$trials()) == 1, length(size$trials()) == 0,
            length(equal$trials()) == 0, usage$samples_breaking == 340,
            .pop_usage_total(usage,"trials") == 1)
})
html <- as.character(pop_app_ui())
ids <- regmatches(html, gregexpr(' id="[^"]+"', html))[[1]]
stopifnot(!anyDuplicated(ids), grepl('id="breaking-apply_selection"', html, fixed=TRUE),
          grepl('Breaking Stuff!',html,fixed=TRUE))
cat("Length selection, fixed kids, immutable trials, diagnostic expectations, bias, and four-tab budgets passed.\n")
