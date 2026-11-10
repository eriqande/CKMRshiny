# Run from project root: Rscript tests/check-hsp.R
for (file in c("pop-simulation.R", "pop-plots.R", "pop-app.R", "hsp-simulation.R",
               "hsp-plots.R", "hsp-app.R", "ckmr-app.R")) source(file.path("R", file))
expect_error <- function(code) stopifnot(inherits(tryCatch({force(code); NULL}, error = identity), "error"))
expect_error(hsp_settings(Sk1 = 1))
expect_error(hsp_settings(scenario = "cross_lucky", Sk2 = 3000))

# Direct enumeration gives an independent check of both match definitions and
# full-sibling bookkeeping; the production function never enumerates pairs.
first <- tibble::tibble(kid_id = 1:5, mother_id = c(1L,1L,1L,2L,3L), father_id = c(1L,1L,2L,2L,1L))
second <- tibble::tibble(kid_id = 1:4, mother_id = c(1L,2L,2L,3L), father_id = c(1L,1L,2L,2L))
direct <- function(a, b = NULL) {
  grid <- if (is.null(b)) t(combn(nrow(a), 2)) else as.matrix(expand.grid(seq_len(nrow(a)), seq_len(nrow(b))))
  if (is.null(b)) b <- a
  maternal <- a$mother_id[grid[,1]] == b$mother_id[grid[,2]]
  paternal <- a$father_id[grid[,1]] == b$father_id[grid[,2]]
  c(sum(maternal), sum(paternal), sum(maternal & paternal), nrow(grid))
}
within <- hsp_identify_matches(first, 4, 3)
across <- hsp_identify_matches(first, 4, 3, second)
stopifnot(identical(as.numeric(c(within$matches, within$full_siblings, within$comparisons)), as.numeric(direct(first))),
          identical(as.numeric(c(across$matches, across$full_siblings, across$comparisons)), as.numeric(direct(first,second))),
          sum(within$families$k1[within$families$sex=="Female"]) == nrow(first))
# Double arithmetic, not overflowing R integers, for large family counts.
big <- tibble::tibble(mother_id=rep(1L,70000), father_id=rep(1L,70000))
stopifnot(all(hsp_identify_matches(big,1,1)$matches == choose(70000,2)))

for (scenario in names(.hsp_titles())) {
  s <- hsp_settings(Nf=40, Nm=60, Nk1=300, Nk2=400, Sk1=70, Sk2=90,
                    scenario=scenario, female_shape=.3, male_shape=2)
  set.seed(300); original_rng <- .Random.seed
  r <- hsp_simulate_rep(s, 17)
  stopifnot(identical(.Random.seed, original_rng), identical(r, hsp_regenerate_rep(r)),
            is.null(r$population_1), is.null(r$pedigree_1), nrow(r$samples$cohort_1)==70)
  full <- hsp_regenerate_rep(r, TRUE)
  for (key in names(r)) stopifnot(identical(r[[key]], full[[key]]))
  if (s$cross) stopifnot(nrow(r$samples$cohort_2)==90,
    identical(full$population_1$parent_id, full$population_2$parent_id), r$comparisons == 70*90)
  else stopifnot(is.null(r$samples$cohort_2), r$comparisons==choose(70,2))
  if (scenario == "cross_lucky") stopifnot(!identical(full$population_1$reproductive_weight, full$population_2$reproductive_weight))
  if (scenario == "size") {
    stopifnot(all(full$population_2$length >= full$population_1$length),
      all(full$population_2$length <= full$population_2$individual_Linf),
      identical(full$population_2$age, full$population_1$age + 1),
      identical(full$population_1$individual_K, full$population_2$individual_K),
      s$size_model$Female$output_cv==0, s$size_model$Male$output_cv==0)
    for (sex in c("Female", "Male")) {
      bins <- r$size_summary[r$size_summary$sex==sex & r$size_summary$variable=="length", ]
      stopifnot(sum(bins$match_sum)==r$matches[[sex]], sum(bins$sample_sum)==160)
      a <- full$population_1[full$population_1$sex==sex, ]
      b <- full$population_2[full$population_2$sex==sex, ]
      stopifnot(abs(sum(bins$sibling_weight_sum) - sum(a$reproductive_weight*b$reproductive_weight)) < 1e-10,
                abs(sum(bins$sample_weight_sum) - sum(70*a$reproductive_weight+90*b$reproductive_weight)) < 1e-9)
    }
    for (episode in c(FALSE,TRUE)) {
      view <- .hsp_episode_traits(r,episode)
      ggplot2::ggplot_build(pop_plot_size_relationship(view$adults,view$model,"growth",view$summary))
      ggplot2::ggplot_build(pop_plot_size_relationship(view$adults,view$model,"output",view$summary))
    }
    ggplot2::ggplot_build(hsp_plot_persistence(r,lengths=TRUE))
  }
  # Census samples reproduce finite-population effective numbers exactly.
  census_settings <- s; census_settings$Sk1 <- s$Nk1
  if (s$cross) census_settings$Sk2 <- s$Nk2
  census <- hsp_simulate_rep(census_settings,17)
  stopifnot(max(abs(census$effective$same_parent_probability - census$matches/census$comparisons)) < 1e-12)
  trial <- hsp_run_trial(s,num_reps=3,seed=17)
  stopifnot(nrow(trial$summary)==12, nrow(hsp_trial_table(list(trial)))==1,
            identical(trial$reps[[1]],r))
  for (plot in list(hsp_plot_reproduction(r),hsp_plot_reproduction(r,TRUE,"log10"),
                   hsp_plot_pairs(trial),pop_plot_estimates(trial),pop_plot_intervals(trial))) ggplot2::ggplot_build(plot)
  if (s$cross) ggplot2::ggplot_build(hsp_plot_persistence(r))
  if (scenario=="size") for (matched in c(TRUE,FALSE)) for (variable in c("length","reproductive_weight"))
    ggplot2::ggplot_build(hsp_plot_traits(r,matched,variable))
}

# Independent seasonal Dirichlet weights have average sharing probability 1/N;
# within-season weights have (alpha+1)/(N*alpha+1). Check the simulation's own
# full-population propensities rather than a noisy inverse estimator.
simulate_prob <- function(scenario) {
  trial <- hsp_run_trial(hsp_settings(Nf=200,Nm=300,Nk1=300,Nk2=300,Sk1=60,Sk2=60,
    scenario=scenario,female_shape=.2,male_shape=1),num_reps=250,seed=101)
  rowMeans(vapply(trial$reps,function(r) r$effective$weight_probability,numeric(2)))
}
equal_p <- simulate_prob("equal"); lucky_p <- simulate_prob("lucky"); cross_p <- simulate_prob("cross_lucky")
stopifnot(max(abs(equal_p * c(200,300)-1))<1e-12,
          max(abs(lucky_p / c(1.2/41,2/301)-1))<.08,
          max(abs(cross_p * c(200,300)-1))<.08,
          lucky_p[1] > 4 * cross_p[1])

# All-unbounded intervals and a one-parent census remain plottable.
zero <- hsp_run_trial(hsp_settings(Nf=100000,Nm=100000,Nk1=2,Sk1=2),num_reps=2)
stopifnot(all(!zero$summary$finite))
invisible(ggplot2::ggplot_build(pop_plot_estimates(zero)))
invisible(ggplot2::ggplot_build(pop_plot_intervals(zero)))
one <- hsp_simulate_rep(hsp_settings(Nf=1,Nm=1,Nk1=5,Sk1=5),1)
stopifnot(all(one$matches==10), one$full_siblings==10, all(one$estimates$estimate==1))
stopifnot(all(hsp_estimate(c(0,0),10)$estimate==Inf))
stopifnot(!anyDuplicated(hsp_simulate_rep(hsp_settings(),.Machine$integer.max)$seeds))

# Zero annual growth-rate variation is exactly a one-year von B update on each
# individual's curve; large variation still cannot shrink or overshoot Linf.
model <- list(Female=pop_size_parameters(output_cv=0),Male=pop_size_parameters(output_cv=0))
adults <- .pop_with_seed(19,pop_generate_size_population(500,400,model))
deterministic <- hsp_grow_adults(adults,model,0)
theory <- adults$individual_Linf * (-expm1(-adults$individual_K * (adults$age + 1)))
stopifnot(max(abs(deterministic$length - theory))<1e-12)
variable <- .pop_with_seed(92,hsp_grow_adults(adults,model,3))
stopifnot(all(variable$length >= adults$length),all(variable$length <= adults$individual_Linf),
          any(abs(variable$length - deterministic$length)>1e-6))
fixed <- hsp_simulate_rep(hsp_settings(scenario="size",growth=FALSE),20,TRUE)
stopifnot(identical(fixed$population_1,fixed$population_2))
# Legacy schema-one archives used no growth. Their sampled pedigrees replay.
legacy <- fixed; legacy$schema_version <- 1L
legacy$settings$growth <- legacy$settings$growth_cv <- NULL
replayed <- hsp_regenerate_rep(legacy)
stopifnot(identical(fixed$samples,replayed$samples))
expect_error(hsp_settings(scenario="size",growth_cv=-.1))

flat <- hsp_simulate_rep(hsp_settings(Nf=30,Nm=50,scenario="size",
  size_model=list(Female=pop_size_parameters(exponent=0,output_cv=2),
                  Male=pop_size_parameters(exponent=0,output_cv=2))),31)
stopifnot(all(flat$display_adults$weight_1==1),all(flat$display_adults$weight_2==1),
          max(abs(flat$effective$weight_probability*c(30,50)-1))<1e-12)

if (requireNamespace("shiny", quietly=TRUE)) {
  # Boundary validation must not allocate the requested populations.
  shiny::isolate({
    usage <- .pop_session_usage()
    boundary <- hsp_settings(Nf=4e6,Nm=4e6,Nk1=1e6,Nk2=1e6,Sk1=1e6,Sk2=1e6,scenario="cross_lucky")
    stopifnot(.hsp_check_budget(boundary,10,usage)==2e7)
    expect_error(.hsp_check_budget(boundary,11,usage))
    boundary$Nk2 <- boundary$Nk2+1
    expect_error(.hsp_check_budget(boundary,1,usage))
  })
  shiny::testServer(ckmr_app_server, {
    stopifnot(inherits(tryCatch(hsp$equal$current_trial(),error=identity),"shiny.silent.error"))
    session$setInputs(Nf=40,Nm=60,Nk=300,Sf=10,Sm=10,Sk=60,num_reps=2,seed=10,run_trial=1)
    session$setInputs(`hs_equal-Nf`=40,`hs_equal-Nm`=60,`hs_equal-Nk1`=300,`hs_equal-Sk1`=60,
      `hs_equal-num_reps`=2,`hs_equal-seed`=10,`hs_equal-run_trial`=1)
    stopifnot(length(pop$equal$trials())==1,length(hsp$equal$trials())==1,
      .pop_usage_total(usage,"trials")==2, .pop_usage_total(usage,"samples")==280)
    saved_pop <- pop$equal$trials(); saved_hsp <- hsp$equal$trials()
    session$setInputs(open_pop=1,app_view="pop")
    session$setInputs(open_hsp=1,app_view="hsp")
    session$setInputs(choose_app=1,app_view="home")
    stopifnot(identical(pop$equal$trials(),saved_pop),identical(hsp$equal$trials(),saved_hsp))
    session$setInputs(`hs_lucky-Nf`=40,`hs_lucky-Nm`=60,`hs_lucky-Nk1`=300,`hs_lucky-Sk1`=60,
      `hs_lucky-use_weights`=TRUE,`hs_lucky-female_shape`=.2,`hs_lucky-male_shape`=1,
      `hs_lucky-num_reps`=2,`hs_lucky-seed`=10,`hs_lucky-run_trial`=1)
    session$setInputs(`hs_lucky-female_shape`=8,`hs_equal-population_y_toggle`=1)
    stopifnot(hsp$lucky$trials()[[1]]$settings$female_shape==.2,
      identical(hsp$equal$trials(),saved_hsp))
    usage$trials_hs_size <- 30
    session$setInputs(`hs_equal-run_trial`=2)
    stopifnot(length(hsp$equal$trials())==1)
  })
  # Every input/output ID is unique across the chooser and eight scenario UIs.
  html <- as.character(ckmr_app_ui())
  ids <- regmatches(html, gregexpr(' id="[^"]+"', html))[[1]]
  stopifnot(!anyDuplicated(ids))
}
cat("Sibling counting, independent luck, persistence, regeneration, plotting, switching, and shared budgets passed.\n")
