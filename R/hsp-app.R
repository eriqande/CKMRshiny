#' Titles for the four sibling teaching scenarios
#' @return Named character vector in teaching order.
#' @keywords internal
.hsp_titles <- function() c(equal = "Equal reproduction", lucky = "Lucky litters",
                            cross_lucky = "Cross-cohort siblings", size = "Persistent size effects")

#' Build the sibling application's four scenario tabs
#' @return Tag list suitable for embedding beside the POP app.
#' @export
hsp_app_ui <- function() {
  shiny::tagList(
    shiny::div(class = "hero", shiny::div(class = "eyebrow", "CKMR \u00b7 sibling pair fundamentals"),
      shiny::h1("From lucky litters to cross-cohort siblings"),
      shiny::p("Sample offspring, find shared parents, and explore which reproductive differences persist.")),
    do.call(shiny::tabsetPanel, c(list(id = "hsp_scenario"), lapply(names(.hsp_titles()), function(key)
      shiny::tabPanel(.hsp_titles()[[key]], .hsp_scenario_ui(key)))))
  )
}

#' Describe a sibling scenario and its learning goals
#' @param scenario A sibling scenario key.
#' @return Introductory section with the shared-parent counting convention.
#' @keywords internal
.hsp_intro <- function(scenario) {
  text <- switch(scenario,
    equal = "Start with one cohort produced by random mating. Every female has the same expected contribution, as does every male. Sample offspring uniformly and identify pairs sharing a mother or father. No adult survey sample is needed.",
    lucky = "Some parents get lucky: their offspring develop in a nutrient-rich pancake of water and contribute more strongly to the surviving cohort. Gamma reproductive sampling weights represent this seasonal advantage. Compare siblings within that single cohort and see what happens to the naive abundance estimator.",
    cross_lucky = "The same adults reproduce in two episodes and all survive between them. Each adult gets a fresh, independent gamma reproductive sampling weight in episode 2. Sample both cohorts uniformly and compare only offspring from different episodes. Does an advantage that lasts only one season still concentrate sibling matches across seasons?",
    size = "Bigger adults contribute more offspring in both episodes. Growth and mortality settings generate the initial age and length distribution. By default adults then grow for one year, with variation but no shrinking, before the second episode. All survive, and fecundity is recalculated from their new lengths, with no additional lucky-litter effect. Cross-cohort comparisons remove independent seasonal luck on average, but do they remove a persistent size advantage?")
  shiny::div(class = "section", shiny::h3(.hsp_titles()[[scenario]]), shiny::p(text),
    shiny::p("Mothers and fathers are chosen independently for every offspring. A maternal sibling pair contributes one female match; a paternal sibling pair contributes one male match. Full siblings contribute one match to each category. Kin identification is perfect; parental sex and identity are known in this simulation."),
    shiny::p(class = "hint", "We count shared-parent matches rather than strictly half-sibling pairs. A mating system in which one male fathers an entire female's brood would couple the maternal and paternal relationships; it would need separate treatment."))
}

#' Build one namespaced sibling scenario using the POP section layout
#' @param scenario Sibling scenario key.
#' @return Complete scenario UI with settings, sections 00-04, and saved trials.
#' @keywords internal
.hsp_scenario_ui <- function(scenario) {
  ns <- shiny::NS(paste0("hs_", scenario))
  cross <- scenario %in% c("cross_lucky", "size")
  lucky <- scenario %in% c("lucky", "cross_lucky")
  limits <- .pop_app_limits()
  shiny::tagList(.hsp_intro(scenario), shiny::sidebarLayout(
    shiny::sidebarPanel(width = 3, shiny::h3("Prepare a trial"),
      if (lucky) shiny::checkboxInput(ns("use_weights"), "Use random reproductive sampling weights", TRUE),
      shiny::p(class = "hint", "Controls apply to the next trial. Each replicate has its own population and offspring."),
      shiny::numericInput(ns("Nf"), "Females (Nf)", 1000, min = 1, max = 1e7, step = 100),
      shiny::numericInput(ns("Nm"), "Males (Nm)", 1000, min = 1, max = 1e7, step = 100),
      shiny::numericInput(ns("Nk1"), if (cross) "Episode 1 offspring (Nk1)" else "Offspring (Nk)", 2000, min = 2, step = 100),
      if (cross) shiny::numericInput(ns("Nk2"), "Episode 2 offspring (Nk2)", 2000, min = 1, step = 100),
      shiny::h4("Sampling offspring"),
      shiny::numericInput(ns("Sk1"), if (cross) "Episode 1 sample (Sk1)" else "Sampled offspring (Sk)", 200, min = if (cross) 1 else 2, step = 10),
      if (cross) shiny::numericInput(ns("Sk2"), "Episode 2 sample (Sk2)", 200, min = 1, step = 10),
      shiny::h4("Replication"),
      shiny::numericInput(ns("num_reps"), "Replicates", 100, min = 1, max = limits$reps, step = 10),
      shiny::numericInput(ns("seed"), "First replicate seed", 1001, min = 1, step = 1),
      shiny::actionButton(ns("run_trial"), "Run and save trial", class = "btn-primary"),
      shiny::p(class = "hint", "The next seed advances after each successful trial. Samples, summaries, settings, and seeds stay in this session; populations can be regenerated."),
      shiny::hr(), shiny::p(class = "hint", "Shared POP + HSP session limits: 30 trials, 40 million sampled individuals. Each trial allows 20 million samples across reps; each replicate allows 10 million adults plus offspring across its episodes."),
      shiny::downloadButton(ns("download_trials"), "Save complete trials (RDS)"), shiny::br(), shiny::br(),
      shiny::downloadButton(ns("download_register"), "Save trial table (CSV)")),
    shiny::mainPanel(width = 9,
      if (lucky) .pop_gamma_section(ns),
      if (scenario == "size") .pop_size_section(ns, fixed_output_cv = TRUE, episode_growth = TRUE),
      shiny::div(class = "section", shiny::h3("Explore a saved trial"),
        shiny::p("Changing preparation controls does not alter a saved trial. Select a trial and replicate below."),
        shiny::selectInput(ns("trial_id"), "Saved trial", choices = character()),
        shiny::uiOutput(ns("saved_settings")), shiny::uiOutput(ns("rep_control")), shiny::textOutput(ns("rep_details"))),
      shiny::div(class = "section", shiny::h3("01 \u00b7 Population / Reproduction"),
        shiny::p("Realized offspring counts include every adult, including those with no offspring."),
        if (scenario == "size") shiny::tagList(
          shiny::radioButtons(ns("trait_episode"), "Show adult traits in", c("Episode 1" = "1", "Episode 2" = "2"), inline = TRUE),
          shiny::plotOutput(ns("population_growth"), height = "300px"),
          shiny::plotOutput(ns("population_output"), height = "300px"),
          shiny::actionButton(ns("erro_y_toggle"), "ERRO density y-axis: linear (switch to log10)"),
          shiny::plotOutput(ns("population_erro"), height = "300px")),
        if (cross) shiny::tagList(shiny::h4("Does an adult's advantage persist?"),
          shiny::radioButtons(ns("persistence_view"), "Compare episodes using",
            c("ERRO" = "weights", "Realized offspring" = "counts", if (scenario == "size") c("Length" = "lengths")), inline = TRUE),
          shiny::plotOutput(ns("persistence_plot"), height = "310px")),
        shiny::actionButton(ns("population_y_toggle"), "Y-axis: linear (switch to log10)"),
        shiny::plotOutput(ns("population_plot"), height = if (cross) "470px" else "290px"),
        shiny::h4(if (cross) "Realized cross-cohort effective numbers" else "Realized effective numbers of parents"),
        shiny::tableOutput(ns("effective_table")),
        shiny::p(class = "hint", if (cross) "From the selected replicate's full offspring populations: h = sum(k_i1 * k_i2) / (Nk1 * Nk2), Effective N = 1/h. This is the reciprocal cross-cohort sharing probability, not an adult survey estimate. The ratio can exceed one; no shared parents gives Infinity." else "From the selected replicate's full offspring population: h = sum(k_i * (k_i - 1)) / (Nk * (Nk - 1)), Effective N = 1/h. These are population quantities, not sample estimates. The ratio can exceed one; no sibling pairs gives Infinity.")),
      shiny::div(class = "section", shiny::h3("02 \u00b7 Sampling"),
        shiny::p("Offspring are sampled uniformly without replacement within each cohort. No adults are survey-sampled. Family-size histograms use simulation truth to show where the sampled offspring came from, including census adults with zero sampled offspring."),
        shiny::tableOutput(ns("sampling_table")),
        shiny::actionButton(ns("sampling_y_toggle"), "Y-axis: linear (switch to log10)"),
        shiny::plotOutput(ns("sampling_plot"), height = if (cross) "470px" else "290px"),
        if (scenario == "size") shiny::tagList(shiny::h4("Parents represented by sampled offspring"),
          shiny::plotOutput(ns("sample_traits"), height = "310px"))),
      shiny::div(class = "section", shiny::h3("03 \u00b7 Pairs"),
        shiny::p(if (cross) "Only comparisons between episodes are used: C = Sk1 * Sk2. Within-cohort pairs are not counted in this trial." else "Each unordered pair of distinct sampled offspring is compared once: C = Sk * (Sk - 1) / 2."),
        shiny::p("A family with k sampled offspring contributes k * (k - 1) / 2 within-cohort matches. Across cohorts it contributes k1 * k2. Full siblings contribute one maternal match and one paternal match."),
        shiny::tableOutput(ns("match_table")), shiny::plotOutput(ns("pairs_plot"), height = "330px"),
        if (scenario == "size") shiny::tagList(shiny::h4("Which adults connect the sibling pairs?"),
          shiny::radioButtons(ns("trait_variable"), "Compare parent traits", c("ERRO" = "reproductive_weight", "Length" = "length"), inline = TRUE),
          shiny::plotOutput(ns("match_traits"), height = "330px"),
          shiny::p(class = "hint", "Traits in these diagnostic histograms refer to episode 1. Each parent is represented once per sibling match. Expected enrichment is proportional to the product of its two episodes' ERROs; without growth this reduces to ERRO squared."))),
      shiny::div(class = "section", shiny::h3("04 \u00b7 Estimation"),
        shiny::p("The naive model uses E(matches) = C / N for each sex. Both the moment estimator and continuous Poisson pseudo-MLE give Nhat = C / observed matches. This targets census adult abundance; zero matches give an unbounded estimate."),
        shiny::p("For equal reproduction the sharing probability is 1/N. Within-cohort reproductive inequality can elevate it. Independent seasonal luck averages out in cross-cohort comparisons; persistent size-based differences do not. Finite-sample inverse estimates can also be biased even when the expected match count is correct."),
        shiny::plotOutput(ns("estimates_plot"), height = "460px"),
        shiny::h4("Approximate 95% confidence intervals"),
        shiny::p("The same Poisson Wald approximation as in the POP app is used: SE = Nhat / sqrt(matches), limits Nhat plus or minus 1.96 SE. Lines are sorted by lower then upper limit; the selected replicate is thicker. These working intervals can fail when matches cluster in families. Negative lower limits are retained; zero matches give no interval."),
        shiny::plotOutput(ns("intervals_plot"), height = "420px"),
        shiny::div(style = "overflow-x:auto;", shiny::tableOutput(ns("estimation_table"))),
        shiny::p(class = "hint", "Median, mean, and CV exclude unbounded estimates. CI inclusion counts finite intervals containing census N, conditional on finite estimates. This is not unconditional coverage.")),
      shiny::div(class = "section", shiny::h3("Saved trials"), shiny::textOutput(ns("storage_details")),
        shiny::div(class = "trial-scroll", shiny::tableOutput(ns("trial_register")))),
      .hsp_questions(scenario))))
}

#' Questions specific to the sibling teaching progression
#' @param scenario Scenario key.
#' @return A prose section beneath the saved-trial table.
#' @keywords internal
.hsp_questions <- function(scenario) {
  questions <- switch(scenario,
    equal = c("Increase census abundance by 100 and offspring sample size by 10. Does the expected number of sibling matches stay approximately the same?",
      "Aim for about 50 matches per sex. How do finite CV and confidence-interval inclusion compare with the POP app?",
      "Inspect full-sibling counts. Can you explain why a full sibling contributes one match to each sex's count?"),
    lucky = c("Change the gamma shape for only one sex. Which direction creates more unequal family sizes?",
      "Compare estimates with census N and realized Effective N. What quantity is the naive sibling estimator tracking?",
      "Do the Poisson intervals still include census abundance as often as you would expect?"),
    cross_lucky = c("Use a small gamma shape. Compare the episode-to-episode ERRO cloud with the Lucky litters tab. Does an adult's advantage persist?",
      "Compare census estimates from within-cohort and cross-cohort trials with similar expected match counts.",
      "Does independent luck remove variability as well as its systematic effect on the sharing probability? Explore population size and sample fraction."),
    size = c("Increase the female length exponent while leaving the male exponent alone. Do cross-cohort comparisons still recover census abundance?",
      "Compare adults represented by sampled offspring with those represented by sibling matches. Why does the second group show stronger enrichment?",
      "Set both length exponents to zero. What happens to ERRO, Effective N, and the estimates?"))
  shiny::div(class = "section", shiny::h3("Questions to Ask / Things to Try"),
    do.call(shiny::tags$ol, lapply(questions, shiny::tags$li)))
}

#' Validate sibling trial requests against shared POP/HSP storage limits
#' @param settings Validated sibling settings.
#' @param num_reps Number of replicates.
#' @param usage Shared reactive session counters.
#' @return Invisibly returns the new trial's retained offspring sample count.
#' @keywords internal
.hsp_check_budget <- function(settings, num_reps, usage) {
  .pop_check_budget(list(Nf = settings$Nf, Nm = settings$Nm, Nk = settings$Nk1 + as.double(settings$Nk2),
    Sf = 0, Sm = 0, Sk = settings$Sk1 + as.double(settings$Sk2)), num_reps,
    .pop_usage_total(usage, "samples"), .pop_usage_total(usage, "trials"))
}

#' Bind one sibling scenario's immutable trial registry and plots
#' @param id Shiny module ID.
#' @param scenario Scenario key.
#' @param shared_usage Session-wide counters shared with POP and other sibling tabs.
#' @return Reactive handles for trials and selected replicate, useful in tests.
#' @keywords internal
.hsp_scenario_server <- function(id, scenario, shared_usage) {
  shiny::moduleServer(id, function(input, output, session) {
    trials <- shiny::reactiveVal(list())
    current_trial <- shiny::reactive({
      saved <- trials(); shiny::req(length(saved) > 0)
      key <- input$trial_id
      if (length(key) != 1L || is.na(key) || !key %in% names(saved)) key <- tail(names(saved), 1)
      saved[[key]]
    })
    rep_id <- shiny::reactive({
      n <- current_trial()$num_reps
      if (is.null(input$rep_id)) 1L else max(1L, min(n, as.integer(input$rep_id)))
    })
    current_rep <- shiny::reactive(current_trial()$reps[[rep_id()]])
    cross <- scenario %in% c("cross_lucky", "size")
    lucky <- scenario %in% c("lucky", "cross_lucky")
    shiny::observeEvent(input$run_trial, {
      tryCatch({
        s <- hsp_settings(Nf = input$Nf, Nm = input$Nm, Nk1 = input$Nk1,
          Nk2 = if (cross) input$Nk2 else 0, Sk1 = input$Sk1, Sk2 = if (cross) input$Sk2 else 0,
          scenario = scenario, use_weights = if (lucky) isTRUE(input$use_weights) else FALSE,
          female_shape = if (lucky) input$female_shape else 1, male_shape = if (lucky) input$male_shape else 1,
          size_model = if (scenario == "size") .pop_size_inputs(input) else NULL,
          growth = if (scenario == "size") isTRUE(input$growth) else FALSE,
          growth_cv = if (scenario == "size") input$growth_cv else 0)
        .hsp_check_budget(s, input$num_reps, shared_usage)
        id <- length(trials()) + 1L
        trial <- shiny::withProgress(message = paste("Simulating sibling trial", id), value = 0,
          hsp_run_trial(s, input$num_reps, input$seed, id,
            progress = function(done, total) shiny::setProgress(value = done / total)))
        saved <- trials(); saved[[as.character(id)]] <- trial; trials(saved)
        shared_usage[[paste0("trials_hs_", scenario)]] <- length(saved)
        shared_usage[[paste0("samples_hs_", scenario)]] <- sum(vapply(saved,
          function(t) (t$settings$Sk1 + as.double(t$settings$Sk2)) * t$num_reps, numeric(1)))
        labels <- vapply(saved, function(t) paste("Trial", t$trial_id, "|", t$num_reps, "reps | seeds", t$seed,
          "to", as.double(t$seed) + t$num_reps - 1), character(1))
        shiny::updateSelectInput(session, "trial_id", choices = stats::setNames(names(saved), labels), selected = as.character(id))
        shiny::updateNumericInput(session, "seed", value = min(.Machine$integer.max, as.double(trial$seed) + trial$num_reps))
      }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 12))
    })
    output$rep_control <- shiny::renderUI({
      n <- current_trial()$num_reps
      if (n == 1) shiny::numericInput(session$ns("rep_id"), "Replicate to inspect", 1, min = 1, max = 1)
      else shiny::sliderInput(session$ns("rep_id"), "Replicate to inspect", min = 1, max = n, value = 1, step = 1)
    })
    output$saved_settings <- shiny::renderUI({
      t <- current_trial(); s <- t$settings
      shiny::div(class = "saved-settings", shiny::strong(paste("Saved trial", t$trial_id)),
        shiny::p(paste("Nf", s$Nf, "| Nm", s$Nm, "| offspring", s$Nk1,
          if (cross) paste("and", s$Nk2) else "", "| sampled", s$Sk1,
          if (cross) paste("and", s$Sk2) else "")),
        shiny::p(if (scenario == "size") paste(
          if (isTRUE(s$growth)) paste0("One year of growth; growth-rate CV = ", s$growth_cv, ".") else "Adult traits fixed between episodes.",
          "Output CV = 0 in both episodes. Full settings are in Saved trials.") else
          if (s$use_weights) paste("Gamma shapes: female", s$female_shape, ", male", s$male_shape,
            if (cross) "| independent redraw in episode 2" else "") else "Equal expected reproduction."))
    })
    output$rep_details <- shiny::renderText({
      r <- current_rep(); paste("Replicate", rep_id(), "| seed", r$seed,
        "| maternal matches", r$matches[1], "| paternal matches", r$matches[2])
    })
    scale <- function(x) if (!is.null(x) && x %% 2 == 1) "log10" else "linear"
    shiny::observe({
      for (key in c("population_y_toggle", "sampling_y_toggle")) {
        value <- scale(input[[key]])
        shiny::updateActionButton(session, key, label = paste0("Y-axis: ", value, " (switch to ", if (value == "linear") "log10" else "linear", ")"))
      }
    })
    output$population_plot <- shiny::renderPlot(hsp_plot_reproduction(current_rep(), y_scale = scale(input$population_y_toggle)))
    output$sampling_plot <- shiny::renderPlot(hsp_plot_reproduction(current_rep(), sampled = TRUE, y_scale = scale(input$sampling_y_toggle)))
    output$persistence_plot <- shiny::renderPlot(hsp_plot_persistence(current_rep(),
      realized = identical(input$persistence_view, "counts"), lengths = identical(input$persistence_view, "lengths")))
    output$pairs_plot <- shiny::renderPlot(hsp_plot_pairs(current_trial(), rep_id()))
    output$estimates_plot <- shiny::renderPlot(pop_plot_estimates(current_trial(), rep_id()))
    output$intervals_plot <- shiny::renderPlot(pop_plot_intervals(current_trial(), rep_id()))
    output$effective_table <- shiny::renderTable({
      x <- current_rep()$effective
      data.frame(Sex = x$sex, `Census N` = x$census_n, `Effective N` = signif(x$effective_n, 5),
        `Effective N / N` = signif(x$effective_ratio, 4), check.names = FALSE)
    }, striped = TRUE, digits = 4)
    output$sampling_table <- shiny::renderTable({
      s <- current_trial()$settings
      data.frame(Episode = seq_len(if (cross) 2 else 1),
        Population = if (cross) c(s$Nk1, s$Nk2) else s$Nk1,
        Sample = if (cross) c(s$Sk1, s$Sk2) else s$Sk1,
        `Sample fraction` = sprintf("%.1f%%", 100 * if (cross) c(s$Sk1/s$Nk1, s$Sk2/s$Nk2) else s$Sk1/s$Nk1), check.names = FALSE)
    }, striped = TRUE)
    output$match_table <- shiny::renderTable({
      r <- current_rep()
      data.frame(Comparisons = r$comparisons, `Maternal matches` = r$matches[1],
        `Paternal matches` = r$matches[2], `Full-sibling pairs (included in both)` = r$full_siblings, check.names = FALSE)
    }, striped = TRUE)
    output$estimation_table <- shiny::renderTable({
      x <- pop_summarize_estimates(current_trial()$summary)
      data.frame(Sex = x$sex, Method = x$method, Truth = x$truth,
        `Finite median` = round(x$finite_median, 1), `Finite mean` = round(x$finite_mean, 1),
        `Finite CV` = ifelse(is.finite(x$finite_cv), sprintf("%.1f%%", 100*x$finite_cv), "Unavailable"),
        `Finite 95% CI inclusion` = ifelse(x$method == "Poisson pseudo-MLE",
          paste0(x$num_ci_included, " / ", x$num_finite, " (",
            ifelse(is.finite(x$finite_ci_inclusion), sprintf("%.1f%%", 100*x$finite_ci_inclusion), "unavailable"), ")"), "-"),
        Unbounded = paste(x$num_unbounded, "/", x$num_reps), check.names = FALSE)
    }, striped = TRUE)
    output$trial_register <- shiny::renderTable(hsp_trial_table(trials()), striped = TRUE, digits = 3)
    output$storage_details <- shiny::renderText(paste(length(trials()), "trials in this tab |",
      .pop_usage_total(shared_usage, "trials"), "of 30 shared trial slots |",
      format(.pop_usage_total(shared_usage, "samples"), big.mark = ","), "sampled individuals across POP + HSP"))
    output$download_trials <- shiny::downloadHandler(filename = function() paste0("ckmr-siblings-", scenario, ".rds"),
      content = function(file) saveRDS(trials(), file))
    output$download_register <- shiny::downloadHandler(filename = function() paste0("ckmr-siblings-", scenario, ".csv"),
      content = function(file) utils::write.csv(hsp_trial_table(trials()), file, row.names = FALSE))
    if (lucky) .pop_gamma_server(input, output, session)
    if (scenario == "size") {
      .pop_size_server(input, output, session, current_rep)
      episode_traits <- shiny::reactive(.hsp_episode_traits(current_rep(), identical(input$trait_episode, "2")))
      output$population_growth <- shiny::renderPlot({
        x <- episode_traits(); pop_plot_size_relationship(x$adults, x$model, "growth", x$summary)
      })
      output$population_output <- shiny::renderPlot({
        x <- episode_traits(); pop_plot_size_relationship(x$adults, x$model, "output", x$summary)
      })
      output$population_erro <- shiny::renderPlot(pop_plot_erro(episode_traits()$summary, scale(input$erro_y_toggle)))
      output$sample_traits <- shiny::renderPlot(hsp_plot_traits(current_rep(), matched = FALSE, variable = "length"))
      output$match_traits <- shiny::renderPlot(hsp_plot_traits(current_rep(), variable =
        if (identical(input$trait_variable, "length")) "length" else "reproductive_weight"))
    }
    list(trials = trials, current_trial = current_trial, current_rep = current_rep)
  })
}

#' Start all sibling registries with independent trial stores and shared budgets
#' @param input,output,session Standard Shiny server arguments.
#' @param shared_usage Optional shared POP/HSP budget counters.
#' @return Named reactive handles for the four sibling scenarios.
#' @export
hsp_app_server <- function(input, output, session, shared_usage = NULL) {
  if (is.null(shared_usage)) shared_usage <- .pop_session_usage()
  stats::setNames(lapply(names(.hsp_titles()), function(key)
    .hsp_scenario_server(paste0("hs_", key), key, shared_usage)), names(.hsp_titles()))
}
