#' Storage and simulation limits shared by the UI and server
#' @return A named list of count limits. Population size applies to one replicate;
#' sample limits count Sf + Sm + Sk across replicates and then across trials.
#' @keywords internal
.pop_app_limits <- function() {
  list(trial_samples = 2e7, session_samples = 4e7, replicate_population = 1e7,
       trials = 30, reps = 500)
}

#' Validate an app request against transient and retained-data budgets
#' @param settings Validated simulation settings.
#' @param num_reps Number of requested replicates.
#' @param stored_samples Sampled individuals already retained in the session.
#' @param num_trials Number of saved trials.
#' @param limits Limits from .pop_app_limits(); injectable for boundary checks.
#' @return Invisibly returns the number of sampled individuals in the new trial.
#' No allocation or simulation occurs when checking budgets.
#' @keywords internal
.pop_check_budget <- function(settings, num_reps, stored_samples = 0,
                               num_trials = 0, limits = .pop_app_limits()) {
  .pop_integer(num_reps, "Replicates")
  label <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (num_reps > limits$reps) stop("Use at most ", limits$reps, " replicates per trial.", call. = FALSE)
  population_size <- sum(unlist(settings[c("Nf", "Nm", "Nk")]))
  samples_size <- sum(unlist(settings[c("Sf", "Sm", "Sk")])) * num_reps
  if (population_size > limits$replicate_population) {
    stop("One replicate would exceed ", label(limits$replicate_population),
         " population individuals. Reduce population sizes.", call. = FALSE)
  }
  if (samples_size > limits$trial_samples) {
    stop("This trial exceeds ", label(limits$trial_samples),
         " sampled individuals across replicates. Reduce sample sizes or replicates.", call. = FALSE)
  }
  if (num_trials >= limits$trials) {
    stop("This session already has ", limits$trials,
         " trials. Save the results before starting a new session.", call. = FALSE)
  }
  if (stored_samples + samples_size > limits$session_samples) {
    stop("This trial would exceed the session limit of ", label(limits$session_samples),
         " sampled individuals. Save the results before starting a new session.", call. = FALSE)
  }
  invisible(samples_size)
}

#' Build the interface for the parent-offspring scenarios
#' @return A Shiny UI object.
#' @keywords internal
pop_app_ui <- function() {
  shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML("
      body { background:#f5f5ef; color:#233a3d; font-size:15px; }
      .container-fluid { max-width:1480px; }
      .hero { padding:24px 8px 18px; border-bottom:3px solid #237e83; margin-bottom:22px; }
      .hero h1 { font-size:32px; font-weight:700; margin:5px 0 10px; }
      .eyebrow { text-transform:uppercase; letter-spacing:2px; color:#237e83; font-size:12px; font-weight:700; }
      .well { background:#fff; border:1px solid #dce3df; box-shadow:none; border-radius:10px; }
      .section { background:white; border:1px solid #dce3df; border-radius:10px; padding:18px 22px; margin:18px 0; }
      .section h3 { margin-top:0; font-size:22px; }
      .section p, .hint { color:#55696c; }
      .btn-primary { background:#237e83; border-color:#237e83; }
      .btn-primary:hover { background:#195e62; border-color:#195e62; }
      .saved-settings { background:#e9f3f1; border-radius:6px; padding:14px; margin-bottom:16px; }
      .trial-scroll { overflow:auto; max-height:340px; font-size:13px; }
      .trial-scroll table { white-space:nowrap; }
      .trial-scroll th { position:sticky; top:0; background:#e9f3f1; }
      .table { margin-bottom:8px; }
      .effective-numbers .table.shiny-table.spacing-xs > thead > tr > th,
      .effective-numbers .table.shiny-table.spacing-xs > tbody > tr > td { padding:4px; }
      .nav-tabs { margin-bottom:16px; }
      .shiny-output-error-validation { color:#55696c; }
    "))),
    shiny::div(class = "hero", shiny::div(class = "eyebrow", "CKMR \u00b7 parent-offspring pair fundamentals"),
      shiny::h1("An increasingly interesting series of simple parent-offspring scenarios"),
      shiny::p("Simulate a population, take a sample, and explore what kin matches tell us about abundance.")),
    shiny::tabsetPanel(id = "scenario",
      shiny::tabPanel("1 \u00b7 Equal reproduction", .pop_scenario_ui()),
      shiny::tabPanel("Randomly unequal reproduction",
        .pop_scenario_ui(ns = shiny::NS("unequal"), unequal = TRUE)),
      shiny::tabPanel("Size-based reproduction",
        .pop_scenario_ui(ns = shiny::NS("size"), size_based = TRUE))
    )
  )
}

#' Build the shared section layout with scenario-specific additions
#' @param ns Namespace function; identity for tab one, a Shiny namespace for tab two.
#' @param unequal Whether to include the gamma introduction, controls, and section 00.
#' @param size_based Include growth-based reproduction controls and trait figures.
#' @return A tag list with the scenario introduction, sidebar, and sections 01-04.
#' @keywords internal
.pop_scenario_ui <- function(ns = identity, unequal = FALSE, size_based = FALSE) {
  limits <- .pop_app_limits()
  shiny::tagList(
        if (size_based) .pop_size_intro() else if (unequal) .pop_unequal_intro() else shiny::div(class = "section",
          shiny::h3("Start with the simplest possible population"),
          shiny::p("Imagine a population of females and males producing a single cohort of offspring. Each offspring has one mother and one father, chosen independently at random. All females have the same chance of being a mother, and all males have the same chance of being a father. Actual family sizes vary by chance, even though expected reproductive output is equal within each sex. There are no age or size differences, and no adults die between reproduction and sampling."),
          shiny::p("Take random samples of females, males, and offspring, then identify every parent-offspring pair within those samples. Use those matches to estimate female and male abundance separately. With sample sizes held fixed, a match is more likely in a smaller population."),
          shiny::p("Choose population and sample sizes to create a trial, then repeat the experiment in independent replicates. Follow one replicate through reproduction and sampling, and use the trial-level figures to explore how chance and sampling effort affect pair counts and the precision of abundance estimates.")
        ),
        shiny::sidebarLayout(
          shiny::sidebarPanel(width = 3,
            shiny::h3("Prepare a trial"),
            if (unequal) shiny::checkboxInput(ns("use_weights"), "Use random reproductive sampling weights", value = TRUE),
            shiny::p(class = "hint", "Each replicate creates a new population. These controls apply to the next trial."),
            shiny::h4("Population"),
            shiny::numericInput(ns("Nf"), "Females (Nf)", 1000, min = 1, max = limits$replicate_population, step = 100),
            shiny::numericInput(ns("Nm"), "Males (Nm)", 1000, min = 1, max = limits$replicate_population, step = 100),
            shiny::numericInput(ns("Nk"), "Offspring (Nk)", 2000, min = 1, max = limits$replicate_population, step = 100),
            shiny::h4("Sampling"),
            shiny::numericInput(ns("Sf"), "Sampled females (Sf)", 100, min = 1, step = 10),
            shiny::numericInput(ns("Sm"), "Sampled males (Sm)", 100, min = 1, step = 10),
            shiny::numericInput(ns("Sk"), "Sampled offspring (Sk)", 200, min = 1, step = 10),
            shiny::h4("Replication"),
            shiny::numericInput(ns("num_reps"), "Replicates", 100, min = 1, max = limits$reps, step = 10),
            shiny::numericInput(ns("seed"), "First replicate seed", 1001, min = 1, step = 1),
            shiny::p(class = "hint", "Replicate r uses seed + r \u2212 1. The next seed advances after a successful trial."),
            shiny::actionButton(ns("run_trial"), "Run and save trial", class = "btn-primary"),
            shiny::hr(),
            shiny::p(class = "hint", paste0("Complete samples and population summaries are kept in this session. Limits: ",
              limits$trials, " trials; ", limits$trial_samples / 1e6,
              " million sampled individuals across reps per trial; ", limits$session_samples / 1e6,
              " million sampled individuals per session; ", limits$replicate_population / 1e6,
              " million population individuals in any one replicate. Session limits are shared across all tabs.")),
            shiny::downloadButton(ns("download_trials"), "Save complete trials (RDS)"),
            shiny::br(), shiny::br(),
            shiny::downloadButton(ns("download_register"), "Save trial table (CSV)")
          ),
          shiny::mainPanel(width = 9,
            if (unequal) .pop_gamma_section(ns),
            if (size_based) .pop_size_section(ns),
            shiny::div(class = "section",
              shiny::h3("Explore a saved trial"),
              shiny::p("Run a trial to begin. Select a saved trial and a replicate to inspect; changing the controls on the left does not change saved results."),
              shiny::selectInput(ns("trial_id"), "Saved trial", choices = character()),
              shiny::uiOutput(ns("saved_settings")),
              shiny::uiOutput(ns("rep_control")),
              shiny::textOutput(ns("rep_details"))
            ),
            shiny::div(class = "section",
              shiny::h3("01 \u00b7 Population / Reproduction"),
              shiny::p(if (size_based) "Lengths and reproductive propensities vary among adults. Points show a random display subset of at most 500 adults per sex. ERRO histograms summarize every adult; realized offspring counts show the subsequent reproduction draw." else if (unequal) "Each offspring independently receives one mother and father, chosen using their reproductive sampling weights when enabled. These plots show realized offspring counts in the saved replicate, including adults with no offspring." else "Every adult has the same expected reproductive output. Each offspring independently receives one randomly chosen mother and father. Realized family sizes still vary."),
              if (size_based) shiny::tagList(
                shiny::plotOutput(ns("population_growth"), height = "300px"),
                shiny::plotOutput(ns("population_output"), height = "300px"),
                shiny::h4("Assigned expected relative reproductive outputs"),
                shiny::actionButton(ns("erro_y_toggle"), "ERRO density y-axis: linear (switch to log10)"),
                shiny::plotOutput(ns("population_erro"), height = "300px"),
                shiny::h4("Realized offspring counts")
              ),
              shiny::actionButton(ns("population_y_toggle"), "Y-axis: linear (switch to log10)"),
              shiny::plotOutput(ns("population_plot"), height = "290px"),
              shiny::h4("Realized effective numbers of parents"),
              shiny::div(class = "effective-numbers", style = "overflow-x:auto; font-size:13px;", shiny::tableOutput(ns("effective_numbers_table"))),
              shiny::p(class = "hint",
                shiny::strong("Calculated from all offspring in the selected replicate's population, not just the sample. "),
                "Effective N is the reciprocal of the probability that two distinct offspring drawn without replacement share a mother (female row) or a father (male row). The ratio divides Effective N by census N for that sex. This finite-cohort effective number can exceed census N. Infinity means there are no same-parent offspring pairs; values are unavailable when the population has fewer than two offspring.")
            ),
            shiny::div(class = "section",
              shiny::h3("02 \u00b7 Sampling"),
              shiny::p("Females, males, and offspring are sampled independently without replacement. Each plot counts only sampled offspring belonging to sampled parents, including parents with no matches."),
              shiny::tableOutput(ns("sampling_table")),
              if (size_based) shiny::tagList(
                shiny::p("These are the actual sampled adults, drawn uniformly regardless of age, length, or reproductive sampling weight. Axes match the population views above."),
                shiny::plotOutput(ns("sample_growth"), height = "300px"),
                shiny::plotOutput(ns("sample_output"), height = "300px")
              ),
              shiny::actionButton(ns("sampling_y_toggle"), "Y-axis: linear (switch to log10)"),
              shiny::plotOutput(ns("sampling_plot"), height = "290px")

            ),
            shiny::div(class = "section",
              shiny::h3("03 \u00b7 Pairs"),
              shiny::p("These distributions summarize all replicates in the saved trial. A parent with three sampled offspring contributes three pairs. Kin identification is assumed perfect."),
              shiny::plotOutput(ns("pairs_plot"), height = "310px"),
              if (size_based) shiny::tagList(
                shiny::h4("Which adults appear in identified pairs?"),
                shiny::p("For the selected replicate, compare normalized histograms of population adults, sampled adults, and parents represented once per identified pair. A parent with five matched offspring contributes five observations. Adult survey sampling is uniform; enrichment arises through reproduction."),
                shiny::actionButton(ns("pair_trait_toggle"), "Compare ERRO (switch to length)"),
                shiny::plotOutput(ns("pair_traits_plot"), height = "auto"),
                shiny::p(class = "hint", "The left axes show probability density. Each right axis gives identified pairs per bin for the solid outline only: density times that sex's number of pairs times bin width. Whole-number ticks show discrete pair counts. With no pairs, the count axis is omitted."),
                shiny::p(class = "hint", "The dark dashed line shows the model expectation for pair parents, averaging over reproduction and uniform adult and offspring surveys. For mean-one ERRO r with population distribution f(r), this is r * f(r). In each bin we sum actual population ERROs and divide by total ERRO and bin width; the same reproductive weighting gives the expectation for length. All distributions use identical bins within each sex. With zero identified pairs, there is no observed-pair outline.")
              )
            ),
            shiny::div(class = "section",
              shiny::h3("04 \u00b7 Estimation"),
              shiny::p("For each sex, E(#pairs) = sampled-adults \u00d7 sampled-offspring / num-adults. The moment-based estimator is num-adults = sampled-adults \u00d7 sampled-offspring / observed-number-of-pairs."),
              shiny::p("For each sex, abundance = sampled parents \u00d7 sampled offspring / observed pairs. The Poisson pseudo-likelihood uses that same expected pair count, so its continuous maximum gives the same point estimate."),
              shiny::p("Zero pairs give an unbounded estimate. Those replicates are counted in every panel, but cannot appear on a finite histogram axis. The likelihood is a working model."),
              if (unequal || size_based) shiny::p("Uniform adult sampling preserves the expected pair count, but unequal reproductive sampling weights can increase variability between replicates. These Poisson intervals do not account for that additional variability; explore their finite 95% CI inclusion."),
              shiny::plotOutput(ns("estimates_plot"), height = "460px"),
              shiny::h4("Approximate 95% confidence intervals"),
              shiny::p("Each thin line is one replicate, sorted from bottom to top by lower limit, then upper limit. The selected replicate has a thicker line and a point at its estimate. Both estimators give the same intervals here. Wald limits use the estimate plus or minus 1.96 standard errors; with few matches this approximation can be poor and lower limits may be negative. Zero matches have no Wald interval."),
              shiny::plotOutput(ns("intervals_plot"), height = "420px"),
              shiny::div(style = "overflow-x:auto;", shiny::tableOutput(ns("estimation_table"))),
              shiny::p(class = "hint", "Median, mean, and CV summarize finite estimates only. CV is the sample standard deviation divided by the mean, shown as a percentage; it is unavailable with fewer than two finite estimates. For Poisson pseudo-MLE rows, finite 95% CI inclusion reports how many finite-replicate confidence intervals contain the true abundance, divided by the number of finite replicates. A single replicate is an illustration, not a distribution of estimation performance.")
            ),
            shiny::div(class = "section",
              shiny::h3(paste0("Saved trials \u00b7 ", if (size_based) "Size-based reproduction" else if (unequal) "Randomly unequal reproduction" else "Equal reproduction")),
              shiny::p("This scenario's register preserves the settings for each trial. Scroll across for seeds and summaries; select a trial above to explore its replicates."),
              shiny::div(class = "trial-scroll", shiny::tableOutput(ns("trial_register"))),
              shiny::textOutput(ns("storage_details"))
            )
          )
        )
  )
}

#' Coordinate immutable trials and reactive views for one scenario
#'
#' The only trial mutation is appending a complete, successfully generated trial.
#' Input edits and view changes never modify an existing object. Complete samples
#' and population frequency summaries are retained in a session-local list. The
#' trial register is a derived tibble. Budgets distinguish retained sample rows
#' from the transient population generated for one replicate.
#' @param input,output,session Standard Shiny server arguments.
#' @param unequal Whether this is the gamma-reproduction scenario.
#' @param size_based Whether this is the length-based reproduction scenario.
#' @param shared_usage Session usage shared across independent scenario registries.
#' @return Used for its side effects on the Shiny session.
#' @keywords internal
.pop_scenario_server <- function(input, output, session, unequal = FALSE, shared_usage = NULL, size_based = FALSE) {
  if (is.null(shared_usage)) shared_usage <- .pop_session_usage()
  scenario_key <- if (size_based) "size" else if (unequal) "unequal" else "equal"
  ns <- session$ns
  trials <- shiny::reactiveVal(list())
  limits <- .pop_app_limits()
  # Count sampled individuals only; unsampled population rows are discarded.
  stored_individuals <- shiny::reactive({
    sum(vapply(trials(), function(x) {
      sum(unlist(x$settings[c("Sf", "Sm", "Sk")])) * x$num_reps
    }, numeric(1)))
  })
  current_trial <- shiny::reactive({
    shiny::req(length(trials()), input$trial_id)
    trial <- trials()[[as.character(input$trial_id)]]
    shiny::req(!is.null(trial))
    trial
  })
  current_rep_id <- shiny::reactive({
    shiny::req(input$rep_id)
    max(1L, min(as.integer(input$rep_id), current_trial()$num_reps))
  })
  current_rep <- shiny::reactive(current_trial()$reps[[current_rep_id()]])

  shiny::observeEvent(input$run_trial, {
    tryCatch({
      settings <- pop_settings(input$Nf, input$Nm, input$Nk, input$Sf, input$Sm, input$Sk,
        use_weights = unequal && isTRUE(input$use_weights),
        female_shape = if (unequal) input$female_shape else 1,
        male_shape = if (unequal) input$male_shape else 1,
        size_model = if (size_based) .pop_size_inputs(input) else NULL)
      .pop_integer(input$seed, "Seed")
      .pop_check_budget(settings, input$num_reps,
        shared_usage$samples_equal + shared_usage$samples_unequal + shared_usage$samples_size,
        shared_usage$trials_equal + shared_usage$trials_unequal + shared_usage$trials_size, limits)
      id <- length(trials()) + 1L
      trial <- shiny::withProgress(message = paste("Simulating trial", id), value = 0, {
        pop_run_trial(settings, input$num_reps, input$seed, id,
                       progress = function(done, total) shiny::setProgress(value = done / total))
      })
      saved <- trials()
      saved[[as.character(id)]] <- trial
      trials(saved)
      shared_usage[[paste0("samples_", scenario_key)]] <- sum(vapply(saved, function(x)
        sum(unlist(x$settings[c("Sf", "Sm", "Sk")])) * x$num_reps, numeric(1)))
      shared_usage[[paste0("trials_", scenario_key)]] <- length(saved)
      labels <- vapply(saved, function(x) {
        paste0("Trial ", x$trial_id, " \u00b7 ", x$num_reps, " reps \u00b7 Nf ", x$settings$Nf,
               " / Nm ", x$settings$Nm, " / Nk ", x$settings$Nk)
      }, character(1))
      shiny::updateSelectInput(session, "trial_id",
                                 choices = stats::setNames(names(saved), labels),
                                 selected = as.character(id))
      shiny::updateNumericInput(session, "seed", value = min(.Machine$integer.max,
                                                            as.double(trial$seed) + trial$num_reps))
    }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 10))
  })

  output$rep_control <- shiny::renderUI({
    n <- current_trial()$num_reps
    if (n == 1) {
      shiny::numericInput(ns("rep_id"), "Replicate to inspect (one replicate saved)",
                           value = 1, min = 1, max = 1, step = 1)
    } else {
      shiny::sliderInput(ns("rep_id"), "Replicate to inspect", min = 1, max = n,
                          value = 1, step = 1)
    }
  })
  output$saved_settings <- shiny::renderUI({
    x <- current_trial()
    s <- x$settings
    shiny::div(class = "saved-settings",
      shiny::strong(paste("Trial", x$trial_id, "\u00b7 saved settings")),
      shiny::p(paste0("Population: ", s$Nf, " females, ", s$Nm, " males, ", s$Nk,
                      " offspring. Samples: ", s$Sf, " females, ", s$Sm,
                      " males, ", s$Sk, " offspring.")),
      if (size_based) shiny::tagList(lapply(c("Female", "Male"), function(sex) {
        p <- s$size_model[[sex]]
        shiny::p(paste0(sex, ": ", paste(paste(names(p), unlist(p), sep = "="), collapse = "; "), "."))
      })),
      if (unequal) shiny::p(if (s$use_weights)
        paste0("Random reproductive sampling weights: Gamma(scale = 1), female shape ",
               s$female_shape, ", male shape ", s$male_shape, ".") else
        paste0("Equal reproduction: random reproductive sampling weights disabled. Saved preview shapes: female ",
               s$female_shape, ", male ", s$male_shape, ".")),
      paste0(x$num_reps, " independent population replicates \u00b7 seeds ", x$seed,
             "\u2013", as.double(x$seed) + x$num_reps - 1))
  })
  output$rep_details <- shiny::renderText({
    x <- current_rep()
    paste0("Replicate ", current_rep_id(), " | seed ", x$seed,
           " | mother-offspring pairs: ", sum(x$pairs$sex == "Female"),
           " | father-offspring pairs: ", sum(x$pairs$sex == "Male"))
  })
  # Each button affects its own view only. Its click count selects the scale;
  # changing scales never changes a saved simulation or its sampling results.
  population_y_scale <- shiny::reactive({
    if (!is.null(input$population_y_toggle) && input$population_y_toggle %% 2 == 1) "log10" else "linear"
  })
  sampling_y_scale <- shiny::reactive({
    if (!is.null(input$sampling_y_toggle) && input$sampling_y_toggle %% 2 == 1) "log10" else "linear"
  })
  shiny::observe({
    label <- if (population_y_scale() == "log10") "Y-axis: log10 (switch to linear)" else "Y-axis: linear (switch to log10)"
    shiny::updateActionButton(session, "population_y_toggle", label = label)
  })
  shiny::observe({
    label <- if (sampling_y_scale() == "log10") "Y-axis: log10 (switch to linear)" else "Y-axis: linear (switch to log10)"
    shiny::updateActionButton(session, "sampling_y_toggle", label = label)
  })
  output$population_plot <- shiny::renderPlot(pop_plot_reproduction(current_rep(), y_scale = population_y_scale()))
  output$sampling_plot <- shiny::renderPlot(pop_plot_reproduction(current_rep(), sampled = TRUE, y_scale = sampling_y_scale()))
  output$pairs_plot <- shiny::renderPlot(pop_plot_pairs(current_trial(), current_rep_id()))
  output$estimates_plot <- shiny::renderPlot(pop_plot_estimates(current_trial(), current_rep_id()))
  output$intervals_plot <- shiny::renderPlot(pop_plot_intervals(current_trial(), current_rep_id()))
  output$sampling_table <- shiny::renderTable({
    s <- current_trial()$settings
    data.frame(Group = c("Females", "Males", "Offspring"),
               Population = c(s$Nf, s$Nm, s$Nk), Sample = c(s$Sf, s$Sm, s$Sk),
               `Sample fraction` = sprintf("%.1f%%", 100 * c(s$Sf / s$Nf, s$Sm / s$Nm, s$Sk / s$Nk)),
               check.names = FALSE)
  }, striped = TRUE)
  output$effective_numbers_table <- shiny::renderTable({
    x <- pop_effective_numbers(current_rep())
    # Preserve small ratios with significant digits, rather than rounding them
    # to zero. Infinity and undefined two-offspring draws have distinct labels.
    format_value <- function(values, ratio = FALSE) vapply(values, function(value) {
      if (is.na(value)) "Unavailable"
      else if (is.infinite(value)) "Infinity"
      else format(if (ratio) signif(value, 4) else round(value, 1),
                  big.mark = ",", scientific = FALSE, trim = TRUE)
    }, character(1))
    data.frame(Sex = x$sex, `Census N` = format_value(x$census_n),
               `Effective N` = format_value(x$effective_n),
               `Effective N / N` = format_value(x$effective_ratio, ratio = TRUE),
               check.names = FALSE)
  }, striped = TRUE, spacing = "xs")
  output$estimation_table <- shiny::renderTable({
    x <- pop_summarize_estimates(current_trial()$summary)
    # Format only for display; the reusable summary preserves numeric precision.
    abundance <- function(values) vapply(values, function(value) {
      if (is.finite(value)) format(round(value, 1), big.mark = ",", scientific = FALSE)
      else "No finite estimates"
    }, character(1))
    data.frame(Sex = x$sex, Method = x$method, Truth = x$truth,
               `Finite median` = abundance(x$finite_median),
               `Finite mean` = abundance(x$finite_mean),
               `Finite CV` = ifelse(is.finite(x$finite_cv), sprintf("%.1f%%", 100 * x$finite_cv), "Unavailable"),
               `Finite 95% CI inclusion` = ifelse(x$method == "Poisson pseudo-MLE",
                 paste0(x$num_ci_included, " / ", x$num_finite, " (",
                   ifelse(is.finite(x$finite_ci_inclusion),
                          sprintf("%.1f%%", 100 * x$finite_ci_inclusion), "unavailable"), ")"), "-"),
               Unbounded = paste0(x$num_unbounded, " / ", x$num_reps), check.names = FALSE)
  }, striped = TRUE)
  output$trial_register <- shiny::renderTable({
    shiny::req(length(trials()))
    if (size_based) pop_size_trial_table(trials()) else pop_trial_table(trials(), include_weights = unequal)
  }, striped = TRUE, digits = 2)
  output$storage_details <- shiny::renderText({
    paste(length(trials()), "trials in this tab |",
          shared_usage$trials_equal + shared_usage$trials_unequal + shared_usage$trials_size, "of", limits$trials, "session trial slots used |",
          format(shared_usage$samples_equal + shared_usage$samples_unequal + shared_usage$samples_size, big.mark = ",", scientific = FALSE),
          "of", format(limits$session_samples, big.mark = ",", scientific = FALSE),
          "sampled individuals retained across all replicates.")
  })
  output$download_trials <- shiny::downloadHandler(
    filename = function() if (size_based) "ckmr-size-based-reproduction-trials.rds" else if (unequal) "ckmr-randomly-unequal-reproduction-trials.rds" else "ckmr-equal-reproduction-trials.rds",
    content = function(file) saveRDS(trials(), file)
  )
  output$download_register <- shiny::downloadHandler(
    filename = function() if (size_based) "ckmr-size-based-reproduction-trials.csv" else if (unequal) "ckmr-randomly-unequal-reproduction-trials.csv" else "ckmr-equal-reproduction-trials.csv",
    content = function(file) utils::write.csv(if (size_based) pop_size_trial_table(trials()) else pop_trial_table(trials(), include_weights = unequal), file, row.names = FALSE)
  )
  if (unequal) .pop_gamma_server(input, output, session)
  if (size_based) .pop_size_server(input, output, session, current_rep)
  invisible(list(trials = trials, current_trial = current_trial, current_rep = current_rep,
                 stored_individuals = stored_individuals))
}

#' Explain the random-reproduction scenario before its controls
#' @return A Shiny section with the learning goals and sampling distinction.
#' @keywords internal
.pop_unequal_intro <- function() {
  shiny::div(class = "section",
    shiny::h3("Some adults get lucky"),
    shiny::p("Here, adults receive randomly assigned reproductive sampling weights before producing offspring. Some adults have a much greater chance of being chosen as parents than others. There is no length, age, or other observable characteristic that tells us who will have a large family."),
    shiny::p("The reproductive sampling weights are drawn independently from gamma distributions, with a separate shape for females and males and scale fixed at 1. Small shapes produce greater reproductive inequality. Each replicate gets a new set of reproductive sampling weights, held fixed while that replicate's offspring are generated."),
    shiny::p("These are sampling weights for choosing parents during reproduction, not biological body mass. Sampling adults for the observed data is still uniform and independent of reproductive sampling weights. Explore whether similar average pair counts can conceal greater variation in pair counts, estimates, and confidence interval inclusion.")
  )
}

#' Build section 00 with inputs next to immediate gamma density previews
#' @param ns Scenario namespace function.
#' @return A Shiny section. Preview inputs configure only future trials.
#' @keywords internal
.pop_gamma_section <- function(ns) {
  shiny::div(class = "section",
    shiny::h3("00 \u00b7 Random reproductive sampling weights"),
    shiny::p("Explore the gamma distributions before running a trial. Only relative reproductive sampling weights within each sex matter: multiplying every female's sampling weight, or every male's sampling weight, by the same positive constant leaves reproductive probabilities unchanged. The gamma scale is fixed at 1."),
    shiny::p("With shape a, the raw mean reproductive sampling weight is a and its CV is 1/sqrt(a). Adjust the shapes below; these curves update immediately, while saved trials stay unchanged. Uncheck the sidebar option to use equal reproduction for the next trial."),
    shiny::textOutput(ns("gamma_status")),
    shiny::actionButton(ns("gamma_y_toggle"), "Density y-axis: linear (switch to log10)"),
    shiny::fluidRow(
      shiny::column(6,
        shiny::h4(style = "color:#237e83;", "Females"),
        shiny::numericInput(ns("female_shape"), "Female gamma shape", value = 1, min = 0.05, max = 1000, step = 0.05),
        shiny::textOutput(ns("female_gamma_summary")),
        shiny::plotOutput(ns("female_gamma_plot"), height = "250px")
      ),
      shiny::column(6,
        shiny::h4(style = "color:#bd7330;", "Males"),
        shiny::numericInput(ns("male_shape"), "Male gamma shape", value = 1, min = 0.05, max = 1000, step = 0.05),
        shiny::textOutput(ns("male_gamma_summary")),
        shiny::plotOutput(ns("male_gamma_plot"), height = "250px")
      )
    ),
    shiny::p(class = "hint", "Curves show theoretical raw sampling weight densities, not the realized adults in a saved trial. Each curve extends to its 99.5th percentile. For shape below 1, density diverges at zero; the curve starts just above zero. Density panels may use different axis ranges.")
  )
}

#' Bind live gamma previews without mutating saved trial data
#' @param input,output,session Shiny arguments in the second tab's module scope.
#' @return Used for its side effects on preview outputs.
#' @keywords internal
.pop_gamma_server <- function(input, output, session) {
  shapes <- shiny::reactive({
    valid <- tryCatch({
      .pop_shape(input$female_shape)
      .pop_shape(input$male_shape)
      TRUE
    }, error = function(e) FALSE)
    shiny::validate(shiny::need(valid, "Choose gamma shapes between 0.05 and 1000."))
    list(female = input$female_shape, male = input$male_shape)
  })
  density_scale <- shiny::reactive({
    if (!is.null(input$gamma_y_toggle) && input$gamma_y_toggle %% 2 == 1) "log10" else "linear"
  })
  shiny::observe({
    label <- if (density_scale() == "log10") "Density y-axis: log10 (switch to linear)" else "Density y-axis: linear (switch to log10)"
    shiny::updateActionButton(session, "gamma_y_toggle", label = label)
  })
  output$gamma_status <- shiny::renderText({
    if (isTRUE(input$use_weights)) "Next trial: use the random reproductive sampling weights shown below."
    else "Next trial: equal reproduction. Gamma previews remain available for exploration."
  })
  output$female_gamma_summary <- shiny::renderText({
    sprintf("Raw mean = %g; theoretical CV = %.1f%%", shapes()$female, 100 / sqrt(shapes()$female))
  })
  output$male_gamma_summary <- shiny::renderText({
    sprintf("Raw mean = %g; theoretical CV = %.1f%%", shapes()$male, 100 / sqrt(shapes()$male))
  })
  output$female_gamma_plot <- shiny::renderPlot(pop_plot_gamma(shapes()$female, "Female", density_scale()))
  output$male_gamma_plot <- shiny::renderPlot(pop_plot_gamma(shapes()$male, "Male", density_scale()))
}

#' Create usage counters shared by all scenario registries
#' @return Session-local reactive values for retained samples and trial slots.
#' @keywords internal
.pop_session_usage <- function() {
  shiny::reactiveValues(samples_equal = 0, samples_unequal = 0, samples_size = 0,
                        trials_equal = 0L, trials_unequal = 0L, trials_size = 0L)
}

#' Register the namespaced second scenario
#' @param id Shiny module ID.
#' @param shared_usage Optional session usage counters; new counters in tests.
#' @return Reactive handles returned by the shared scenario server.
#' @keywords internal
.pop_unequal_server <- function(id, shared_usage = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    .pop_scenario_server(input, output, session, unequal = TRUE, shared_usage = shared_usage)
  })
}

#' Run all scenario registries with shared session storage safeguards
#' @param input,output,session Standard Shiny server arguments.
#' @return Reactive handles for each tab, useful for integration checks.
#' @keywords internal
pop_app_server <- function(input, output, session) {
  usage <- .pop_session_usage()
  equal <- .pop_scenario_server(input, output, session, shared_usage = usage)
  unequal <- .pop_unequal_server("unequal", shared_usage = usage)
  size <- .pop_size_scenario_server("size", shared_usage = usage)
  invisible(list(equal = equal, unequal = unequal, size = size, usage = usage))
}

#' Register the namespaced size-based scenario
#' @param id Shiny module ID.
#' @param shared_usage Optional shared session budgets.
#' @return Reactive handles for this scenario's saved trials and selection.
#' @keywords internal
.pop_size_scenario_server <- function(id, shared_usage = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    .pop_scenario_server(input, output, session, shared_usage = shared_usage, size_based = TRUE)
  })
}

#' Introduce length-based reproduction with uniform adult survey sampling
#' @return A scenario introduction section.
#' @keywords internal
.pop_size_intro <- function() {
  shiny::div(class = "section", shiny::h3("Bigger adults have more offspring"),
    shiny::p("Adults now differ in age, growth, and length. A simple von Bertalanffy growth model and a length-based reproductive relationship turn these differences into expected relative reproductive outputs (ERROs). Larger adults tend to contribute more offspring, with separate settings for females and males."),
    shiny::p("Adult ages follow an exponential survival distribution between the chosen minimum and maximum ages. Mortality shapes this standing age distribution; no adults die between reproduction and sampling. Individual growth parameters can vary together, and reproductive output can vary even among fish of the same length."),
    shiny::p("Adult survey sampling is still uniform and independent of size, age, and reproductive sampling weight. Explore how reproductive inequality affects kin pairs and estimation, and why parents in identified pairs can be larger than typical adults even under a random survey. Size-selective sampling will come in a later scenario."),
    shiny::p("Although the sampled adults span different ages within an age-structured population, we still focus on a single cohort of offspring, all produced together and sampled at the same time. Sampling is non-lethal, and we consider only parent-offspring pairs."))
}

#' Build sex-specific model controls and immediate numerical previews
#' @param ns Shiny namespace function.
#' @return Section 00 UI; previews do not modify saved trials.
#' @keywords internal
.pop_size_section <- function(ns) {
  labels <- c(Linf = "Mean Linf (length units)", K = "Mean K (per year)",
    t0 = "t0 (years)", Z = "Mortality Z (per year)", min_age = "Minimum adult age",
    max_age = "Maximum adult age", Linf_cv = "Linf CV", K_cv = "K CV",
    rho = "Correlation of log(Linf) and log(K)", exponent = "Length exponent",
    output_cv = "Reproductive output CV given length")
  mins <- c(1, 0.01, -20, 0.01, 0, 0.01, 0, 0, -1, 0, 0)
  maxs <- c(500, 3, 100, 5, 100, 200, 3, 3, 1, 8, 3)
  steps <- c(5, .05, .1, .05, 1, 1, .05, .05, .1, .25, .05)
  shiny::div(class = "section",
    shiny::h3("00 \u00b7 Growth, age, and reproductive sampling weights"),
    shiny::p("Length = Linf * (1 - exp(-K * (age - t0))). Individual Linf and K are positive lognormal traits; their controls specify arithmetic means and CVs, with correlation on the log scale. A CV of zero removes that source of variation."),
    shiny::p("Reproductive propensity = length^exponent times independent, mean-one lognormal variation. Divide by the within-sex population mean to obtain ERRO. These relative reproductive sampling weights determine parental assignments; they are not biological body mass. A CV of 0.3 means 30%."),
    shiny::fluidRow(lapply(seq_along(c("Female", "Male")), function(i) {
      sex <- c("Female", "Male")[i]
      defaults <- pop_size_parameters(exponent = if (sex == "Female") 3 else 1)
      shiny::column(6, shiny::h4(style = paste0("color:", if (sex == "Female") "#237e83" else "#bd7330", ";"), sex),
        lapply(seq_along(defaults), function(j) {
          name <- names(defaults)[j]
          shiny::numericInput(ns(paste(tolower(sex), name, sep = "_")), labels[[name]],
            defaults[[name]], min = mins[j], max = maxs[j], step = steps[j])
        }))
    })),
    shiny::textOutput(ns("size_preview_status")),
    shiny::plotOutput(ns("preview_growth"), height = "300px"),
    shiny::plotOutput(ns("preview_output"), height = "300px"),
    shiny::actionButton(ns("preview_erro_toggle"), "Preview density y-axis: linear (switch to log10)"),
    shiny::plotOutput(ns("preview_erro"), height = "300px"),
    shiny::p(class = "hint", "ERRO distributions are numerical model approximations using 12,000 fixed draws per sex, including age, growth, and reproductive variation. They are not analytic densities. Fixed preview seeds are independent of trial seeds. The maximum age truncates the distribution; it does not create a pile-up of old fish. Controls configure the next trial; all figures below use saved settings."))
}

#' Read and validate nested size-model settings from a scenario's inputs
#' @param input Shiny input list.
#' @return Named Female and Male validated parameter lists.
#' @keywords internal
.pop_size_inputs <- function(input) {
  stats::setNames(lapply(c("female", "male"), function(s) {
    p <- lapply(names(pop_size_parameters()), function(name) input[[paste(s, name, sep = "_")]])
    names(p) <- names(pop_size_parameters())
    do.call(pop_size_parameters, p)
  }), c("Female", "Male"))
}

#' Bind size previews, saved-population/sample views, and pair enrichment
#' @param input,output,session Standard Shiny module arguments.
#' @param current_rep Reactive accessor for the selected saved replicate.
#' @return Used for output bindings; all view toggles leave saved data unchanged.
#' @keywords internal
.pop_size_server <- function(input, output, session, current_rep) {
  model <- shiny::reactive({
    result <- tryCatch(.pop_size_inputs(input), error = function(e) e)
    shiny::validate(shiny::need(!inherits(result, "error"),
      if (inherits(result, "error")) conditionMessage(result) else ""))
    result
  })
  # Debounce expensive previews while a participant edits several model controls.
  preview_model <- shiny::debounce(model, 350)
  preview <- shiny::reactive(pop_size_preview(preview_model()))
  preview_summary <- shiny::reactive(pop_summarize_size(preview()))
  output$size_preview_status <- shiny::renderText({
    preview_model()
    "Preview: current controls. Trials and the population/sample figures retain their saved settings."
  })
  output$preview_growth <- shiny::renderPlot(pop_plot_size_relationship(preview(), preview_model(), "growth", curves_only = TRUE))
  output$preview_output <- shiny::renderPlot(pop_plot_size_relationship(preview(), preview_model(), "output", curves_only = TRUE))
  scale_for <- function(value) if (!is.null(value) && value %% 2 == 1) "log10" else "linear"
  shiny::observe({
    for (key in c("preview_erro_toggle", "erro_y_toggle")) {
      scale <- scale_for(input[[key]])
      shiny::updateActionButton(session, key, label = paste0(
        if (key == "preview_erro_toggle") "Preview" else "ERRO",
        " density y-axis: ", scale, " (switch to ", if (scale == "linear") "log10" else "linear", ")"))
    }
    shiny::updateActionButton(session, "pair_trait_toggle", label =
      if (scale_for(input$pair_trait_toggle) == "linear") "Compare ERRO (switch to length)" else "Compare length (switch to ERRO)")
  })
  output$preview_erro <- shiny::renderPlot(pop_plot_erro(preview_summary(), scale_for(input$preview_erro_toggle)))
  output$population_growth <- shiny::renderPlot({
    r <- current_rep()
    pop_plot_size_relationship(r$display_adults, r$settings$size_model, "growth", r$size_summary)
  })
  output$population_output <- shiny::renderPlot({
    r <- current_rep()
    pop_plot_size_relationship(r$display_adults, r$settings$size_model, "output", r$size_summary)
  })
  sampled <- shiny::reactive(rbind(current_rep()$samples$females, current_rep()$samples$males))
  output$sample_growth <- shiny::renderPlot(pop_plot_size_relationship(sampled(), current_rep()$settings$size_model, "growth", current_rep()$size_summary))
  output$sample_output <- shiny::renderPlot(pop_plot_size_relationship(sampled(), current_rep()$settings$size_model, "output", current_rep()$size_summary))
  output$population_erro <- shiny::renderPlot(pop_plot_erro(current_rep()$size_summary, scale_for(input$erro_y_toggle)))
  output$pair_traits_plot <- shiny::renderPlot({
    plot <- pop_plot_pair_traits(current_rep(),
      if (scale_for(input$pair_trait_toggle) == "linear") "reproductive_weight" else "length")
    # app.R sources helpers into its local environment. Draw explicitly rather
    # than relying on Shiny finding a globally registered S3 print method.
    print.pop_pair_traits(plot)
    invisible(NULL)
  },
    height = function() {
      width <- session$clientData[[paste0("output_", session$ns("pair_traits_plot"), "_width")]]
      if (!is.null(width) && width < 562) 700 else 390
    })
}
