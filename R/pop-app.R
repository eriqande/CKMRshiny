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

#' Build the interface for the first parent-offspring scenario
#' @return A Shiny UI object.
#' @keywords internal
pop_app_ui <- function() {
  limits <- .pop_app_limits()
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
      .nav-tabs { margin-bottom:16px; }
      .shiny-output-error-validation { color:#55696c; }
    "))),
    shiny::div(class = "hero", shiny::div(class = "eyebrow", "CKMR \u00b7 parent-offspring pair fundamentals"),
      shiny::h1("An increasingly interesting series of simple parent-offspring scenarios"),
      shiny::p("Simulate a population, take a sample, and explore what kin matches tell us about abundance.")),
    shiny::tabsetPanel(id = "scenario",
      shiny::tabPanel("1 \u00b7 Equal reproduction",
        shiny::div(class = "section",
          shiny::h3("Start with the simplest possible population"),
          shiny::p("Imagine a population of females and males producing a single cohort of offspring. Each offspring has one mother and one father, chosen independently at random. All females have the same chance of being a mother, and all males have the same chance of being a father. Actual family sizes vary by chance, even though expected reproductive output is equal within each sex. There are no age or size differences, and no adults die between reproduction and sampling."),
          shiny::p("Take random samples of females, males, and offspring, then identify every parent-offspring pair within those samples. Use those matches to estimate female and male abundance separately. With sample sizes held fixed, a match is more likely in a smaller population."),
          shiny::p("Choose population and sample sizes to create a trial, then repeat the experiment in independent replicates. Follow one replicate through reproduction and sampling, and use the trial-level figures to explore how chance and sampling effort affect pair counts and the precision of abundance estimates.")
        ),
        shiny::sidebarLayout(
          shiny::sidebarPanel(width = 3,
            shiny::h3("Prepare a trial"),
            shiny::p(class = "hint", "Each replicate creates a new population. These controls apply to the next trial."),
            shiny::h4("Population"),
            shiny::numericInput("Nf", "Females (Nf)", 1000, min = 1, max = limits$replicate_population, step = 100),
            shiny::numericInput("Nm", "Males (Nm)", 1000, min = 1, max = limits$replicate_population, step = 100),
            shiny::numericInput("Nk", "Offspring (Nk)", 2000, min = 1, max = limits$replicate_population, step = 100),
            shiny::h4("Sampling"),
            shiny::numericInput("Sf", "Sampled females (Sf)", 100, min = 1, step = 10),
            shiny::numericInput("Sm", "Sampled males (Sm)", 100, min = 1, step = 10),
            shiny::numericInput("Sk", "Sampled offspring (Sk)", 200, min = 1, step = 10),
            shiny::h4("Replication"),
            shiny::numericInput("num_reps", "Replicates", 100, min = 1, max = limits$reps, step = 10),
            shiny::numericInput("seed", "First replicate seed", 1001, min = 1, step = 1),
            shiny::p(class = "hint", "Replicate r uses seed + r \u2212 1. The next seed advances after a successful trial."),
            shiny::actionButton("run_trial", "Run and save trial", class = "btn-primary"),
            shiny::hr(),
            shiny::p(class = "hint", paste0("Complete samples and population summaries are kept in this session. Limits: ",
              limits$trials, " trials; ", limits$trial_samples / 1e6,
              " million sampled individuals across reps per trial; ", limits$session_samples / 1e6,
              " million sampled individuals per session; ", limits$replicate_population / 1e6,
              " million population individuals in any one replicate.")),
            shiny::downloadButton("download_trials", "Save complete trials (RDS)"),
            shiny::br(), shiny::br(),
            shiny::downloadButton("download_register", "Save trial table (CSV)")
          ),
          shiny::mainPanel(width = 9,
            shiny::div(class = "section",
              shiny::h3("Explore a saved trial"),
              shiny::p("Run a trial to begin. Select a saved trial and a replicate to inspect; changing the controls on the left does not change saved results."),
              shiny::selectInput("trial_id", "Saved trial", choices = character()),
              shiny::uiOutput("saved_settings"),
              shiny::uiOutput("rep_control"),
              shiny::textOutput("rep_details")
            ),
            shiny::div(class = "section",
              shiny::h3("01 \u00b7 Population / Reproduction"),
              shiny::p("Every adult has the same expected reproductive output. Each offspring independently receives one randomly chosen mother and father. Realized family sizes still vary."),
              shiny::plotOutput("population_plot", height = "290px")
            ),
            shiny::div(class = "section",
              shiny::h3("02 \u00b7 Sampling"),
              shiny::p("Females, males, and offspring are sampled independently without replacement. Each plot counts only sampled offspring belonging to sampled parents, including parents with no matches."),
              shiny::tableOutput("sampling_table"),
              shiny::plotOutput("sampling_plot", height = "290px")
            ),
            shiny::div(class = "section",
              shiny::h3("03 \u00b7 Pairs"),
              shiny::p("These distributions summarize all replicates in the saved trial. A parent with three sampled offspring contributes three pairs. Kin identification is assumed perfect."),
              shiny::plotOutput("pairs_plot", height = "310px")
            ),
            shiny::div(class = "section",
              shiny::h3("04 \u00b7 Estimation"),
              shiny::p("For each sex, E(#pairs) = sampled-adults \u00d7 sampled-offspring / num-adults. The moment-based estimator is num-adults = sampled-adults \u00d7 sampled-offspring / observed-number-of-pairs."),
              shiny::p("For each sex, abundance = sampled parents \u00d7 sampled offspring / observed pairs. The Poisson pseudo-likelihood uses that same expected pair count, so its continuous maximum gives the same point estimate."),
              shiny::p("Zero pairs give an unbounded estimate. Those replicates are counted in every panel, but cannot appear on a finite histogram axis. The likelihood is a working model."),
              shiny::plotOutput("estimates_plot", height = "460px"),
              shiny::h4("Approximate 95% confidence intervals"),
              shiny::p("Each thin line is one replicate, sorted from bottom to top by lower limit, then upper limit. The selected replicate has a thicker line and a point at its estimate. Both estimators give the same intervals here. Wald limits use the estimate plus or minus 1.96 standard errors; with few matches this approximation can be poor and lower limits may be negative. Zero matches have no Wald interval."),
              shiny::plotOutput("intervals_plot", height = "420px"),
              shiny::div(style = "overflow-x:auto;", shiny::tableOutput("estimation_table")),
              shiny::p(class = "hint", "Median, mean, and CV summarize finite estimates only. CV is the sample standard deviation divided by the mean, shown as a percentage; it is unavailable with fewer than two finite estimates. For Poisson pseudo-MLE rows, finite 95% CI inclusion reports how many finite-replicate confidence intervals contain the true abundance, divided by the number of finite replicates. A single replicate is an illustration, not a distribution of estimation performance.")
            ),
            shiny::div(class = "section",
              shiny::h3("Saved trials \u00b7 Equal reproduction"),
              shiny::p("This scenario's register preserves the settings for each trial. Scroll across for seeds and summaries; select a trial above to explore its replicates."),
              shiny::div(class = "trial-scroll", shiny::tableOutput("trial_register")),
              shiny::textOutput("storage_details")
            )
          )
        )
      )
    )
  )
}

#' Coordinate immutable trials and reactive views for scenario one
#'
#' The only trial mutation is appending a complete, successfully generated trial.
#' Input edits and view changes never modify an existing object. Complete samples
#' and population frequency summaries are retained in a session-local list. The
#' trial register is a derived tibble. Budgets distinguish retained sample rows
#' from the transient population generated for one replicate.
#' @param input,output,session Standard Shiny server arguments.
#' @return Used for its side effects on the Shiny session.
#' @keywords internal
pop_app_server <- function(input, output, session) {
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
      settings <- pop_settings(input$Nf, input$Nm, input$Nk, input$Sf, input$Sm, input$Sk)
      .pop_integer(input$seed, "Seed")
      .pop_check_budget(settings, input$num_reps, stored_individuals(), length(trials()), limits)
      id <- length(trials()) + 1L
      trial <- shiny::withProgress(message = paste("Simulating trial", id), value = 0, {
        pop_run_trial(settings, input$num_reps, input$seed, id,
                       progress = function(done, total) shiny::setProgress(value = done / total))
      })
      saved <- trials()
      saved[[as.character(id)]] <- trial
      trials(saved)
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
      shiny::numericInput("rep_id", "Replicate to inspect (one replicate saved)",
                           value = 1, min = 1, max = 1, step = 1)
    } else {
      shiny::sliderInput("rep_id", "Replicate to inspect", min = 1, max = n,
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
      paste0(x$num_reps, " independent population replicates \u00b7 seeds ", x$seed,
             "\u2013", as.double(x$seed) + x$num_reps - 1))
  })
  output$rep_details <- shiny::renderText({
    x <- current_rep()
    paste0("Replicate ", current_rep_id(), " | seed ", x$seed,
           " | mother-offspring pairs: ", sum(x$pairs$sex == "Female"),
           " | father-offspring pairs: ", sum(x$pairs$sex == "Male"))
  })
  output$population_plot <- shiny::renderPlot(pop_plot_reproduction(current_rep()))
  output$sampling_plot <- shiny::renderPlot(pop_plot_reproduction(current_rep(), sampled = TRUE))
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
    pop_trial_table(trials())
  }, striped = TRUE, digits = 2)
  output$storage_details <- shiny::renderText({
    paste(length(trials()), "of", limits$trials, "trial slots used |",
          format(stored_individuals(), big.mark = ",", scientific = FALSE),
          "of", format(limits$session_samples, big.mark = ",", scientific = FALSE),
          "sampled individuals retained across all replicates.")
  })
  output$download_trials <- shiny::downloadHandler(
    filename = function() "ckmr-equal-reproduction-trials.rds",
    content = function(file) saveRDS(trials(), file)
  )
  output$download_register <- shiny::downloadHandler(
    filename = function() "ckmr-equal-reproduction-trials.csv",
    content = function(file) utils::write.csv(pop_trial_table(trials()), file, row.names = FALSE)
  )
}
