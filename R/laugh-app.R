#' Build the curve editor and microscoping design-calculation interface
#' @return Embedded Shiny interface with manual calculation and saved snapshots.
#' @export
laugh_app_ui <- function() {
  info <- .laugh_curve_info()
  shiny::tagList(
    shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "laugh-editor.css"),
                     shiny::tags$script(src = "laugh-editor.js")),
    shiny::div(class = "hero", shiny::div(class = "eyebrow", "CKMR \u00b7 sampling design"),
      shiny::h1("Sample design---laugh test"),
      shiny::p("Shape age schedules, inspect the numbers, and see how a proposed sample is allocated."),
      shiny::div(class = "laugh-stage", shiny::strong("Explore, then calculate"),
        " \u00b7 Curves and allocation previews update immediately. Calculate expected kin pairs when your design is ready.")),
    shiny::div(class = "section",
      shiny::h3("What are we describing?"),
      shiny::p("One age-structured population, with a stable age distribution, equal numbers of females and males, and the same demographic, reproductive, and sampling schedules for both sexes. Population and sampling totals include both sexes. The laugh test is a rough feasibility calculation: it predicts expected kin-pair counts under these assumptions, rather than fitting an abundance estimator."),
      shiny::p("All schedules start at age 0. The final input age is a plus group, not a maximum lifespan. For example, 12+ means every animal aged 12 or older. Each curve below has one value for that entire final category."),
      shiny::p("The default curves are an illustrative example, not a stock assessment. For a real scoping exercise, population numbers should come from a defensible assessment or another appropriate source.")),
    shiny::div(class = "section",
      shiny::h3("01 \u00b7 Population age grid and abundance"),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput("laugh_plus_age", "Plus group begins at age", 12, min = 4, max = 60, step = 1),
          shiny::actionButton("laugh_apply_age", "Set age grid"),
          shiny::p(class = "hint", "This replaces the population curves and every active or cached sampling curve with fresh examples on the new grid. Undo restores the previous grid and curves. Merely typing an age does not change them.")),
        shiny::column(6, shiny::numericInput("laugh_abundance", "Set total population abundance to", 1e5, min = 1, step = 1000),
          shiny::actionButton("laugh_scale_abundance", "Set total abundance"),
          shiny::p(class = "hint", "Rescales the population curve while preserving its shape. Includes both sexes and the entire plus group. Drawing individual ages can change the total again."))),
      shiny::numericInput("laugh_juvenile", "Oldest age treated as juvenile for comparisons", 2, min = 0, max = 11, step = 1),
      shiny::p("This is MAX_JUVE_AGE. Ages 0 through this cutoff are included in juvenile-juvenile sibling and grandparent comparisons, and identify juvenile offspring in POP summaries. It is a comparison and age-precision cutoff, not a declaration that every individual below it is immature. The dashed line on each editor marks its boundary."),
      shiny::p(class = "hint", "Choose an exact age below the plus group. Older ages inside a plus group cannot be separated using one aggregate entry."),
      shiny::tags$details(shiny::tags$summary("More about the calculation"),
        shiny::p("The package excludes same-birth-cohort sibling comparisons. Grandparent-grandchild matches can resemble half-sibling matches genetically, so the calculation reports both. The juvenile cutoff helps limit that ambiguity and restricts attention to offspring that can be aged adequately."),
        shiny::numericInput("laugh_catch_min", "Youngest age used to fit the catch curve", 1, min = 0, max = 9, step = 1),
        shiny::p("MIN_CATCH_CURVE_AGE controls which population ages are used to smooth numbers at age and infer survival. It can exclude the youngest ages if their mortality differs strongly. This does not remove those ages from your population or samples. The input plus group is excluded from fitting. This editor leaves at least three exact-age classes for the fit."),
        shiny::checkboxInput("laugh_after_spawning", "Allow a parent sampled in its offspring's birth year", FALSE),
        shiny::p("This sets beqy_POP. The default FALSE excludes these same-year parent/offspring comparisons because parents sampled during spawning may not yet have completed reproduction. TRUE is appropriate only if adult sampling is understood to happen after spawning in the year convention being used."),
        shiny::p("The calculation keeps the supplied plus-group abundance (extraplus = FALSE) and expands it into older ages internally. Body weight and maturity stay at their last supplied values. The default expansion tolerance, etriv = 0.01, targets an old-age tail accounting for about 1% of fecundity. This processing runs only when you calculate; it does not change the editor curves."))),
    shiny::div(class = "laugh-toolbar",
      shiny::actionButton("laugh_undo", "Undo last curve edit"),
      shiny::actionButton("laugh_reset", "Reset all curves"),
      shiny::tags$a(href="#laugh_results_section",class="btn btn-primary","Expected kin counts"),
      shiny::span(shiny::textOutput("laugh_edit_status", inline = TRUE))),
    lapply(c("nata","wata","pmatata"), function(key) {
      x <- info[[key]]
      shiny::div(class = "section laugh-curve-section",
        shiny::h3(x$title), shiny::p(class = "laugh-argument", paste0("Input: ", key)),
        shiny::p(x$definition), shiny::p(x$consequence),
        shiny::div(id = paste0("laugh_editor_", key), class = "laugh-editor", `data-key` = key),
        shiny::tags$details(shiny::tags$summary("How the package uses this schedule"), shiny::p(x$package)))
    }),
    shiny::div(class = "section",
      shiny::h3("02 \u00b7 Sampling years, numbers, and age composition"),
      shiny::textInput("laugh_years", "Samples in consecutive years", "500, 500, 500"),
      shiny::p("Enter total samples of both sexes combined, one count per year, separated by commas or spaces. The model assumes half are female and half male. Add another number at the end to add a year: 500, 500, 500, 800 gives four consecutive years. Use zero to retain a gap year; 500, 0, 1000 means sampling in years 1 and 3. Up to 12 years can be previewed."),
      shiny::p(class="hint","Year numbers refer to positions in this list. Removing entries from the end hides those years' curves without deleting them. Re-adding those years restores their previous curves; genuinely new years copy the selected year's curve. Changing counts alone leaves curves unchanged."),
      shiny::fluidRow(shiny::column(6,shiny::selectInput("laugh_sampling_year","Year to edit",
        choices=c("Year 1"="1","Year 2"="2","Year 3"="3"),selected="1",selectize=FALSE)),
        shiny::column(6,shiny::actionButton("laugh_copy_years","Copy this curve to all active years"))),
      shiny::textOutput("laugh_year_cache_status"), shiny::uiOutput("laugh_year_legend"),
      shiny::p(info$samp_pata$definition), shiny::p(info$samp_pata$consequence),
      shiny::p(class="hint","Only the selected year's curve changes when you draw, edit numbers, paste, apply a shape, or multiply values. Other active years appear as muted lines behind it. Copy-to-all affects the current design's years; hidden curves are retained unchanged. All curve changes can be undone."),
      shiny::div(id="laugh_editor_samp_pata",class="laugh-editor",`data-key`="samp_pata"),
      shiny::tags$details(shiny::tags$summary("How the sampling matrix will be used"),shiny::p(info$samp_pata$package)),
      shiny::h4("Normalized proportions and allocated sample numbers"),
      shiny::p("Each year's scores are divided by that row's sum. For example, a score of 2 among scores summing to 10 allocates 20% of that year's samples to the age category. The matrix has one row per active year and one column per age category. Even a zero-sample year retains a normalized row; multiplying it by zero gives zero allocated samples."),
      shiny::textOutput("laugh_preview_status"),
      shiny::plotOutput("laugh_sample_mix", height = "260px"),
      shiny::plotOutput("laugh_allocation_plot", height = "300px"),
      shiny::p(class = "hint", "Allocations can be fractional: these are expected numbers allocated to age categories, not a random or rounded realization. The plus-group allocation includes every age in that group; it is not a per-age allocation within the group."),
      shiny::div(class = "trial-scroll", shiny::tableOutput("laugh_allocations"))),
    shiny::div(class = "section",
      shiny::h3("03 \u00b7 Inspect and export the inputs"),
      shiny::textOutput("laugh_readiness"),
      shiny::p(class = "hint", "Rows below are grouped by sampling year. relative_output is biological weight times proportion mature, and has not been normalized to mean one. sampling_score is that year's drawn value; samp_pata is its row-normalized sample proportion. Hidden years are excluded from previews and exports."),
      shiny::div(class = "trial-scroll", shiny::tableOutput("laugh_input_table")),
      shiny::uiOutput("laugh_export_controls"),
      shiny::p("The CSV includes year-specific scores, normalized proportions, and expected allocations. The R file contains the actual package arguments: samp_nya is the year-by-age matrix of allocated sample NUMBERS. It is the proportion matrix multiplied by annual sample counts, not a matrix of proportions. Its optional example call is commented out."),
      shiny::tags$details(shiny::tags$summary("Show the actual R argument list passed to microscoping"), shiny::verbatimTextOutput("laugh_exact_inputs"))),
    shiny::div(id="laugh_results_section",class="section",
      shiny::h3("04 \u00b7 Expected kin counts"),
      shiny::p("Run the current design through microscoping::ckmr_laugh_test(). Each successful calculation saves its exact inputs and package version. Editing curves afterward does not change saved calculations; choose an earlier one to compare results."),
      shiny::textOutput("laugh_engine_status"),
      shiny::actionButton("laugh_calculate","Calculate expected kin pairs",class="btn-primary"),
      shiny::textOutput("laugh_calculation_error"),
      shiny::selectInput("laugh_result_id","Saved calculation",choices=character(),selectize=FALSE),
      shiny::uiOutput("laugh_result_status"),
      shiny::plotOutput("laugh_result_plot",height="340px"),
      shiny::div(style="overflow-x:auto;",shiny::tableOutput("laugh_result_table")),
      shiny::textOutput("laugh_ggp_ratio"),
      shiny::p(class="hint","Counts include both sexes and all sampling years. HSPs and GGPs are restricted to juvenile-juvenile comparisons and exclude pairs born in the same cohort. GGPs can resemble HSPs genetically; their ratio helps assess that ambiguity. These are rough expected counts, not estimates of precision or proof that a full CKMR model will be identifiable."),
      shiny::uiOutput("laugh_result_downloads"),
      shiny::tags$details(shiny::tags$summary("Population smoothing and plus-group expansion for this saved calculation"),
        shiny::textOutput("laugh_preparation_status"),
        shiny::plotOutput("laugh_preparation_plot",height="300px"),
        shiny::p(class="hint","Points show supplied exact-age numbers; the line shows the smoothed population with the plus group spread over older ages. The aggregate plus-group input is reported above rather than plotted as though it were one exact age.")),
      shiny::tags$details(shiny::tags$summary("Exact inputs for the selected saved calculation"),
        shiny::verbatimTextOutput("laugh_saved_inputs"))),
    shiny::div(class = "section", shiny::h3("Things to try"), shiny::tags$ol(
      shiny::tags$li("Draw a sample concentrated on young animals, then add some older ages. Watch both the age proportions and the year-by-age allocations."),
      shiny::tags$li("Multiply all sampling scores by two. Why does the normalized sample composition stay the same? Now double an annual sample count instead."),
      shiny::tags$li("Change body weights or the maturity schedule. Inspect their product in the input table: it is the relative reproductive output at age assumed by the laugh test."),
      shiny::tags$li("Try a zero-sample gap year. Notice that it remains a year in the allocation table.")))
  )
}

#' Bind the input editor and manually triggered microscoping calculations
#' @param input,output,session Standard Shiny server arguments.
#' @return Reactive editor state and prepared-input handles, for integration checks.
#' @export
laugh_app_server <- function(input, output, session) {
  state <- shiny::reactiveVal(laugh_default_curves())
  history <- shiny::reactiveVal(list())
  status <- shiny::reactiveVal("Draw across ages, or open a curve's numeric table. Edits stay here when you switch apps.")
  revision <- shiny::reactiveVal(0L)
  info <- .laugh_curve_info()
  selected_year <- shiny::reactiveVal(1L)
  active_counts <- shiny::reactiveVal(c(500,500,500))
  calculations <- shiny::reactiveVal(list())
  calculation_error <- shiny::reactiveVal("")
  commit <- function(next_state, message) {
    previous <- state()
    if (!identical(previous, next_state)) {
      history(tail(c(history(), list(previous)), 25))
      state(next_state); revision(revision()+1L)
    }
    status(message)
  }
  shiny::observeEvent(input$laugh_years,{
    counts <- tryCatch(laugh_parse_years(input$laugh_years),error=identity)
    if(inherits(counts,"error")) { status(conditionMessage(counts)); return() }
    s <- state(); previous <- length(active_counts()); cached <- nrow(s$sample_scores)
    n <- length(counts); template <- selected_year()
    next_state <- laugh_ensure_years(s,n,template)
    state(next_state); active_counts(counts)
    selected_year(min(template,n))
    revision(revision()+1L)
    labels <- paste0("Year ",seq_len(n)," \u00b7 ",format(counts,trim=TRUE,scientific=FALSE)," samples")
    shiny::updateSelectInput(session,"laugh_sampling_year",choices=stats::setNames(as.character(seq_len(n)),labels),
      selected=as.character(selected_year()))
    status(if(n>previous && n<=cached) paste0("Previous curves restored for year",if(n-previous>1)"s " else " ",
      paste((previous+1):n,collapse=", "),".") else if(n>cached)
        paste("New years copy year",template,"; any previously hidden curves have been restored.") else
          "Annual counts updated. Existing age curves are unchanged; removed years stay cached in this session.")
  },priority=50)
  shiny::observeEvent(input$laugh_sampling_year,{
    y <- suppressWarnings(as.integer(input$laugh_sampling_year))
    if(length(y)==1 && !is.na(y) && y>=1 && y<=length(active_counts()) && y!=selected_year()) {
      selected_year(y); revision(revision()+1L)
      status(paste("Editing year",y,"only. Other years' curves are unchanged."))
    }
  })
  edit_year <- function(payload) {
    y <- if(is.null(payload$year)) selected_year() else as.integer(payload$year)
    if(length(y)!=1 || is.na(y) || y<1 || y>length(active_counts())) stop("Choose an active sampling year.")
    y
  }
  # Server messages carry the full canonical state; a completed stroke sends one
  # edit, not hundreds of mouse-move events. Painting never calls the package.
  send <- function() {
    s <- laugh_ensure_years(state(),length(active_counts()),selected_year())
    juvenile <- input$laugh_juvenile
    if (is.null(juvenile) || !is.finite(juvenile)) juvenile <- 2
    session$sendCustomMessage("laugh_editor_state", list(plus_age = s$plus_age,
      juvenile = min(juvenile, s$plus_age-1), revision = revision(),
      curves = lapply(names(info), function(key) {
        values <- if(key=="samp_pata") unname(s$sample_scores[selected_year(),]) else s$curves[[key]]
        metadata <- if(key=="samp_pata") list(year=selected_year(),series=lapply(seq_along(active_counts()),function(y)
          list(year=y,label=paste("Year",y),color=.laugh_year_colors()[y],values=unname(s$sample_scores[y,])))) else list()
        c(list(key=key,values=values),metadata,info[[key]])
      })))
  }
  shiny::observe({ state(); revision(); selected_year(); active_counts(); input$laugh_juvenile; send() })
  shiny::observeEvent(input$laugh_editor_ready, send(), ignoreInit = FALSE)
  shiny::observeEvent(input$laugh_curve_edit, {
    edit <- input$laugh_curve_edit
    tryCatch({
      if (is.null(edit$revision) || edit$revision != revision()) { send(); return() }
      s <- state(); key <- edit$key
      values <- laugh_validate_curve(s, key, as.numeric(unlist(edit$values)))
      if(key=="samp_pata") s$sample_scores[edit_year(edit),] <- values else s$curves[[key]] <- values
      commit(s, paste(info[[key]]$title,if(key=="samp_pata")paste("for year",edit_year(edit)) else "",
        "updated. Sampling previews use the current values."))
    }, error = function(e) { status(conditionMessage(e)); send() })
  })
  shiny::observeEvent(input$laugh_curve_action, {
    action <- input$laugh_curve_action
    tryCatch({
      if (is.null(action$revision) || action$revision != revision()) { send(); return() }
      s <- state(); key <- action$key
      if (!key %in% names(info)) stop("Unknown curve.")
      y <- if(key=="samp_pata") edit_year(action) else 1
      current <- if(key=="samp_pata") s$sample_scores[y,] else s$curves[[key]]
      values <- if (action$action == "preset") laugh_curve_preset(s,key,action$preset,year=y) else {
        factor <- as.numeric(action$factor)
        if (!is.finite(factor) || !factor %in% c(.5, 2)) stop("Unknown scale action.")
        result <- current * factor
        if (key == "pmatata") result <- pmin(result, 1)
        result
      }
      values <- laugh_validate_curve(s,key,values)
      if(key=="samp_pata") s$sample_scores[y,] <- values else s$curves[[key]] <- values
      commit(s,paste(info[[key]]$title,if(key=="samp_pata")paste("for year",y) else "","updated."))
    }, error = function(e) { status(conditionMessage(e)); send() })
  })
  shiny::observeEvent(input$laugh_copy_years,{
    s <- state(); n <- length(active_counts()); y <- selected_year()
    s$sample_scores[seq_len(n),] <- matrix(rep(s$sample_scores[y,],n),nrow=n,byrow=TRUE)
    commit(s,paste("Year",y,"copied to all",n,"active years. Hidden curves are unchanged; Undo is available."))
  })
  shiny::observeEvent(input$laugh_apply_age, {
    tryCatch({
      next_state <- laugh_default_curves(input$laugh_plus_age,max(1,sum(state()$curves$nata)),nrow(state()$sample_scores))
      commit(next_state, "New age grid and example curves, including cached years. Undo restores the previous grid and curves.")
    }, error = function(e) status(conditionMessage(e)))
  })
  shiny::observeEvent(input$laugh_scale_abundance, {
    tryCatch({
      value <- input$laugh_abundance
      if (length(value) != 1L || !is.finite(value) || value <= 0) stop("Enter a positive total abundance.")
      s <- state(); total <- sum(s$curves$nata)
      if (total <= 0) stop("Choose or draw a nonzero population shape first.")
      s$curves$nata <- s$curves$nata * value / total
      commit(s, "Population shape rescaled to the requested total abundance.")
    }, error = function(e) status(conditionMessage(e)))
  })
  shiny::observeEvent(input$laugh_reset, {
    commit(laugh_default_curves(state()$plus_age,num_years=nrow(state()$sample_scores)),
      "All curves, including cached years, reset to the illustrative example. Undo is available.")
  })
  shiny::observeEvent(input$laugh_undo, {
    h <- history()
    if (!length(h)) { status("No earlier curve edit to undo."); return() }
    restored <- h[[length(h)]]; current <- state()
    # A year added since this snapshot must not disappear when an earlier curve
    # edit is undone. Preserve its cached row when the age grid is unchanged.
    old_rows <- nrow(restored$sample_scores); new_rows <- nrow(current$sample_scores)
    if(restored$plus_age==current$plus_age && new_rows>old_rows)
      restored$sample_scores <- rbind(restored$sample_scores,current$sample_scores[(old_rows+1):new_rows,,drop=FALSE])
    restored <- laugh_ensure_years(restored,length(active_counts()),min(selected_year(),nrow(restored$sample_scores)))
    state(restored); history(head(h, -1)); revision(revision()+1L)
    shiny::updateNumericInput(session, "laugh_plus_age", value = state()$plus_age)
    status("Previous age grid and curves restored. Annual sample counts and comparison settings were not changed.")
  })
  shiny::observe({
    a <- state()$plus_age
    shiny::updateNumericInput(session, "laugh_juvenile", max = a-1,
      value = if (is.null(input$laugh_juvenile)) 2 else min(input$laugh_juvenile, a-1))
    shiny::updateNumericInput(session, "laugh_catch_min", max = a-3,
      value = if (is.null(input$laugh_catch_min)) 1 else min(input$laugh_catch_min, a-3))
  })
  years <- shiny::reactive(laugh_parse_years(input$laugh_years))
  prepared <- shiny::reactive(tryCatch(laugh_prepare_inputs(state(), years(),
    juvenile = input$laugh_juvenile, catch_min = input$laugh_catch_min,
    after_spawning = isTRUE(input$laugh_after_spawning)), error = identity))
  arguments <- shiny::reactive({
    x <- prepared()
    shiny::validate(shiny::need(!inherits(x, "error"), if (inherits(x,"error")) conditionMessage(x) else ""))
    x
  })
  age_labels <- shiny::reactive(c(as.character(0:(state()$plus_age-1)), paste0(state()$plus_age,"+")))
  allocation <- shiny::reactive(laugh_sample_allocation(arguments()$samp_pata, arguments()$samp_ny))
  output$laugh_edit_status <- shiny::renderText(status())
  output$laugh_preview_status <- shiny::renderText({
    x <- arguments()
    paste("Current grid: 0 to", paste0(state()$plus_age,"+."),
      "Total population:", format(sum(x$nata),big.mark=",",trim=TRUE,scientific=FALSE),
      "|", length(x$samp_ny), "consecutive sampling years |", format(sum(x$samp_ny),big.mark=",",scientific=FALSE), "total samples.")
  })
  output$laugh_sample_mix <- shiny::renderPlot({
    x <- arguments()
    y <- selected_year()
    frame <- data.frame(age = factor(age_labels(), levels = age_labels()), proportion = x$samp_pata[y,])
    ggplot2::ggplot(frame,ggplot2::aes(age,proportion)) + ggplot2::geom_col(fill=.laugh_year_colors()[y],width=.7) +
      ggplot2::scale_y_continuous(labels=function(x) paste0(round(100*x),"%")) +
      ggplot2::labs(x="Age category",y="Share of sample",title=paste("Year",y,"normalized age composition"),
        caption="Each year's sampling scores are normalized separately. The heatmap below shows all active years.") + .pop_plot_theme()
  })
  output$laugh_allocation_plot <- shiny::renderPlot({
    x <- allocation()
    x$age <- factor(age_labels()[x$age+1],levels=age_labels())
    ggplot2::ggplot(x,ggplot2::aes(age,factor(year),fill=expected_samples)) +
      ggplot2::geom_tile(color="white") + ggplot2::scale_fill_gradient(low="#edf3f8",high="#3c6ea8",name="Allocated\nsamples") +
      ggplot2::labs(x="Age category",y="Sampling year",caption="Zero-sample years stay in the timeline. Values are expected allocations, not simulated counts.") +
      .pop_plot_theme() + ggplot2::theme(legend.position="right")
  })
  output$laugh_allocations <- shiny::renderTable({
    x <- arguments()
    result <- data.frame(Age=age_labels(),check.names=FALSE)
    for (i in seq_along(x$samp_ny)) result[[paste("Year",i)]] <- round(x$samp_pata[i,]*x$samp_ny[i],2)
    result
  },striped=TRUE,rownames=FALSE)
  output$laugh_readiness <- shiny::renderText({
    x <- prepared()
    if (inherits(x,"error")) paste("Inputs need attention:",conditionMessage(x)) else
      "Inputs are ready. Allocation previews update automatically; click Calculate expected kin pairs to run and save this design."
  })
  table_data <- shiny::reactive({
    s <- state(); c <- s$curves; x <- arguments(); n <- length(x$samp_ny); a <- s$plus_age+1
    data.frame(year=rep(seq_len(n),each=a),age=rep(age_labels(),n),
      nata=rep(c$nata,n),wata=rep(c$wata,n),pmatata=rep(c$pmatata,n),
      relative_output=rep(c$wata*c$pmatata,n),sampling_score=as.vector(t(s$sample_scores[seq_len(n),,drop=FALSE])),
      samp_pata=as.vector(t(x$samp_pata)),annual_samples=rep(x$samp_ny,each=a),
      allocated_samples=as.vector(t(x$samp_pata*x$samp_ny)))
  })
  output$laugh_input_table <- shiny::renderTable(table_data(),striped=TRUE,digits=5,rownames=FALSE)
  output$laugh_exact_inputs <- shiny::renderPrint(laugh_function_arguments(arguments()))
  output$laugh_export_controls <- shiny::renderUI({
    if (inherits(prepared(),"error")) shiny::p(class="hint","Complete the inputs above to enable export.") else
      shiny::tagList(shiny::downloadButton("laugh_download_csv","Download age schedules (CSV)"),
                    " ",shiny::downloadButton("laugh_download_r","Download R inputs"))
  })
  output$laugh_download_csv <- shiny::downloadHandler(filename=function()"ckmr-laugh-test-age-inputs.csv",
    content=function(file) { arguments(); utils::write.csv(table_data(),file,row.names=FALSE) })
  output$laugh_download_r <- shiny::downloadHandler(filename=function()"ckmr-laugh-test-inputs.R",
    content=function(file) writeLines(c("# Actual package arguments for the current editor design.",
      "# samp_nya contains sample NUMBERS, not proportions. First column is age 0.",
      "inputs <-",capture.output(dput(laugh_function_arguments(arguments()))),"",
      "# With microscoping installed and samp_nya support available:",
      "# result <- do.call(microscoping::ckmr_laugh_test, inputs)"),file))
  shiny::observeEvent(input$laugh_calculate,{
    tryCatch({
      design <- arguments()
      editor_snapshot <- state()
      result <- shiny::withProgress(message="Calculating expected kin pairs",value=.1,{
        computed <- laugh_run_design(design)
        shiny::setProgress(value=1)
        computed
      })
      saved <- calculations(); id <- length(saved)+1L
      result$calculation_id <- id
      result$editor_state <- editor_snapshot
      saved[[as.character(id)]] <- result
      calculations(saved); calculation_error("")
      labels <- vapply(saved,function(x) paste0("Calculation ",x$calculation_id," | ",
        length(x$design$samp_ny)," years | ",format(sum(x$design$samp_ny),big.mark=",",scientific=FALSE,trim=TRUE)," samples"),character(1))
      shiny::updateSelectInput(session,"laugh_result_id",choices=stats::setNames(names(saved),labels),selected=as.character(id))
    },error=function(e) calculation_error(paste("Calculation not completed:",conditionMessage(e),
      "Previous saved calculations are unchanged.")))
  })
  selected_calculation <- shiny::reactive({
    saved <- calculations(); shiny::req(length(saved)>0)
    id <- input$laugh_result_id
    if(length(id)!=1L || is.na(id) || !id %in% names(saved)) id <- tail(names(saved),1)
    saved[[id]]
  })
  output$laugh_engine_status <- shiny::renderText(tryCatch(
    paste("Ready: installed microscoping",laugh_engine_version(),"with year-by-age sample-number support."),
    error=function(e) conditionMessage(e)))
  output$laugh_calculation_error <- shiny::renderText(calculation_error())
  output$laugh_result_status <- shiny::renderUI({
    r <- selected_calculation(); current <- prepared()
    same <- tryCatch(!inherits(current,"error") && identical(laugh_function_arguments(current),r$function_arguments),error=function(e)FALSE)
    shiny::div(class="saved-settings",shiny::strong(paste("Saved calculation",r$calculation_id)),
      shiny::p(paste("microscoping",r$package_version,"|",length(r$design$samp_ny),"sampling years | plus group",length(r$design$nata)-1)),
      shiny::p(if(same) "These saved results match the current inputs." else
        "The editor inputs differ from this saved calculation. Calculate again to save updated results; the results below still show the selected saved design."),
      if(length(r$warnings)) shiny::tagList(shiny::strong("Package/preparation warnings:"),
        do.call(shiny::tags$ul,lapply(r$warnings,shiny::tags$li))))
  })
  output$laugh_result_plot <- shiny::renderPlot(laugh_plot_results(selected_calculation()))
  output$laugh_result_table <- shiny::renderTable({
    x <- laugh_result_table(selected_calculation())
    data.frame(Output=x$label,Expected=vapply(x$expected,function(v)
      format(signif(v,6),scientific=FALSE,trim=TRUE),character(1)),Meaning=x$meaning,check.names=FALSE)
  },striped=TRUE,rownames=FALSE)
  output$laugh_ggp_ratio <- shiny::renderText({
    x <- selected_calculation()$expected
    if(x[["HSP_JJ"]]>0) paste("Expected GGP/HSP ratio:",
      paste0(format(signif(100*x[["GGP_JJ"]]/x[["HSP_JJ"]],4),scientific=FALSE,trim=TRUE),"%.")) else
        "GGP/HSP ratio is undefined because no juvenile HSPs are predicted."
  })
  output$laugh_preparation_status <- shiny::renderText({
    p <- selected_calculation()$preparation
    paste0("The supplied ",p$input_plus_age,"+ group contains ",
      format(signif(p$input_plus_n,6),scientific=FALSE,trim=TRUE)," animals. Its inferred annual survival is ",
      signif(p$plus_survival,5),". The expanded calculation runs from age 0 to ",max(p$ages),
      ". The edited population curve has not been overwritten.")
  })
  output$laugh_preparation_plot <- shiny::renderPlot({
    r <- selected_calculation(); p <- r$preparation; a <- p$input_plus_age
    model <- data.frame(age=p$ages,n=p$expanded_nata)
    raw <- data.frame(age=0:(a-1L),n=r$design$nata[seq_len(a)])
    ggplot2::ggplot(model,ggplot2::aes(age,n)) + ggplot2::geom_line(color="#237e83",linewidth=.9) +
      ggplot2::geom_point(data=raw,color="#8065a6",size=2) +
      ggplot2::geom_vline(xintercept=a-.5,linetype="dashed",color="grey50") +
      ggplot2::labs(x="Age",y="Population numbers (both sexes)",
        caption="Purple points: supplied exact-age numbers. Teal line: smoothed and expanded population.\nDashed: start of the input plus group; the plus-group aggregate is not plotted as an exact-age point.") + .pop_plot_theme()
  })
  output$laugh_saved_inputs <- shiny::renderPrint(selected_calculation()$function_arguments)
  output$laugh_result_downloads <- shiny::renderUI({
    selected_calculation()
    shiny::tagList(shiny::downloadButton("laugh_download_result_csv","Save expected counts (CSV)")," ",
      shiny::downloadButton("laugh_download_result_rds","Save calculation and inputs (RDS)"))
  })
  output$laugh_download_result_csv <- shiny::downloadHandler(filename=function()"ckmr-laugh-test-expected-counts.csv",
    content=function(file) {
      r <- selected_calculation(); x <- laugh_result_table(r)
      x$calculation_id <- r$calculation_id; x$microscoping_version <- r$package_version
      utils::write.csv(x,file,row.names=FALSE)
    })
  output$laugh_download_result_rds <- shiny::downloadHandler(filename=function()"ckmr-laugh-test-calculation.rds",
    content=function(file) saveRDS(selected_calculation(),file))
  output$laugh_year_cache_status <- shiny::renderText({
    hidden <- nrow(state()$sample_scores)-length(active_counts())
    paste("Editing year",selected_year(),"|",length(active_counts()),if(length(active_counts())==1)"active year |" else "active years |",hidden,
      "hidden curves retained in this session. Previews and exports include active years only.")
  })
  output$laugh_year_legend <- shiny::renderUI(shiny::div(class="laugh-year-legend",
    lapply(seq_along(active_counts()),function(y) shiny::span(
      class=if(y==selected_year())"laugh-year-active" else "laugh-year-muted",
      style=paste0("border-color:",.laugh_year_colors()[y],";"),
      shiny::span(style=paste0("color:",.laugh_year_colors()[y],";"),"\u25cf "),paste("Year",y),
      if(y==selected_year())" (editing)"))))
  invisible(list(state=state,history=history,arguments=arguments,prepared=prepared,allocation=allocation,
    revision=revision,selected_year=selected_year,active_counts=active_counts,
    calculations=calculations,selected_calculation=selected_calculation,calculation_error=calculation_error))
}
