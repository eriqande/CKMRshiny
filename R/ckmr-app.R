#' Build the POP/HSP/design chooser while retaining all apps in the session
#' @return Shiny page with static hidden tabs. Switching visibility never
#' destroys a scenario server, input, or saved-trial registry.
#' @export
ckmr_app_ui <- function() {
  shiny::fluidPage(
    shiny::div(class = "app-switcher", style = "padding:14px 8px; border-bottom:1px solid #dce3df; display:flex; align-items:center; gap:10px; flex-wrap:wrap;",
      shiny::strong("CKMR fundamentals"),
      shiny::actionButton("choose_app", "Switch app"),
      shiny::actionButton("open_pop", "Parent-offspring (POP)"),
      shiny::actionButton("open_hsp", "Siblings (HSP)"),
      shiny::actionButton("open_laugh", "Sample design---laugh test"),
      shiny::span(class = "hint", "Trials and input edits stay here when you switch.")),
    shiny::tabsetPanel(id = "app_view", type = "hidden", selected = "home",
      shiny::tabPanel("Choose an app", value = "home",
        shiny::div(class = "hero", shiny::div(class = "eyebrow", "CKMR \u00b7 population fundamentals"),
          shiny::h1("Explore what close kin tell us about abundance"),
          shiny::p("Choose a kin-pair experiment, or sketch age schedules and a sampling design.")),
        shiny::fluidRow(
          shiny::column(4, shiny::div(class = "section", shiny::h3("Parent-offspring pairs"),
            shiny::p("Sample adults and offspring. Explore equal reproduction, random reproductive differences, size-based fecundity, and selective adult sampling."),
            shiny::actionButton("start_pop", "Explore parent-offspring pairs", class = "btn-primary"))),
          shiny::column(4, shiny::div(class = "section", shiny::h3("Sibling pairs"),
            shiny::p("Sample offspring. Compare within-cohort and cross-cohort matches, seasonal luck, and persistent size effects."),
            shiny::actionButton("start_hsp", "Explore sibling pairs", class = "btn-primary"))),
          shiny::column(4, shiny::div(class = "section", shiny::h3("Sample design---laugh test"),
            shiny::p("Paint age schedules, set sampling effort, and calculate expected kin counts with Mark's microscoping package."),
            shiny::actionButton("start_laugh", "Sample design---laugh test", class = "btn-primary"))))),
      shiny::tabPanel("POP", value = "pop", pop_app_ui(embedded = TRUE)),
      shiny::tabPanel("HSP", value = "hsp", hsp_app_ui()),
      shiny::tabPanel("Sample design---laugh test", value = "laugh", design_app_ui()))
  )
}

#' Bind the app chooser and persistent simulation/input-editor servers
#' @param input,output,session Standard Shiny server arguments.
#' @param design_inputs Optional microscoping argument list to import at startup.
#' @return Handles for the two simulations, both design editors, and session budgets.
#' @export
ckmr_app_server <- function(input, output, session, design_inputs = NULL) {
  usage <- .pop_session_usage()
  session$userData$ckmr_usage <- usage
  pop <- pop_app_server(input, output, session)
  hsp <- hsp_app_server(input, output, session, shared_usage = usage)
  laugh <- laugh_app_server(input, output, session)
  clt2 <- clt2_app_server(input, output, session)
  import_status <- shiny::reactiveVal("Choose an RDS argument list, then click Import inputs. Existing saved calculations are retained.")
  apply_design <- function(x) {
    if(x$kind=="laugh") laugh$import_design(x) else clt2$import_design(x)
    shiny::updateTabsetPanel(session,"app_view",selected="laugh")
    shiny::updateTabsetPanel(session,"design_view",selected=x$kind)
    import_status(paste(c(paste("Loaded",if(x$kind=="laugh") "one-sex laugh test" else "clt2 two-sex design","inputs. Undo is available on that tab."),x$notices),collapse=" "))
  }
  shiny::observeEvent(input$design_import,{
    tryCatch({shiny::req(input$design_file);x <- ckmr_import_design(readRDS(input$design_file$datapath));apply_design(x)},
      error=function(e) import_status(paste("Import not applied:",conditionMessage(e),"The editors and saved calculations are unchanged.")))
  })
  output$design_import_status <- shiny::renderText(import_status())
  if(!is.null(design_inputs)) {
    initial <- ckmr_import_design(design_inputs)
    session$onFlushed(function() shiny::isolate(apply_design(initial)),once=TRUE)
  }
  shiny::observeEvent(input$choose_app, shiny::updateTabsetPanel(session, "app_view", selected = "home"))
  shiny::observeEvent(input$open_pop, shiny::updateTabsetPanel(session, "app_view", selected = "pop"))
  shiny::observeEvent(input$start_pop, shiny::updateTabsetPanel(session, "app_view", selected = "pop"))
  shiny::observeEvent(input$open_hsp, shiny::updateTabsetPanel(session, "app_view", selected = "hsp"))
  shiny::observeEvent(input$start_hsp, shiny::updateTabsetPanel(session, "app_view", selected = "hsp"))
  shiny::observeEvent(input$open_laugh, shiny::updateTabsetPanel(session, "app_view", selected = "laugh"))
  shiny::observeEvent(input$start_laugh, shiny::updateTabsetPanel(session, "app_view", selected = "laugh"))
  invisible(list(pop = pop, hsp = hsp, laugh = laugh, clt2 = clt2, import_status=import_status, usage = usage))
}

#' Build the sample-design tabs and shared data-only import panel
#' @return UI retaining both editors and their saved calculations in one session.
#' @export
design_app_ui <- function() shiny::tagList(
  shiny::div(class="section",shiny::h3("Start from an existing design"),
    shiny::p("Import a named list of arguments for ckmr_laugh_test() or clt2(), saved with saveRDS(inputs, 'design.rds'). The matching tab opens automatically. Validation happens before any inputs are replaced. No calculation is run on import."),
    shiny::fileInput("design_file","Argument list (RDS)",accept=c(".rds",".RDS")),
    shiny::actionButton("design_import","Import inputs"),shiny::textOutput("design_import_status"),
    shiny::p(class="hint","Supported editor grids: plus groups from 2 to 200 and 1-12 sampling years. Full counts take precedence over sampling margins. Unsupported arguments are reported. The one-sex editor retains the supplied plus group (extraplus = FALSE).")),
  shiny::tabsetPanel(id="design_view",selected="laugh",
    shiny::tabPanel("One-sex laugh test",value="laugh",laugh_app_ui()),
    shiny::tabPanel("clt2: a two-sex laugh test",value="clt2",clt2_app_ui())))

#' Launchable CKMR app initialized from a microscoping argument list
#' @param inputs Optional named argument list for ckmr_laugh_test() or clt2().
#' @return Shiny app object for shiny::runApp(). Inputs are validated immediately;
#' the matching editor is populated on session startup without calculating kin.
#' @export
ckmr_design_app <- function(inputs=NULL) {
  if(!is.null(inputs)) ckmr_import_design(inputs)
  shiny::shinyApp(ui=ckmr_app_ui(),server=function(input,output,session)
    ckmr_app_server(input,output,session,design_inputs=inputs))
}
