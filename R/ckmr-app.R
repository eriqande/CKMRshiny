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
      shiny::tabPanel("Sample design---laugh test", value = "laugh", laugh_app_ui()))
  )
}

#' Bind the app chooser and persistent simulation/input-editor servers
#' @param input,output,session Standard Shiny server arguments.
#' @return Handles for the two simulations, design inputs, and session budgets.
#' @export
ckmr_app_server <- function(input, output, session) {
  usage <- .pop_session_usage()
  session$userData$ckmr_usage <- usage
  pop <- pop_app_server(input, output, session)
  hsp <- hsp_app_server(input, output, session, shared_usage = usage)
  laugh <- laugh_app_server(input, output, session)
  shiny::observeEvent(input$choose_app, shiny::updateTabsetPanel(session, "app_view", selected = "home"))
  shiny::observeEvent(input$open_pop, shiny::updateTabsetPanel(session, "app_view", selected = "pop"))
  shiny::observeEvent(input$start_pop, shiny::updateTabsetPanel(session, "app_view", selected = "pop"))
  shiny::observeEvent(input$open_hsp, shiny::updateTabsetPanel(session, "app_view", selected = "hsp"))
  shiny::observeEvent(input$start_hsp, shiny::updateTabsetPanel(session, "app_view", selected = "hsp"))
  shiny::observeEvent(input$open_laugh, shiny::updateTabsetPanel(session, "app_view", selected = "laugh"))
  shiny::observeEvent(input$start_laugh, shiny::updateTabsetPanel(session, "app_view", selected = "laugh"))
  invisible(list(pop = pop, hsp = hsp, laugh = laugh, usage = usage))
}
