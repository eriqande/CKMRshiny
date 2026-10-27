# Run from this directory with shiny::runApp(). The exploratory file
# R/pop-fundamentals.R is deliberately not sourced.
for (file in c("pop-simulation.R", "pop-plots.R", "pop-app.R")) {
  source(file.path("R", file), local = TRUE, encoding = "UTF-8")
}
shiny::shinyApp(ui = pop_app_ui(), server = pop_app_server)
