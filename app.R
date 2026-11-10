# Run from this directory with shiny::runApp(). The exploratory file
# R/pop-fundamentals.R is deliberately not sourced.
for (file in c("pop-simulation.R", "pop-plots.R", "pop-app.R",
               "hsp-simulation.R", "hsp-plots.R", "hsp-app.R", "ckmr-app.R")) {
  source(file.path("R", file), local = TRUE, encoding = "UTF-8")
}
shiny::shinyApp(ui = ckmr_app_ui(), server = ckmr_app_server)
