# Run from the project root: Rscript tests/check-laugh-calculation.R
for(f in c("pop-simulation.R","pop-plots.R","laugh-inputs.R","laugh-calculation.R","laugh-app.R"))
  source(file.path("R",f))
expect_error <- function(code,pattern=NULL) {
  e <- tryCatch({force(code);NULL},error=identity)
  stopifnot(inherits(e,"error"))
  if(!is.null(pattern)) stopifnot(grepl(pattern,conditionMessage(e),fixed=TRUE))
}

# Separate proportions from counts, retain zero-effort gap years, and remove
# display dimnames only in the package call (needed for its age-zero offarray).
s <- laugh_default_curves()
s$sample_scores[,] <- 0
s$sample_scores[cbind(1:3,c(1,6,13))] <- c(2,3,9)
design <- laugh_prepare_inputs(s,c(100,0,300))
args <- laugh_function_arguments(design)
stopifnot(identical(dim(args$samp_nya),c(3L,13L)),is.null(dimnames(args$samp_nya)),
  args$samp_nya[1,1]==100,args$samp_nya[3,13]==300,all(args$samp_nya[2,]==0),
  !"samp_pata" %in% names(args),!"samp_ny" %in% names(args),
  !is.null(dimnames(design$samp_pata)))
bad <- design; bad$samp_pata[1,1] <- .5
expect_error(laugh_function_arguments(bad))

if(!requireNamespace("microscoping",quietly=TRUE) ||
   !"samp_nya" %in% names(formals(microscoping::ckmr_laugh_test))) {
  cat("Matrix conversion passed; native package checks skipped (samp_nya support unavailable).\n")
} else {
  # Equal age profiles must reproduce the legacy margins interface numerically.
  shared <- laugh_prepare_inputs(laugh_default_curves(),c(100,500,100))
  legacy <- shared; legacy$samp_pata <- unname(shared$samp_pata[1,])
  native_legacy <- do.call(microscoping::ckmr_laugh_test,legacy)
  run <- laugh_run_design(shared)
  stopifnot(max(abs(run$expected - as.numeric(native_legacy)))<1e-9,
    abs(sum(run$expected[c("nsamp_J","nsamp_A")])-700)<1e-9,
    identical(names(attributes(run$expected)),"names"),!is.environment(run$preparation))
  # Compare truly different row profiles, including plus-group allocation,
  # against a direct native call. In this design all effort goes to juveniles.
  varying <- laugh_run_design(design)
  direct <- do.call(microscoping::ckmr_laugh_test,args)
  stopifnot(max(abs(varying$expected-as.numeric(direct)))<1e-9,
    abs(sum(varying$expected[c("nsamp_J","nsamp_A")])-400)<1e-9)
  # Preflight population expansion must match the package's own retained data.
  gory_args <- run$function_arguments; gory_args$want_gory_details <- TRUE
  gory <- do.call(microscoping::ckmr_laugh_test,gory_args)
  stopifnot(max(abs(run$preparation$expanded_nata-as.numeric(attr(gory,"gory")$nata)))<1e-9)
  rm(gory)
  # Documented POP cartoon: 100 offspring and 100 adults, 400 total adults,
  # all sampled after spawning, predict approximately 50 POPs.
  cartoon <- list(nata=c(4e6,400,0),wata=rep(1,3),pmatata=c(0,1,1),
    samp_pata=matrix(c(.5,.5,0),nrow=1),samp_ny=200,MAX_JUVE_AGE=0L,
    MIN_CATCH_CURVE_AGE=0L,beqy_POP=TRUE,etriv=.01,extraplus=FALSE)
  answer <- laugh_run_design(cartoon)
  stopifnot(abs(answer$expected[["POP_offJ"]]-50)<.001,answer$expected[["HSP_JJ"]]==0,
    answer$expected[["nsamp_J"]]==100,answer$expected[["nsamp_A"]]==100)
  # Large computation requests fail before allocating the package's 5D arrays.
  expect_error(laugh_preflight(args,max_work_cells=1),"interactive calculation budget")
  rising <- shared; rising$nata <- 100*1.3^(0:12)
  expect_error(laugh_run_design(rising),"catch-curve")
  invisible(ggplot2::ggplot_build(laugh_plot_results(run)))

  shiny::testServer(laugh_app_server,{
    session$setInputs(laugh_years="100,0,300",laugh_juvenile=2,laugh_catch_min=1,laugh_calculate=1)
    stopifnot(length(calculations())==1,nchar(calculation_error())==0)
    first <- calculations()[[1]]
    stopifnot(abs(sum(first$function_arguments$samp_nya)-400)<1e-9,
      first$package_version==as.character(utils::packageVersion("microscoping")))
    session$setInputs(laugh_sampling_year="3")
    session$setInputs(laugh_curve_action=list(key="samp_pata",year=3,action="preset",preset="increasing",revision=revision(),sequence=1))
    stopifnot(identical(first,calculations()[[1]]),
      grepl("differ",output$laugh_result_status$html,fixed=TRUE))
    session$setInputs(laugh_calculate=2)
    stopifnot(length(calculations())==2,identical(first,calculations()[[1]]))
    session$setInputs(laugh_result_id="1")
    stopifnot(identical(selected_calculation(),first))
    # A failed new fit leaves saved results available and clearly identified.
    session$setInputs(laugh_curve_edit=list(key="nata",values=as.list(100*1.3^(0:12)),revision=revision(),sequence=2))
    session$setInputs(laugh_calculate=3)
    stopifnot(length(calculations())==2,nchar(calculation_error())>0,
      identical(selected_calculation(),first))
  })
  cat("Native matrix/margins agreement, POP cartoon, plus-group expansion, budgets, and immutable saved calculations passed.\n")
}
