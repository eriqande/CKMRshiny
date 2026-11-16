# Run from the project root with the repaired microscoping installation.
source("app.R")
fail <- function(x) stopifnot(inherits(tryCatch({force(x);NULL},error=identity),"error"))
same <- function(x,y) stopifnot(isTRUE(all.equal(x,y,check.attributes=FALSE,tolerance=1e-10)))
stopifnot(identical(laugh_parse_years(.laugh_count_text(c(1/3,0,12345.678901234567))),c(1/3,0,12345.678901234567)))
one <- laugh_function_arguments(laugh_prepare_inputs(laugh_default_curves(6,1e5,3),c(100.25,0,300.5),juvenile=1,catch_min=0,after_spawning=TRUE))
one$etriv <- .05
loaded <- ckmr_import_design(one)
stopifnot(loaded$kind=="laugh",identical(loaded$states$laugh$curves$nata,one$nata),loaded$counts[2]==0)
x <- laugh_prepare_inputs(loaded$states$laugh,loaded$counts,1,0,TRUE); x$etriv <- .05
same(one,laugh_function_arguments(x))
legacy <- one;legacy$samp_nya <- NULL;legacy$samp_pata <- c(1,3,rep(0,5));legacy$samp_ny <- c(100.25,0,300.5)
z <- ckmr_import_design(legacy);same(z$counts,legacy$samp_ny)
stopifnot(z$states$laugh$sample_scores[1,2]==.75)
fail(ckmr_import_design(list(nata=1:5)))
fail(ckmr_import_design(within(one,{samp_nya[1,1] <- NA})))
fail(ckmr_import_design(within(one,{extraplus <- TRUE})))
fail(ckmr_import_design(within(one,{samp_nya <- t(samp_nya)})))

sex_state <- list(F=loaded$states$laugh,M=loaded$states$laugh,
  counts=rbind(F=c(100.25,0,300.5),M=c(50.5,200.75,0)),settings=loaded$settings)
sex_state$M$curves$nata <- sex_state$M$curves$nata * .6
sex_state$M$curves$wata <- sex_state$M$curves$wata^.7
args <- clt2_function_arguments(sex_state)
two <- ckmr_import_design(args)
restored <- list(F=two$states$F,M=two$states$M,counts=two$counts,settings=two$settings)
same(args,clt2_function_arguments(restored))
# Every sex/year must keep its effort, including zero rows and fractional totals.
same(apply(args$nsamp_sya,c(1,2),sum),sex_state$counts)
rev_args <- args
for(k in c("n_sa","w_sa","pmat_sa")) rev_args[[k]] <- args[[k]][2:1,]
rev_args$nsamp_sya <- args$nsamp_sya[2:1,,,drop=FALSE]
rev <- ckmr_import_design(rev_args)
same(clt2_function_arguments(list(F=rev$states$F,M=rev$states$M,counts=rev$counts,settings=rev$settings)),args)
margin <- args;margin$nsamp_sya <- NULL;margin$samp_p_sa <- matrix(seq_len(14),2,7,dimnames=list(c("F","M"),NULL));margin$nsamp_y <- c(123.25,0,55)
m <- ckmr_import_design(margin)
same(colSums(m$counts),margin$nsamp_y)
for(i in 1:2) same(m$counts[i,],margin$nsamp_y*sum(margin$samp_p_sa[i,])/sum(margin$samp_p_sa))
fail(ckmr_import_design(within(args,{nsamp_sya <- nsamp_sya[,,1:3]})))
fail(ckmr_import_design(within(args,{nsamp_sya[1,1,1] <- -1})))
extra <- one;extra$weird_argument <- 123;stopifnot(any(grepl("weird_argument",ckmr_import_design(extra)$notices)))
file <- tempfile(fileext=".rds");saveRDS(args,file);same(ckmr_import_design(readRDS(file)),two);unlink(file)

# Actual package agreement and guard before its large arrays are allocated.
r <- clt2_run_design(args)
same(r$expected,as.numeric(do.call(microscoping::clt2,args)))
g <- args;g$want_gory_details <- TRUE;g <- do.call(microscoping::clt2,g)
for(i in 1:2) same(r$preparation$sexes[[i]]$expanded_nata,as.numeric(attr(g,"gory")$n_sa[c("F","M")[i],]))
fail(clt2_preflight(args,max_work_cells=1))
ex <- args;ex$extraplus <- TRUE
r_ex <- clt2_run_design(ex);g_ex <- ex;g_ex$want_gory_details <- TRUE;g_ex <- do.call(microscoping::clt2,g_ex)
for(i in 1:2) same(r_ex$preparation$sexes[[i]]$expanded_nata,as.numeric(attr(g_ex,"gory")$n_sa[c("F","M")[i],]))

shiny::testServer(laugh_app_server,{
  session$setInputs(laugh_years="100,200,300",laugh_juvenile=2,laugh_catch_min=1)
  before <- state(); old <- arguments()
  import_design(loaded);session$flushReact()
  same(laugh_function_arguments(arguments()),one)
  session$setInputs(laugh_calculate=1)
  stopifnot(length(calculations())==1)
  snapshot <- calculations()[[1]]
  session$setInputs(laugh_curve_edit=list(key="wata",values=rep(1,7),revision=revision()))
  stopifnot(identical(snapshot,calculations()[[1]]))
  session$setInputs(laugh_undo=1)
  same(laugh_function_arguments(arguments()),one)
  session$setInputs(laugh_undo=2)
  same(state(),before);same(arguments(),old)
})
shiny::testServer(clt2_app_server,{
  before <- state();import_design(two);session$flushReact();same(arguments(),args)
  nh <- length(history())
  # Browser echoes of updated controls must not create spurious Undo entries.
  session$setInputs(clt2_juvenile=as.numeric(two$settings$MAX_JUVE_AGE),clt2_catch_min=as.numeric(two$settings$MIN_CATCH_CURVE_AGE),
    clt2_after_spawning=two$settings$beqy_POP,clt2_etriv=two$settings$etriv,clt2_extraplus=FALSE,
    clt2_years_F=.laugh_count_text(two$counts[1,]),clt2_years_M=.laugh_count_text(two$counts[2,]))
  stopifnot(length(history())==nh)
  session$setInputs(clt2_calculate=1)
  stopifnot(length(calculations())==1,calculation_error()=="")
  saved <- calculations()[[1]]
  old_m <- state()$M
  session$setInputs(clt2_F_curve_edit=list(key="wata",values=rep(2,7),revision=revision()))
  same(state()$M,old_m);stopifnot(identical(saved,calculations()[[1]]))
  session$setInputs(clt2_undo=1);same(arguments(),args)
  session$setInputs(clt2_undo=2);same(state(),before)
  # Hidden curves survive a reduction followed by re-adding the same years.
  import_design(two);session$flushReact()
  session$setInputs(clt2_years_F="10,20",clt2_years_M="0,30")
  stopifnot(nrow(state()$F$sample_scores)==3,ncol(state()$counts)==2)
  session$setInputs(clt2_years_F="10,20,40",clt2_years_M="0,30,50")
  same(state()$F$sample_scores[3,],two$states$F$sample_scores[3,])
})
# Shared upload handling must route correctly, and failures must be atomic.
f <- tempfile(fileext=".rds");saveRDS(args,f)
shiny::testServer(ckmr_app_server,{
  session$setInputs(design_file=list(datapath=f),design_import=1)
  same(clt2$arguments(),args);old_clt <- clt2$state();old_one <- laugh$state()
  saveRDS(list(n_sa="bad"),f);session$setInputs(design_import=2)
  same(clt2$state(),old_clt);same(laugh$state(),old_one)
  stopifnot(grepl("not applied",import_status()))
  saveRDS(one,f);session$setInputs(design_import=3)
  same(laugh_function_arguments(laugh$arguments()),one);same(clt2$state(),old_clt)
})
unlink(f)
startup_server <- function(input,output,session) { handles <- ckmr_app_server(input,output,session,design_inputs=args) }
shiny::testServer(startup_server,{
  session$flushReact();same(handles$clt2$arguments(),args)
})
cat("One/two-sex import round trips, package agreement, preflight, independent edits, immutable results, undo, and atomic upload handling passed.\n")

# Optional package extension must preserve totals and saved parent-sex metadata.
if(.clt2_parent_sex_supported()) {
  r <- clt2_run_design(args)
  m <- r$parent_sex_expected
  stopifnot(identical(dim(m),c(2L,3L)),isTRUE(r$function_arguments$by_parent_sex))
  same(colSums(m),r$expected[colnames(m)])
  tab <- clt2_parent_sex_table(r)
  stopifnot(nrow(tab)==6,identical(unique(tab$parent_sex),c("F","M")))
  invisible(ggplot2::ggplot_build(clt2_plot_parent_sex(r)))
  # Saved-calculation import may include the reporting flag; it is not a new
  # demographic input and must not produce a spurious 'unknown argument' notice.
  restored <- ckmr_import_design(r)
  stopifnot(!any(grepl("by_parent_sex",restored$notices,fixed=TRUE)))
  shiny::testServer(clt2_app_server,{
    import_design(two);session$flushReact();session$setInputs(clt2_calculate=1)
    saved <- selected_calculation()
    stopifnot(!is.null(saved$parent_sex_expected),grepl("clt2_parent_sex_plot",output$clt2_parent_sex_controls$html,fixed=TRUE))
    session$setInputs(clt2_F_curve_edit=list(key="wata",values=rep(4,7),revision=revision()))
    stopifnot(identical(saved,selected_calculation()))
  })
  cat("Optional native parent-sex accounting, column totals, plotting and immutable metadata passed.\n")
} else {
  r <- clt2_run_design(args)
  stopifnot(is.null(r$parent_sex_expected),is.null(clt2_parent_sex_table(r)))
  cat("Combined-count fallback passed (installed package lacks parent-sex reporting).\n")
}
