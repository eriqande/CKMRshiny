# Run from the project root: Rscript tests/check-laugh-editor.R
for (f in c("pop-simulation.R","pop-plots.R","pop-app.R","hsp-simulation.R",
            "hsp-plots.R","hsp-app.R","laugh-inputs.R","laugh-calculation.R","laugh-app.R","ckmr-app.R")) source(file.path("R",f))
expect_error <- function(code) stopifnot(inherits(tryCatch({force(code);NULL},error=identity),"error"))

s <- laugh_default_curves(12,1e6)
stopifnot(s$plus_age==12,all(lengths(s$curves)==13),abs(sum(s$curves$nata)-1e6)<1e-8,
  abs(tail(s$curves$nata,1)/s$curves$nata[12] - .75/(1-.75))<1e-12)
years <- laugh_parse_years("500, 0; 1000\n250")
stopifnot(identical(years,c(500,0,1000,250)))
for (bad in c("", "0,0", "1,-1", "1,0.5", "1,what", "Inf")) expect_error(laugh_parse_years(bad))
expect_error(laugh_default_curves(2)); expect_error(laugh_default_curves(100))
expect_error(laugh_validate_curve(s,"pmatata",rep(2,13)))
expect_error(laugh_validate_curve(s,"wata",rep(-1,13)))
expect_error(laugh_validate_curve(s,"nata",1:4))
inputs <- laugh_prepare_inputs(s,years,juvenile=2,catch_min=1)
stopifnot(is.matrix(inputs$samp_pata),identical(dim(inputs$samp_pata),c(4L,13L)),
  all(abs(rowSums(inputs$samp_pata)-1)<1e-12),identical(rownames(inputs$samp_pata),paste("Year",1:4)),
  identical(colnames(inputs$samp_pata),c(as.character(0:11),"12+")),
  inputs$MAX_JUVE_AGE==2,!inputs$beqy_POP,!inputs$extraplus,identical(inputs$samp_ny,years))
scaled <- s; scaled$sample_scores <- scaled$sample_scores*2
stopifnot(identical(inputs$samp_pata,laugh_prepare_inputs(scaled,years)$samp_pata))
a <- laugh_sample_allocation(inputs$samp_pata,years)
totals <- tapply(a$expected_samples,a$year,sum)
stopifnot(max(abs(totals-years))<1e-10, all(a$expected_samples[a$year==2]==0),
  abs(sum(a$expected_samples[a$age==12])-sum(inputs$samp_pata[,13]*years))<1e-10)
for (key in names(.laugh_curve_info())) for (preset in c("example","flat","increasing","declining","bell","clear")) {
  values <- laugh_curve_preset(s,key,preset)
  stopifnot(length(laugh_validate_curve(s,key,values))==13)
  if (key=="nata" && preset!="clear") stopifnot(abs(sum(values)-1e6)<1e-8)
}
zero <- s; zero$sample_scores[2,] <- 0
expect_error(laugh_prepare_inputs(zero,years))
zero <- s; zero$curves$pmatata[] <- 0
expect_error(laugh_prepare_inputs(zero,years))
expect_error(laugh_prepare_inputs(s,years,juvenile=12))
expect_error(laugh_prepare_inputs(s,years,catch_min=11))
# Different profiles in every row catch transposed matrix indexing and accidental
# normalization across all years. Zero effort still keeps a normalized row.
different <- s
different$sample_scores[,] <- 0
different$sample_scores[cbind(1:3,c(1,6,13))] <- c(2,9,4)
p <- laugh_prepare_inputs(different,c(100,0,300))$samp_pata
allocated <- laugh_sample_allocation(p,c(100,0,300))
stopifnot(sum(allocated$expected_samples[allocated$year==1 & allocated$age==0])==100,
  sum(allocated$expected_samples[allocated$year==3 & allocated$age==12])==300,
  p[2,6]==1,all(allocated$expected_samples[allocated$year==2]==0))
cached <- laugh_ensure_years(different,2)
stopifnot(identical(cached,different),identical(laugh_ensure_years(cached,3),different))
extended <- laugh_ensure_years(different,4,copy_year=2)
stopifnot(identical(unname(extended$sample_scores[4,]),unname(different$sample_scores[2,])))

shiny::testServer(laugh_app_server, {
  session$setInputs(laugh_years="100, 0, 300",laugh_juvenile=2,laugh_catch_min=1,laugh_after_spawning=FALSE)
  original <- state()
  stopifnot(length(arguments()$samp_ny)==3,abs(sum(allocation()$expected_samples)-400)<1e-10)
  # One completed stroke is one history entry. It normalizes the new scores.
  session$setInputs(laugh_curve_edit=list(key="samp_pata",values=as.list(0:12),revision=revision(),sequence=1))
  stopifnot(identical(unname(state()$sample_scores[1,]),as.numeric(0:12)),length(history())==1,
    identical(state()$sample_scores[2:3,],original$sample_scores[2:3,]),
    abs(arguments()$samp_pata[1,13]-12/sum(0:12))<1e-12)
  after_stroke <- state()
  session$setInputs(laugh_curve_edit=list(key="samp_pata",values=as.list(rep(99,13)),revision=0,sequence=2))
  stopifnot(identical(state(),after_stroke)) # stale browser updates cannot erase edits
  session$setInputs(laugh_undo=1)
  stopifnot(identical(state(),original),length(history())==0)
  session$setInputs(laugh_curve_action=list(key="pmatata",action="scale",factor=2,revision=revision(),sequence=3))
  stopifnot(all(state()$curves$pmatata<=1))
  before_grid <- state()
  session$setInputs(laugh_plus_age=6,laugh_apply_age=1)
  stopifnot(state()$plus_age==6,all(lengths(state()$curves)==7))
  session$setInputs(laugh_undo=2)
  stopifnot(identical(state(),before_grid))
  session$setInputs(laugh_abundance=234567,laugh_scale_abundance=1)
  stopifnot(abs(sum(state()$curves$nata)-234567)<1e-8)
  session$setInputs(laugh_curve_action=list(key="samp_pata",action="preset",preset="clear",revision=revision(),sequence=4))
  stopifnot(inherits(prepared(),"error"),grepl("Inputs need attention",output$laugh_readiness,fixed=TRUE))
  session$setInputs(laugh_undo=3)
  stopifnot(!inherits(prepared(),"error"))
  session$setInputs(laugh_years="100, wrong")
  stopifnot(inherits(prepared(),"error"))
})

shiny::testServer(laugh_app_server,{
  session$setInputs(laugh_years="100,0,300",laugh_juvenile=2,laugh_catch_min=1)
  session$setInputs(laugh_sampling_year="3")
  custom <- as.numeric((1:13)^2)
  session$setInputs(laugh_curve_edit=list(key="samp_pata",year=3,values=as.list(custom),revision=revision(),sequence=1))
  saved <- state()
  session$setInputs(laugh_years="100,0")
  stopifnot(selected_year()==2,nrow(state()$sample_scores)==3,
    identical(unname(state()$sample_scores[3,]),custom),nrow(arguments()$samp_pata)==2,
    nrow(table_data())==26,all(table_data()$allocated_samples[table_data()$year==2]==0))
  # Copy only the active rows; the hidden third year is not overwritten.
  session$setInputs(laugh_sampling_year="1",laugh_copy_years=1)
  stopifnot(identical(state()$sample_scores[3,],saved$sample_scores[3,]))
  session$setInputs(laugh_years="100,0,300")
  stopifnot(identical(unname(state()$sample_scores[3,]),custom),nrow(arguments()$samp_pata)==3,
    grepl("restored",output$laugh_edit_status,fixed=TRUE))
  # A fourth, previously unseen year copies the active third curve.
  session$setInputs(laugh_sampling_year="3")
  session$setInputs(laugh_years="100,0,300,900")
  stopifnot(identical(unname(state()$sample_scores[4,]),custom))
  before_counts <- state()
  session$setInputs(laugh_years="10,0,30,90")
  stopifnot(identical(state(),before_counts))
  # Copy-to-all is atomic and undo restores each year, not just the active row.
  session$setInputs(laugh_copy_years=2)
  stopifnot(all(state()$sample_scores[1,]==custom),all(state()$sample_scores[2,]==custom))
  session$setInputs(laugh_undo=1)
  stopifnot(identical(state(),before_counts))
  # Undoing an older curve edit must retain a newer cached year.
  session$setInputs(laugh_undo=2)
  stopifnot(nrow(state()$sample_scores)==4,identical(unname(state()$sample_scores[4,]),custom))
  # Resetting the grid includes hidden rows; Undo recovers all prior cached values.
  session$setInputs(laugh_years="10,0")
  before_grid <- state()
  session$setInputs(laugh_plus_age=30,laugh_apply_age=1)
  stopifnot(identical(dim(state()$sample_scores),c(4L,31L)))
  session$setInputs(laugh_undo=3)
  stopifnot(identical(state(),before_grid))
})

shiny::testServer(ckmr_app_server, {
  session$setInputs(laugh_years="100,0,200",laugh_juvenile=2,laugh_catch_min=1)
  session$setInputs(laugh_curve_edit=list(key="wata",values=as.list(seq(.1,6,length.out=13)),revision=laugh$revision(),sequence=1))
  saved <- laugh$state()
  session$setInputs(open_pop=1,app_view="pop")
  session$setInputs(open_laugh=1,app_view="laugh")
  stopifnot(identical(laugh$state(),saved), .pop_usage_total(usage,"samples")==0)
})
ui <- as.character(ckmr_app_ui())
ids <- regmatches(ui,gregexpr(' id="[^"]+"',ui))[[1]]
stopifnot(!anyDuplicated(ids))
cat("Year-by-age matrices, independent curves, cached years, copy/undo, allocations, and launcher state passed.\n")
