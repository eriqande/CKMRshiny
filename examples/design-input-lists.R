# Small illustrative lists, not a stock assessment. Run from the repo root:
# source("app.R")
# source("examples/design-input-lists.R")
# shiny::runApp(ckmr_design_app(inputs = clt2_inputs))
# Or saveRDS(clt2_inputs, "design.rds") and use Import inputs in the app.

a <- 0:6
n <- c(.65^(0:5),.65^6/(1-.65))
w <- .1 + (1-exp(-.3*(a+.2)))^3 * 10
maturity <- plogis((a-2.5)/.7); maturity[1] <- 0

# Rows are consecutive years, columns are ages 0,...,5,6+.
samples <- rbind(c(20,30,30,10,5,3,2),rep(0,7),c(5,10,15,20,20,20,10))
laugh_inputs <- list(nata=1e5*n/sum(n),wata=w,pmatata=maturity,
  samp_nya=samples,MAX_JUVE_AGE=2,MIN_CATCH_CURVE_AGE=0,
  beqy_POP=FALSE,etriv=.01,extraplus=FALSE)

# Sex x sampling year x age COUNTS. Each sex can have its own sampling mix.
sex_samples <- array(0,c(2,3,7),dimnames=list(SEX=c("F","M"),YEAR=NULL,AGE=NULL))
sex_samples["F",,] <- samples/2
sex_samples["M",,] <- samples/2
sex_samples["M",2,] <- c(0,0,0,10,10,10,20)
clt2_inputs <- list(n_sa=rbind(F=laugh_inputs$nata/2,M=laugh_inputs$nata/2),
  w_sa=rbind(F=w,M=w^.7),pmat_sa=rbind(F=maturity,M=maturity),
  nsamp_sya=sex_samples,MAX_JUVE_AGE=2,MIN_CATCH_CURVE_AGE=0,
  beqy_POP=FALSE,etriv=.01,extraplus=FALSE)

# These same lists also work outside Shiny:
# do.call(microscoping::ckmr_laugh_test, laugh_inputs)
# do.call(microscoping::clt2, clt2_inputs)

# With the parent-sex reporting extension in microscoping:
# result <- do.call(microscoping::clt2, c(clt2_inputs, list(by_parent_sex=TRUE)))
# attr(result, "kin_by_parent_sex") # F/M rows; POP_offJ, POP_offA, HSP_JJ columns
