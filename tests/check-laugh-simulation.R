# Standalone: no app or other repository R files are sourced.
source("R/laugh-simulation.R")
expect_error <- function(expr,text) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,"error"),grepl(text,conditionMessage(e),fixed=TRUE))
}

# Hand-built pedigree deliberately includes full siblings across cohorts,
# multiple grandparent paths, and individuals with several relationship roles.
p <- data.frame(id=1:7,birth_year=c(-4,-4,-2,-1,0,0,1),
  sex=c("F","M","F","M","F","M","F"),
  mother=c(0,0,1,1,3,3,3),father=c(0,0,2,2,2,4,4))
pair <- function(ids,years,mja=10,after=FALSE) {
  s <- data.frame(id=ids,year=years,age=years-p$birth_year[ids])
  count_laugh_kin(p,s,mja,after)
}
stopifnot(pair(c(3,4),c(1,1))["HSP_JJ"]==2,
  pair(c(5,6),c(1,1))["HSP_JJ"]==0,
  pair(c(1,6),c(1,1))["GGP_JJ"]==2,
  pair(c(1,6),c(-2,1))["GGP_JJ"]==0,
  pair(c(1,6),c(-1,1))["GGP_JJ"]==1,
  pair(c(3,6),c(0,1),0,FALSE)["POP_offA"]==0,
  pair(c(3,6),c(0,1),0,TRUE)["POP_offA"]==1,
  pair(c(3,7),c(0,2),0,TRUE)["POP_offA"]==0,
  pair(c(6,6),c(1,2))["HSP_JJ"]==0,
  pair(c(1,6),c(1,1),0)["GGP_JJ"]==0)
expect_error(pair(c(6,6),c(1,1)),"repeated samples")
s <- data.frame(id=c(1,3,4,5,6,7),year=2)
s$age <- s$year-p$birth_year[s$id]
x <- count_laugh_kin(p,s,10)
backwards <- count_laugh_kin(p,s[nrow(s):1,],10)
stopifnot(identical(x,backwards))
cat("Hand-built pedigree counting and timing conventions passed.\n")

if(!requireNamespace("microscoping",quietly=TRUE)) {
  cat("Demographic and Monte Carlo checks skipped: microscoping is not installed.\n")
} else {
  d <- list(nata=c(512,256,128,64,64),wata=c(0,1,1,1,1),
    pmatata=c(0,1,1,1,1),samp_nya=rbind(c(10,6,6,6,0),c(10,6,6,6,0),c(10,6,6,6,0)),
    MAX_JUVE_AGE=2)
  prep <- do.call(prepare_laugh_simulation,d)
  stopifnot(max(abs(prep$ages$package_annual_survival-.5))<1e-10,
    max(abs(2*prep$per_sex-prep$ages$target_both_sexes))<1e-8,
    tail(prep$ages$simulated_annual_survival,1)==0,
    all(prep$samples$expected_per_sex>=0))
  # Independently compare demographic preparation with the package's retained
  # intermediate objects. This is validation only, never simulation input.
  args <- prep$arguments; args$want_gory_details <- TRUE
  direct <- do.call(microscoping::ckmr_laugh_test,args)
  gory <- attr(direct,"gory")
  stopifnot(max(abs(as.numeric(gory$nata)-prep$ages$target_both_sexes))<1e-9,
    max(abs(as.numeric(gory$fecata)-prep$ages$relative_output))<1e-9,
    max(abs(as.numeric(gory$fsamp_ya)-as.numeric(prep$allocation)))<1e-9)
  rm(gory,direct)
  # Different profiles, a zero-effort year, fractional allocations, and plus
  # samples: expansion preserves the supplied totals and splits the plus group.
  gap <- d; gap$samp_nya <- rbind(c(5,0,0,0,2),rep(0,5),c(0,3,2,1,0))
  pg <- do.call(prepare_laugh_simulation,gap)
  stopifnot(max(abs(2*rowSums(pg$allocation)-rowSums(gap$samp_nya)))<1e-10,
    all(pg$allocation[2,]==0),all(pg$allocation[1,5:ncol(pg$allocation)]>0))
  legacy <- d; legacy$samp_nya <- NULL
  legacy$samp_pata <- c(10,6,6,6,0); legacy$samp_ny <- rep(28,3)
  pl <- do.call(prepare_laugh_simulation,legacy)
  stopifnot(identical(pl$allocation,prep$allocation))
  # All budget failures occur before any population is simulated.
  expect_error(do.call(prepare_laugh_simulation,c(d,list(max_population=10))),"max_population")
  expect_error(do.call(prepare_laugh_simulation,c(d,list(max_individuals=1000))),"max_individuals")
  expect_error(do.call(prepare_laugh_simulation,c(d,list(max_pair_comparisons=1))),"max_pair_comparisons")
  expect_error(do.call(simulate_laugh_test,c(d,list(max_total_individuals=1000))),"max_total_individuals")
  expect_error(do.call(simulate_laugh_test,c(d,list(max_total_comparisons=1))),"max_total_comparisons")
  expect_error(check_laugh_test(d,num_reps=2,max_work_cells=1),"max_work_cells")
  expect_error(check_laugh_test(c(d,list(seed=3))),"without simulation controls")
  bad <- d; bad$samp_nya[1,1] <- 1e4
  expect_error(do.call(prepare_laugh_simulation,bad),"requested sample")
  bad <- d; bad$pmatata[1] <- 1; bad$wata[1] <- 1
  expect_error(do.call(prepare_laugh_simulation,bad),"age zero")
  expect_error(do.call(prepare_laugh_simulation,c(d,list(extraplus=TRUE))),"extraplus")
  # Local seed restoration, reproducibility, and fixed census reconstruction.
  set.seed(73); old <- .Random.seed; old_kind <- RNGkind()
  a <- do.call(simulate_laugh_test,c(d,list(num_reps=2,seed=19)))
  stopifnot(identical(old,.Random.seed),identical(old_kind,RNGkind()),
    identical(a,do.call(simulate_laugh_test,c(d,list(num_reps=2,seed=19)))))
  pp <- a$first$pedigree; ss <- a$first$samples
  # Reconstruct every census in the sampling period from births and deaths.
  for(y in 1:3) for(sex in c("F","M")) {
    present <- pp$sex==sex & pp$birth_year<=y & (is.na(pp$death_year) | pp$death_year>y)
    by_age <- tabulate(y-pp$birth_year[present]+1,nbins=nrow(prep$ages))
    stopifnot(identical(as.integer(by_age),prep$per_sex))
  }
  stopifnot(all(pp$mother[ss$id]>0),all(pp$father[ss$id]>0),
    all(pp$mother[pp$mother[ss$id]]>0),all(pp$father[pp$father[ss$id]]>0),
    !anyDuplicated(ss[c("id","year")]),
    all(pp$birth_year[pp$mother[ss$id]]<pp$birth_year[ss$id]),
    all(pp$birth_year[pp$father[ss$id]]<pp$birth_year[ss$id]),
    all(a$counts$nsamp_J+a$counts$nsamp_A==84))
  expect_error(.laugh_sim_seed(3,stop("intentional")),"intentional")
  stopifnot(identical(old,.Random.seed))
  rm(".Random.seed",envir=.GlobalEnv)
  invisible(.laugh_sim_seed(3,runif(1)))
  stopifnot(!exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  assign(".Random.seed",old,.GlobalEnv)
  # Fixed-seed smoke check of all four kin summaries. Small rounding/tail effects
  # remain explicit; broad MC bounds detect major counting or orientation errors.
  result <- check_laugh_test(d,num_reps=500,seed=1401,keep_first=FALSE)
  print(result$comparison)
  z <- result$comparison$z_difference[1:4]
  stopifnot(all(is.finite(z)),all(abs(z)<5),is.null(result$simulation$first))
  # Also exercise non-geometric input abundances, changing output at age,
  # nonzero maturity fractions, and the alternate same-year POP convention.
  nonlinear <- list(nata=c(800,500,300,160,120),wata=c(0,.5,1,2,3),
    pmatata=c(0,.2,.7,1,1),samp_nya=rbind(c(10,8,6,4,0),c(8,8,8,4,0)),
    MAX_JUVE_AGE=2,beqy_POP=TRUE)
  curved <- check_laugh_test(nonlinear,num_reps=150,seed=1703,keep_first=FALSE)
  stopifnot(all(is.finite(curved$comparison$z_difference[1:4])),
    all(abs(curved$comparison$z_difference[1:4])<5))
  # Fractional sample sizes are unbiased in expectation (not rounded once).
  frac <- do.call(simulate_laugh_test,c(gap,list(num_reps=200,seed=92,keep_first=FALSE)))
  totals <- frac$counts$nsamp_J+frac$counts$nsamp_A
  se <- sd(totals)/sqrt(length(totals))
  stopifnot(abs(mean(totals)-sum(gap$samp_nya))<5*se)
  cat("Demographic preparation, limits, RNG isolation, pedigree simulation, and Monte Carlo comparisons passed.\n")
}
