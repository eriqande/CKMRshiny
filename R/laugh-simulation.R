# Standalone verification of microscoping's single-schedule laugh test.
# Source this file alone: no Shiny or other files in this repository are needed.
# Only demographic helpers are borrowed. Pedigrees and relationship counts are
# simulated independently; no microscoping kin probabilities are used here.

#' Prepare a small stationary population for a laugh-test simulation
#'
#' Matches the demographic preparation in microscoping 1.1.2, extraplus=FALSE.
#' Uses its private smooth_n() and get_psurv() helpers. The resulting exact-age
#' abundances are rounded to whole individuals PER SEX. Each year, random subsets
#' survive into fixed age/sex targets; recruitment is fixed. This preserves a
#' constant reproductive denominator and avoids demographic drift. The population
#' is stationary exactly after rounding, rather than only in expectation.
#'
#' Age zero must have zero reproductive output. At the final expanded age all
#' individuals die before the next census. microscoping's survival helper can
#' extrapolate beyond this age, so its HSP tail calculation is not exactly the
#' same finite population. Rounding and terminal-age diagnostics are returned:
#' this is an approximate benchmark, not a claim of exact generative equivalence.
#' @param nata Population numbers, both sexes, ages zero through the plus group.
#' @param wata,pmatata Biological weight and mature proportion at those ages.
#' Reproductive output is their product, shared by all individuals of an age;
#' no individual maturity state or extra reproductive heterogeneity is drawn.
#' @param samp_pata,samp_ny Optional legacy sample-age proportions and annual
#' totals. Ignored when samp_nya is supplied, as in microscoping.
#' @param samp_nya Optional year-by-age matrix of expected sample NUMBERS.
#' Its last column is a plus-group total. Dimension labels are display-only.
#' @param MAX_JUVE_AGE Oldest age in juvenile HSP/GGP and POP offspring summaries.
#' @param MIN_CATCH_CURVE_AGE Youngest age included in the catch-curve fit.
#' @param beqy_POP Include parents sampled in the offspring's birth year.
#' @param etriv Plus-group expansion tolerance, as in microscoping.
#' @param extraplus Only FALSE is supported in this initial verification tool.
#' @param want_gory_details Accepted for compatibility; must be FALSE. The
#' simulator returns ordinary inspectable objects instead of package environments.
#' @param burn_in Years before year 1. NULL uses three expanded lifespans,
#' enough to replace founders in the sampled individuals' grandparent pedigrees.
#' @param max_population Maximum standing population (both sexes), default 20000.
#' @param max_individuals Maximum pedigree records in ONE replicate, default 200000.
#' @param max_pair_comparisons Maximum sample-event pairs per replicate, default
#' 500000. This intentionally small tool counts pairs directly.
#' @return List containing canonical package arguments, inspectable age and
#' sample tables, fixed census targets, work estimates, and diagnostics. The
#' package source is never modified. Internal helper changes need revalidation.
#' @export
prepare_laugh_simulation <- function(nata, wata, pmatata, samp_pata=NULL,
    samp_ny=NULL, samp_nya=NULL, MAX_JUVE_AGE, MIN_CATCH_CURVE_AGE=0,
    beqy_POP=FALSE, etriv=.01, extraplus=FALSE, want_gory_details=FALSE,
    burn_in=NULL, max_population=20000, max_individuals=200000,
    max_pair_comparisons=500000) {
  if (!requireNamespace("microscoping",quietly=TRUE) ||
      !requireNamespace("offarray",quietly=TRUE))
    stop("Install microscoping and offarray before preparing a simulation.")
  if (!"samp_nya" %in% names(formals(microscoping::ckmr_laugh_test)))
    stop("This check requires microscoping with samp_nya support.")
  for (x in c("max_population","max_individuals","max_pair_comparisons"))
    .laugh_sim_integer(get(x),x,1)
  n <- length(nata); a <- n-1L
  if (n<3L || n>61L) stop("Use 3 to 61 age categories, including the plus group.")
  for (x in c("nata","wata","pmatata")) {
    v <- get(x)
    if (!is.numeric(v) || length(v)!=n || any(!is.finite(v)) || any(v<0))
      stop(x," must contain finite nonnegative values at every age.")
  }
  if (any(pmatata>1)) stop("pmatata must be at most one.")
  if (sum(nata)<=0 || sum(nata)>max_population)
    stop("Input population exceeds max_population or is empty; use a small example.")
  if (wata[1]*pmatata[1]!=0)
    stop("This simulator requires zero reproductive output at age zero.")
  .laugh_sim_integer(MAX_JUVE_AGE,"MAX_JUVE_AGE",0,a-1)
  .laugh_sim_integer(MIN_CATCH_CURVE_AGE,"MIN_CATCH_CURVE_AGE",0,a-2)
  if (!identical(extraplus,FALSE) || !identical(want_gory_details,FALSE))
    stop("Use extraplus=FALSE and want_gory_details=FALSE for this check.")
  if (!is.logical(beqy_POP) || length(beqy_POP)!=1 || is.na(beqy_POP))
    stop("beqy_POP must be TRUE or FALSE.")
  if (length(etriv)!=1 || !is.finite(etriv) || etriv<=0 || etriv>=1)
    stop("etriv must lie strictly between zero and one.")
  if (is.null(samp_nya)) {
    if (!is.numeric(samp_pata) || length(samp_pata)!=n ||
        any(!is.finite(samp_pata)) || any(samp_pata<0) || sum(samp_pata)<=0 ||
        !is.numeric(samp_ny) || !length(samp_ny) ||
        any(!is.finite(samp_ny)) || any(samp_ny<0))
      stop("Supply samp_nya, or valid samp_pata and samp_ny margins.")
    samp_nya <- outer(samp_ny,samp_pata/sum(samp_pata))
  }
  if (!is.matrix(samp_nya) || !is.numeric(samp_nya) || ncol(samp_nya)!=n ||
      nrow(samp_nya)<1 || nrow(samp_nya)>12 || any(!is.finite(samp_nya)) ||
      any(samp_nya<0) || sum(samp_nya)<=0)
    stop("samp_nya must be a nonnegative year-by-age count matrix, with 1 to 12 years.")
  samp_nya <- matrix(as.numeric(samp_nya),nrow=nrow(samp_nya),ncol=ncol(samp_nya))
  ny <- nrow(samp_nya)
  adjusted <- as.numeric(nata); adjusted[adjusted==0] <- .001
  indexed <- offarray::offarray(adjusted,
    first=c(AGE=0),dim=n)
  smooth <- as.numeric(utils::getFromNamespace("smooth_n","microscoping")(
    indexed[MIN_CATCH_CURVE_AGE:(a-1L)]))
  plus <- if(nata[n]==0) .001 else nata[n]
  ps <- plus/(tail(smooth,1)+plus)
  f <- as.numeric(wata*pmatata)
  pre <- sum(head(f,-1)*smooth); post <- tail(f,1)*plus
  if (any(!is.finite(smooth)) || any(diff(smooth)>0) ||
      !is.finite(pre+post) || pre+post<=0 || !is.finite(ps) || ps<=0 || ps>=1)
    stop("The fitted population must decline and have positive reproductive output.")
  extra <- floor(1+max(log(etriv*(1+pre/post))/log(ps),0))
  oldest <- a+extra
  if (!is.finite(oldest) || oldest>200)
    stop("The expanded population exceeds 201 age classes; use a shorter-lived example.")
  target <- c(smooth,plus*(1-ps)*ps^(0:extra))
  if(any(!is.finite(target)) || any(target<=0))
    stop("The expanded age schedule contains nonpositive or nonfinite abundances.")
  fec <- c(head(f,-1),rep(tail(f,1),extra+1))
  per_sex <- round(target/2)
  if (sum(target)>max_population || 2*sum(per_sex)>max_population)
    stop("The smoothed/expanded population exceeds max_population.")
  if (per_sex[1]<1 || any(diff(per_sex)>0) || sum(per_sex*fec)<=0)
    stop("Rounded age targets need recruits, nonincreasing numbers, and breeders.")
  if (is.null(burn_in)) burn_in <- 3L*(oldest+1L)
  .laugh_sim_integer(burn_in,"burn_in",3L*(oldest+1L))
  records <- 2*sum(per_sex)+2*per_sex[1]*(burn_in+ny-1)
  if (!is.finite(records) || records>max_individuals || records>.Machine$integer.max)
    stop("Pedigree size exceeds max_individuals after burn-in; use a smaller example.")
  # Retain fractional expected allocations. Stochastic rounding below preserves
  # each sex/year/age expectation, even for odd totals and plus-group fractions.
  allocation <- cbind(samp_nya[,seq_len(a),drop=FALSE],
    as.numeric(samp_nya[,n]) %o% rep(1/(extra+1),extra+1)) / 2
  if (any(sweep(ceiling(allocation),2,per_sex,">") ))
    stop("A requested sample can exceed its rounded sex/age population; reduce sampling or increase that age class.")
  max_samples <- 2*sum(ceiling(allocation))
  pairs <- max_samples*(max_samples-1)/2
  if (pairs>max_pair_comparisons)
    stop("Sample-event pairs exceed max_pair_comparisons; reduce sampling.")
  oo <- offarray::offarray(target,first=c(AGE=0),dim=length(target))
  surv <- as.numeric(utils::getFromNamespace("get_psurv","microscoping")(oo,1)[,1])
  actual <- c(ifelse(head(per_sex,-1)>0,tail(per_sex,-1)/head(per_sex,-1),NA_real_),0)
  ages <- data.frame(age=0:oldest,target_both_sexes=target,
    simulated_per_sex=per_sex,relative_output=fec,
    package_annual_survival=surv,simulated_annual_survival=actual)
  sample_table <- data.frame(year=rep(seq_len(ny),each=oldest+1),
    age=rep(0:oldest,ny),expected_per_sex=as.vector(t(allocation)))
  diagnostics <- c(
    relative_abundance_rounding=(2*sum(per_sex)-sum(target))/sum(target),
    relative_output_rounding=(2*sum(per_sex*fec)-sum(target*fec))/sum(target*fec),
    terminal_output_fraction=tail(target*fec,1)/sum(target*fec))
  args <- list(nata=as.numeric(nata),wata=as.numeric(wata),pmatata=as.numeric(pmatata),
    samp_nya=samp_nya,MAX_JUVE_AGE=MAX_JUVE_AGE,
    MIN_CATCH_CURVE_AGE=MIN_CATCH_CURVE_AGE,beqy_POP=beqy_POP,
    etriv=etriv,extraplus=FALSE,want_gory_details=FALSE)
  list(arguments=args,ages=ages,samples=sample_table,allocation=allocation,
    per_sex=as.integer(per_sex),burn_in=burn_in,num_years=ny,
    records_per_rep=records,max_pairs_per_rep=pairs,diagnostics=diagnostics,
    package_version=as.character(utils::packageVersion("microscoping")))
}

#' Simulate independent small pedigrees for a laugh-test design
#'
#' The biological census is just after reproduction. Individuals surviving from
#' the previous year age by one, then fixed numbers of newborns are assigned
#' independent mothers/fathers among current breeders, proportional to output.
#' Uniform thinning within age/sex maintains the fixed census schedule.
#'
#' Samples are nonlethal observations, without replacement within year/age/sex;
#' recaptures across years are possible. POP counts apply microscoping's birth-
#' year eligibility restriction explicitly. Thus this tests its restricted
#' comparison counts; it does not simulate mortality caused by lethal sampling.
#' GGP counts require the grandparent to have been sampled strictly after the
#' intermediate parent's birth, matching the package's separate GGP convention.
#' Each fractional sex/year/age sample size is independently stochastically
#' rounded, preserving its expectation.
#' Multiple shared-parent or grandparent paths are counted separately.
#'
#' @inheritParams prepare_laugh_simulation
#' @param num_reps Number of independent pedigrees, default 100.
#' @param seed Positive integer seed. Caller RNG kind and state are restored.
#' @param keep_first Retain the first replicate's pedigree and samples for inspection.
#' @param max_total_individuals Maximum pedigree records across all replicates,
#' default twenty million. Guards run before any population is allocated.
#' @param max_total_comparisons Maximum pair checks across all replicates,
#' default twenty million.
#' @return List with replicate counts, prepared demography, seed, and optional
#' first replicate. No analytic kin probabilities are computed.
#' @export
simulate_laugh_test <- function(nata,wata,pmatata,samp_pata=NULL,samp_ny=NULL,
    samp_nya=NULL,MAX_JUVE_AGE,MIN_CATCH_CURVE_AGE=0,beqy_POP=FALSE,
    etriv=.01,extraplus=FALSE,want_gory_details=FALSE,num_reps=100,seed=1401,
    burn_in=NULL,keep_first=TRUE,max_population=20000,max_individuals=200000,
    max_pair_comparisons=500000,max_total_individuals=2e7,
    max_total_comparisons=2e7) {
  .laugh_sim_integer(num_reps,"num_reps",1)
  .laugh_sim_integer(seed,"seed",1)
  .laugh_sim_integer(max_total_individuals,"max_total_individuals",1)
  .laugh_sim_integer(max_total_comparisons,"max_total_comparisons",1)
  if(!is.logical(keep_first) || length(keep_first)!=1 || is.na(keep_first))
    stop("keep_first must be TRUE or FALSE.")
  prep <- prepare_laugh_simulation(nata,wata,pmatata,samp_pata,samp_ny,samp_nya,
    MAX_JUVE_AGE,MIN_CATCH_CURVE_AGE,beqy_POP,etriv,extraplus,want_gory_details,
    burn_in,max_population,max_individuals,max_pair_comparisons)
  if(prep$records_per_rep*num_reps>max_total_individuals)
    stop("Batch exceeds max_total_individuals; reduce population, burn-in, or replicates.")
  if(prep$max_pairs_per_rep*num_reps>max_total_comparisons)
    stop("Batch exceeds max_total_comparisons; reduce sampling or replicates.")
  .laugh_sim_seed(seed,{
    counts <- matrix(0,num_reps,6,dimnames=list(NULL,
      c("POP_offJ","POP_offA","HSP_JJ","GGP_JJ","nsamp_J","nsamp_A")))
    first <- NULL
    for(r in seq_len(num_reps)) {
      pop <- .laugh_sim_pedigree(prep)
      counts[r,] <- count_laugh_kin(pop$pedigree,pop$samples,MAX_JUVE_AGE,beqy_POP)
      if(r==1 && keep_first) first <- pop
    }
    list(counts=data.frame(rep=seq_len(num_reps),counts,check.names=FALSE),
      preparation=prep,seed=seed,first=first)
  })
}

#' Count pedigree connections independently of microscoping
#' @param pedigree Data frame with sequential id, birth_year, sex, mother, father.
#' Unknown parents are zero; known parents must precede their offspring in id order.
#' @param samples Data frame with id, year, age; no repeated individual within year.
#' @param MAX_JUVE_AGE,beqy_POP Comparison conventions from prepare_laugh_simulation().
#' @return Six named counts in microscoping's order. HSPs count shared parents
#' (a full-sibling dyad contributes two); GGPs count grandparent paths. Both are
#' restricted to different birth cohorts and two juvenile sample events. A dyad
#' can contribute to multiple summaries. Self-recaptures contribute no kin paths.
#' @export
count_laugh_kin <- function(pedigree,samples,MAX_JUVE_AGE,beqy_POP=FALSE) {
  .laugh_sim_integer(MAX_JUVE_AGE,"MAX_JUVE_AGE",0)
  if(!is.logical(beqy_POP) || length(beqy_POP)!=1 || is.na(beqy_POP))
    stop("beqy_POP must be TRUE or FALSE.")
  p <- pedigree; s <- samples
  if(!is.data.frame(p) || !all(c("id","birth_year","sex","mother","father") %in% names(p)) ||
     !is.data.frame(s) || !all(c("id","year","age") %in% names(s)))
    stop("Supply pedigree and sample tables with the documented columns.")
  if(!is.numeric(p$id) || any(p$id!=seq_len(nrow(p))) ||
     any(!is.finite(as.matrix(p[c("id","birth_year","mother","father")]))))
    stop("Pedigree IDs must be sequential, with finite numeric ancestry.")
  if(any(!p$sex %in% c("F","M")) || any(p$birth_year!=floor(p$birth_year)))
    stop("Pedigree sex must be F/M and birth years must be integers.")
  for(k in c("mother","father")) {
    v <- p[[k]]
    if(any(v<0 | v>=p$id | v!=floor(v))) stop("Parents must be zero or preceding IDs.")
    known <- v>0
    if(any(p$sex[v[known]]!=if(k=="mother")"F" else "M") ||
       any(p$birth_year[v[known]]>=p$birth_year[known]))
      stop("Known parents must have the correct sex and be born before offspring.")
  }
  if(any(!is.finite(as.matrix(s[c("id","year","age")])) ) ||
     any(s$id<1 | s$id>nrow(p) | s$id!=floor(s$id)) ||
     any(s$age<0 | s$age!=floor(s$age) | s$year!=floor(s$year) |
         s$age!=s$year-p$birth_year[s$id]) ||
     anyDuplicated(s[c("id","year")])) stop("Invalid sample IDs, ages, or repeated samples within year.")
  id <- s$id; b <- p$birth_year[id]; m <- p$mother[id]; f <- p$father[id]
  gp <- cbind(c(0,p$mother)[m+1],c(0,p$father)[m+1],
    c(0,p$mother)[f+1],c(0,p$father)[f+1])
  middle_birth <- cbind(c(NA,p$birth_year)[m+1],c(NA,p$birth_year)[m+1],
    c(NA,p$birth_year)[f+1],c(NA,p$birth_year)[f+1])
  ans <- c(POP_offJ=0,POP_offA=0,HSP_JJ=0,GGP_JJ=0,
    nsamp_J=sum(s$age<=MAX_JUVE_AGE),nsamp_A=sum(s$age>MAX_JUVE_AGE))
  if(length(id)<2) return(ans)
  for(i in seq_len(length(id)-1L)) {
    j <- (i+1L):length(id)
    # Each orientation has its own offspring age and parent sampling year.
    ij <- (m[j]==id[i] | f[j]==id[i]) & s$year[i]>=b[j]+!beqy_POP
    ji <- (m[i]==id[j] | f[i]==id[j]) & s$year[j]>=b[i]+!beqy_POP
    ans["POP_offJ"] <- ans["POP_offJ"]+sum(ij & s$age[j]<=MAX_JUVE_AGE)+sum(ji & s$age[i]<=MAX_JUVE_AGE)
    ans["POP_offA"] <- ans["POP_offA"]+sum(ij & s$age[j]>MAX_JUVE_AGE)+sum(ji & s$age[i]>MAX_JUVE_AGE)
    eligible <- b[j]!=b[i] & s$age[i]<=MAX_JUVE_AGE & s$age[j]<=MAX_JUVE_AGE & id[j]!=id[i]
    h <- (m[i]>0 & m[j]==m[i])+(f[i]>0 & f[j]==f[i])
    g <- rowSums(gp[j,,drop=FALSE]==id[i] & s$year[i]>middle_birth[j,,drop=FALSE],na.rm=TRUE)+
      vapply(j,function(z)sum(gp[i,]==id[z] & s$year[z]>middle_birth[i,],na.rm=TRUE),numeric(1))
    ans["HSP_JJ"] <- ans["HSP_JJ"]+sum(h[eligible])
    ans["GGP_JJ"] <- ans["GGP_JJ"]+sum(g[eligible])
  }
  ans
}

#' Compare independent pedigree counts with the installed laugh test
#' @param design Named list of arguments for microscoping::ckmr_laugh_test().
#' @param num_reps,seed,keep_first,burn_in Simulation settings.
#' @param max_population,max_individuals,max_pair_comparisons,max_total_individuals,max_total_comparisons Simulation budgets.
#' @param max_work_cells Maximum cells in the package's principal age-cubed by
#' year-squared work grid. Checked before invoking its kin calculations.
#' @return List with comparison table, full simulation result, package version,
#' and canonical design. z_difference is (mean-expected)/Monte Carlo SE; when
#' SE is zero it is NA. There is deliberately no automatic pass/fail label:
#' rounding, tail approximations, and finite-population effects need inspection.
#' @export
check_laugh_test <- function(design,num_reps=100,seed=1401,keep_first=TRUE,
    burn_in=NULL,max_population=20000,max_individuals=200000,
    max_pair_comparisons=500000,max_total_individuals=2e7,
    max_total_comparisons=2e7,max_work_cells=5e6) {
  .laugh_sim_integer(num_reps,"num_reps",2)
  if(!is.list(design) || is.null(names(design)) || any(!nzchar(names(design))) ||
     anyDuplicated(names(design)) || !all(names(design) %in% names(formals(prepare_laugh_simulation))[1:12]))
    stop("design must be a named list of laugh-test arguments, without simulation controls.")
  # Check the analytic work grid before spending time on simulation or letting
  # the package allocate its probability arrays. This repeats only demographic
  # preparation; neither this preflight nor the simulator computes kin probabilities.
  .laugh_sim_integer(max_work_cells,"max_work_cells",1)
  preflight <- do.call(prepare_laugh_simulation,c(design,list(burn_in=burn_in,
    max_population=max_population,max_individuals=max_individuals,
    max_pair_comparisons=max_pair_comparisons)))
  if(nrow(preflight$ages)^3*preflight$num_years^2>max_work_cells)
    stop("Package work grid exceeds max_work_cells; reduce age span or sampling years.")
  sim <- do.call(simulate_laugh_test,c(design,list(num_reps=num_reps,seed=seed,
    keep_first=keep_first,burn_in=burn_in,max_population=max_population,
    max_individuals=max_individuals,max_pair_comparisons=max_pair_comparisons,
    max_total_individuals=max_total_individuals,max_total_comparisons=max_total_comparisons)))
  prep <- sim$preparation
  expected <- do.call(microscoping::ckmr_laugh_test,prep$arguments)
  keys <- names(sim$counts)[-1]
  if(!all(keys %in% names(expected)) || any(!is.finite(expected[keys])))
    stop("The installed laugh test returned unexpected or nonfinite outputs.")
  means <- colMeans(sim$counts[keys])
  mcse <- vapply(sim$counts[keys],stats::sd,numeric(1))/sqrt(num_reps)
  difference <- means-as.numeric(expected[keys])
  tab <- data.frame(metric=keys,expected=as.numeric(expected[keys]),
    simulated_mean=unname(means),mcse=unname(mcse),difference=unname(difference),
    z_difference=ifelse(mcse>0,difference/mcse,NA_real_),row.names=NULL)
  list(comparison=tab,simulation=sim,package_version=prep$package_version,
    design=prep$arguments)
}

#' Generate one complete pedigree under fixed age/sex census targets
#' @param prep Validated output of prepare_laugh_simulation().
#' @return Pedigree and sample-event data frames. IDs are row indices; all
#' generations are retained up to the preflighted record budget. Founders have
#' unknown parents. Uniform thinning is independent of ancestry and output.
#' death_year is the first census year an individual is absent; NA means still
#' alive at the final census. This allows reconstruction of every census.
#' @keywords internal
.laugh_sim_pedigree <- function(prep) {
  n <- prep$per_sex; ages <- prep$ages$age; k <- length(n)
  cap <- prep$records_per_rep
  birth <- mother <- father <- numeric(cap); sex <- character(cap)
  death <- rep(NA_real_,cap)
  alive <- lapply(1:2,function(s)lapply(n,function(x)integer(x)))
  used <- 0L; start <- 1L-prep$burn_in
  for(s in 1:2) for(a in seq_len(k)) if(n[a]>0) {
    ids <- used+seq_len(n[a]); used <- used+n[a]
    birth[ids] <- start-ages[a]; sex[ids] <- c("F","M")[s]
    alive[[s]][[a]] <- ids
  }
  observations <- list(); obs <- 0L
  for(y in seq.int(start,prep$num_years)) {
    if(y>start) {
      for(s in 1:2) {
        death[alive[[s]][[k]]] <- y
        for(a in k:2) {
          pool <- alive[[s]][[a-1L]]
          kept <- if(n[a]) pool[sample.int(length(pool),n[a])] else integer()
          death[setdiff(pool,kept)] <- y
          alive[[s]][[a]] <- kept
        }
      }
      # Newborns cannot be parents. Draw independently for every new individual.
      moms <- unlist(alive[[1]][-1],use.names=FALSE)
      dads <- unlist(alive[[2]][-1],use.names=FALSE)
      weights <- rep(prep$ages$relative_output[-1],n[-1])
      for(s in 1:2) {
        ids <- used+seq_len(n[1]); used <- used+n[1]
        birth[ids] <- y; sex[ids] <- c("F","M")[s]
        mother[ids] <- moms[sample.int(length(moms),n[1],replace=TRUE,prob=weights)]
        father[ids] <- dads[sample.int(length(dads),n[1],replace=TRUE,prob=weights)]
        alive[[s]][[1]] <- ids
      }
    }
    if(y>=1) for(s in 1:2) for(a in seq_len(k)) {
      expected <- prep$allocation[y,a]
      take <- floor(expected)+as.integer(runif(1)<expected-floor(expected))
      if(take>0) {
        ids <- alive[[s]][[a]][sample.int(n[a],take)]
        obs <- obs+1L
        observations[[obs]] <- data.frame(id=ids,year=y,age=ages[a],sex=c("F","M")[s])
      }
    }
  }
  stopifnot(used==cap)
  list(pedigree=data.frame(id=seq_len(used),birth_year=birth,death_year=death,
    sex=sex,mother=mother,father=father),
    samples=if(obs)do.call(rbind,observations) else data.frame(id=integer(),year=integer(),age=integer(),sex=character()))
}

#' Check a bounded integer control for the standalone simulator
#' @param x Value to check.
#' @param name Name used in the error message.
#' @param lower,upper Inclusive bounds.
#' @return Invisibly x, or an error before allocating a population.
#' @keywords internal
.laugh_sim_integer <- function(x,name,lower=0,upper=.Machine$integer.max) {
  if(!is.numeric(x) || length(x)!=1 || !is.finite(x) || x!=floor(x) || x<lower || x>upper)
    stop(name," must be a whole number from ",lower," to ",upper,".")
  invisible(x)
}

#' Evaluate simulation with an isolated, reproducible RNG state
#' @param seed Validated positive integer seed.
#' @param code Expression to evaluate after setting the simulation RNG kinds.
#' @return Value of code. Restores caller state and kinds, including on error.
#' @keywords internal
.laugh_sim_seed <- function(seed,code) {
  kinds <- RNGkind(); had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had) old <- get(".Random.seed",.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(kinds))
    if(had) assign(".Random.seed",old,.GlobalEnv) else
      if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)})
  set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
  force(code)
}
