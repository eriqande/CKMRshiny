# Small analytic checks independent of Mark's packages or course datasets.
# Run from the repository root: Rscript tests/check-tidy-delfi-S.R
source("R/tidy-delfi-S.R")

bins <- delfi_s_length_bins(c(1.5, 1.8))
stopifnot(identical(bins$lower, c(-Inf, 1.5, 1.8)),
          identical(bins$upper, c(1.5, 1.8, Inf)),
          identical(findInterval(c(1.4, 1.5, 1.8, 2), c(-Inf, 1.5, 1.8, Inf)), c(1L,2L,3L,3L)))

# Bayes' rule on two ages/two bins has a direct analytic answer.
forward <- tibble::tibble(age = c(2,7,2,7), bin_id = c(1,1,2,2),
                         p_length_given_age = c(.9,.1,.1,.9))
key <- delfi_s_age_key(log(.2), forward)
prior_odds <- exp(-.2 * (7-2))
p_young_short <- .9 / (.9 + .1 * prior_odds)
p_old_long <- .9 * prior_odds / (.1 + .9 * prior_odds)
stopifnot(abs(key$p_age_given_length[1] - p_young_short) < 1e-12,
          abs(key$p_age_given_length[4] - p_old_long) < 1e-12,
          abs(sum(key$p_age_given_length[key$bin_id==1])-1) < 1e-12,
          abs(sum(key$p_age_given_length[key$bin_id==2])-1) < 1e-12)
counts <- tibble::tibble(age=c(2,7,2,7),bin_id=c(1,1,2,2),n_aged=c(3L,0L,0L,1L))
ll <- delfi_s_calibration_loglik(log(.2),counts,forward)
stopifnot(abs(ll - (3*log(p_young_short)+log(p_old_long))) < 1e-12)

# Known ages: female 2 was captured in 2005 at age 4, before a possible
# daughter's 2008 birth. She must survive three years and be mature at birth.
category <- tibble::tibble(L1=1L,Y1=2010L,L2=2L,Y2=2005L)
known <- tibble::tibble(year=c(2010L,2005L),bin_id=c(1L,2L),age=c(2L,4L),
                       p_age_given_length=c(1,1))
x <- delfi_s_age_pair_example(category,known,log(1000),.02,.8,2005,5)
expected <- .8^3 / (1000*exp(.02*3))
stopifnot(abs(x$contribution-expected) < 1e-12,
          x$contribution_1_to_2==0,x$future_years_2_to_1==3,x$mature_2_at_birth_1)

# A capture after that birth supplies evidence of survival up to birth;
# no prospective survival factor belongs in this direction.
category$Y2 <- 2010L
known$year[2] <- 2010L
known$age[2] <- 9L
x <- delfi_s_age_pair_example(category,known,log(1000),.02,.2,2005,5)
y <- delfi_s_age_pair_example(category,known,log(1000),.02,.9,2005,5)
stopifnot(identical(x$contribution,y$contribution),
          abs(x$contribution - 1/(1000*exp(.02*3))) < 1e-12,
          x$future_years_2_to_1==0)

# The reversed relationship is a separate, mutually exclusive direction.
reverse <- tibble::tibble(L1=2L,Y1=2010L,L2=1L,Y2=2010L)
r <- delfi_s_age_pair_example(reverse,known,log(1000),.02,.9,2005,5)
stopifnot(abs(r$contribution-x$contribution) < 1e-12,
          r$contribution_2_to_1==0,r$contribution_1_to_2>0)
cat("Sanctuary bins, Bayes calibration, maturity, both directions, and survival timing passed.\n")
