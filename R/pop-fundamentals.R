
library(MASS) # Used for generating multivariate normal distributions

simulate_vb_growth <- function(n_fish, ages, mean_Linf = 80, mean_K = 0.2, t0 = -0.5,
                               sd_Linf = 5, sd_K = 0.02, corr_Linf_K = -0.8, sd_obs = 2) {

  # 1. Set up the variance-covariance matrix for Linf and K
  cov_matrix <- matrix(c(
    sd_Linf^2,                         sd_Linf * sd_K * corr_Linf_K,
    sd_Linf * sd_K * corr_Linf_K,     sd_K^2
  ), nrow = 2)

  # 2. Draw individual growth traits
  indiv_traits <- mvrnorm(n_fish, mu = c(mean_Linf, mean_K), Sigma = cov_matrix)
  indiv_Linf <- indiv_traits[, 1]
  indiv_K    <- indiv_traits[, 2]

  # Ensure traits remain biologically realistic (strictly positive)
  indiv_Linf <- pmax(indiv_Linf, 1)
  indiv_K    <- pmax(indiv_K, 0.01)

  # 3. Simulate growth over time with observation error
  sim_data <- data.frame(fish_id = 1:n_fish, age = ages)

  # Calculate deterministic true individual length
  sim_data$true_length <- indiv_Linf[sim_data$fish_id] *
    (1 - exp(-indiv_K[sim_data$fish_id] * (sim_data$age - t0)))

  # Add measurement/observation noise
  sim_data$obs_length <- sim_data$true_length + rnorm(nrow(sim_data), mean = 0, sd = sd_obs)

  return(sim_data)
}

# Example usage: Simulate 500 fish across random ages between 1 and 15
set.seed(42)
max_age <- 20
sampled_ages <- rexp(500, 0.2)
sampled_ages <- sampled_ages[sampled_ages <= 20]
growth_dataframe <- simulate_vb_growth(n_fish = length(sampled_ages), ages = sampled_ages)

fecundity_exponent <- 2.7
as_tibble(growth_dataframe) %>% ggplot(aes(x = age, y = true_length^fecundity_exponent)) + geom_point()

as_tibble(growth_dataframe) %>% ggplot(aes(x = true_length^fecundity_exponent)) + geom_histogram()





#' Assign ages, lengths and fecundity to a set of individuals
#'
#' Super quick and dirty
#' @param ann_surv annual survival prob.
#' @param vBpars von Bertalanffy parameters



#' Function to simulate super simple POP scenarios
#'
#' @param Nf number of females in pop
#' @param Nm number of males in pop
#' @param Noffs total number of offspring simulated in pop
#' @param Sf number of females sampled
#' @param Sm number of males sampled
#' @param Soffs number of offspring sampled
simple_pedigree <- function(
    Nf=10000, Nm=10000, Noffs=20000, Sf = 1000, Sm = 1000, Soffs = 1000,
    fem_cv = 1
  ) {
  moms <- 1:Nf
  dads <- 1:Nm
  kids <- 1:Noffs

  # set the weights to just be the same now
  mom_wts <- rbeta(length(moms), shape1 = 0.01, shape2 = 1)
  dad_wts <- rbeta(length(dads), shape1 = 0.01, shape2 = 1)

  # super simple random mating
  pedigree <- tibble(
    kid = kids
  ) %>%
    mutate(
      pa = sample(moms, size = n(), prob = dad_wts, replace = TRUE),
      ma = sample(dads, size = n(), prob = dad_wts, replace = TRUE),
      pa_wts = dad_wts[pa],
      ma_wts = mom_wts[ma]
    )

  # sample the kids by simply sampling the pedigree
  kid_samp <- pedigree %>%
    sample_n(Soffs)

  # then sample the males and females
  mom_samp <- sample(moms, size = Sf)
  dad_samp <- sample(dads, size = Sm)

  # extract the MOPs and FOPs
  MOPs <- kid_samp %>%
    filter(ma %in% mom_samp)
  POPs <- kid_samp %>%
    filter(pa %in% dad_samp)
  MOPs
}


boing <- lapply(1:100, function(x) simple_pedigree()) %>% bind_rows(.id = "rep")
boing_summ <- boing %>% count(rep)
ggplot(boing_summ, aes(x = n)) +
  geom_histogram(binwidth = 2) +
  geom_vline(xintercept = mean(boing_summ$n), colour = "red")



ggplot(boing, aes(x = ma_wts)) + geom_histogram()
