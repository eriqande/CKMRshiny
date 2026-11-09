# Sanctuary Sound helpers, ordered as they are introduced in fit_delfi_S.Rmd.
# Data preparation stays visible in the notebook. These functions handle the
# calculations that are repeated during calibration or the worked age mixture.

#' Describe left-closed length bins, including both open-ended tails
#' @param cuts Strictly increasing finite length boundaries.
#' @return Tibble with bin_id, lower, upper, and a readable length_bin label.
#' @export
delfi_s_length_bins <- function(cuts) {
  stopifnot(is.numeric(cuts), length(cuts) > 0, all(is.finite(cuts)), all(diff(cuts) > 0))
  lower <- c(-Inf, cuts)
  upper <- c(cuts, Inf)
  label <- function(x) formatC(x, format = "f", digits = 2)
  tibble::tibble(bin_id = seq_along(lower), lower = lower, upper = upper,
    length_bin = ifelse(is.infinite(lower), paste0("<", label(upper)),
      ifelse(is.infinite(upper), paste0(">=", label(lower)),
        paste0("[", label(lower), ", ", label(upper), ")"))))
}

#' Construct a normalized exponential age prior on a finite age grid
#' @param log_decline Log of the positive age-abundance decline parameter.
#' This is an age-prior parameter, not a direct estimate of annual survival.
#' @param ages Unique positive integer ages supporting the calibration.
#' @return Tibble with age and p_age, summing to one across the supplied ages.
#' @export
delfi_s_age_prior <- function(log_decline, ages) {
  stopifnot(length(log_decline) == 1L, is.finite(log_decline),
    length(ages) > 0, !anyDuplicated(ages), all(is.finite(ages)),
    all(ages >= 1 & ages == floor(ages)))
  # Subtracting the youngest age cancels in normalization and avoids underflow
  # of every weight if a caller explores a steep decline.
  w <- exp(-exp(log_decline) * (ages - min(ages)))
  tibble::tibble(age = ages, p_age = w / sum(w))
}

#' Apply Bayes' rule to a discretized length-given-age distribution
#' @param log_decline Log decline parameter for one epoch's age prior.
#' @param length_probs Tibble with one row per age and length bin, including
#' age, bin_id, and p_length_given_age. Other descriptive columns are retained.
#' Every bin must contain the same age grid; probabilities sum to one over bins
#' for each age. These are fixed calibration inputs, not estimated here.
#' @return Input rows plus p_age, joint_probability, p_length_bin, and
#' p_age_given_length. The last column sums to one within each length bin.
#' @export
delfi_s_age_key <- function(log_decline, length_probs) {
  prior <- delfi_s_age_prior(log_decline, sort(unique(length_probs$age)))
  length_probs |>
    dplyr::left_join(prior, by = "age", relationship = "many-to-one") |>
    dplyr::mutate(joint_probability = p_length_given_age * p_age) |>
    dplyr::group_by(bin_id) |>
    dplyr::mutate(p_length_bin = sum(joint_probability),
                  p_age_given_length = joint_probability / p_length_bin) |>
    dplyr::ungroup()
}

#' Evaluate the conditional known-age calibration log likelihood for one epoch
#' @param log_decline Log age-prior decline parameter.
#' @param aged_counts Tibble with age, bin_id, and n_aged for one epoch.
#' Counts may include explicit zeros. Unknown ages must be removed beforehand.
#' @param length_probs Fixed length-given-age probabilities for the full age grid.
#' @return Sum of n_aged times log P(age | length bin), omitting only zero-count
#' terms. Positive observations with zero probability return negative infinity.
#' @export
delfi_s_calibration_loglik <- function(log_decline, aged_counts, length_probs) {
  key <- delfi_s_age_key(log_decline, length_probs)
  rows <- aged_counts |>
    dplyr::filter(n_aged > 0) |>
    dplyr::left_join(dplyr::select(key, age, bin_id, p_age_given_length),
                    by = c("age", "bin_id"), relationship = "one-to-one")
  if (anyNA(rows$p_age_given_length)) stop("Observed age/bin is absent from the calibration grid.")
  sum(rows$n_aged * log(rows$p_age_given_length))
}

#' Expand one observed comparison category over both uncertain ages
#'
#' This transparent reference calculation is for one worked example, not the
#' optimized full likelihood. Both relationship directions are included. A
#' possible mother's maturity is evaluated at the possible daughter's birth;
#' survival is required only for births after the mother's sampling year.
#' The two age probabilities are mixed independently, as in Mark's model.
#' @param category One-row tibble with L1, Y1, L2, and Y2; L columns are bin IDs.
#' @param age_key_by_year Tibble with year, bin_id, age, and p_age_given_length.
#' @param log_N0 Log adult female abundance at the reference year y0.
#' @param growth_rate Exponential annual growth rate of adult female abundance.
#' @param annual_survival Constant annual survival probability, between zero and one.
#' @param y0 Reference year for the illustrative abundance trajectory.
#' @param maturity_age Knife-edge age at first reproduction.
#' @return One row per age pair, with sampling years, possible birth years, age-mixture weights,
#' conditional probabilities, and weighted contributions for both directions.
#' @export
delfi_s_age_pair_example <- function(category, age_key_by_year, log_N0,
                                     growth_rate, annual_survival, y0, maturity_age) {
  stopifnot(nrow(category) == 1, annual_survival > 0, annual_survival <= 1,
    is.finite(log_N0), is.finite(growth_rate), maturity_age > 0)
  first <- age_key_by_year |>
    dplyr::filter(year == category$Y1, bin_id == category$L1) |>
    dplyr::transmute(age_1 = age, p_age_1 = p_age_given_length)
  second <- age_key_by_year |>
    dplyr::filter(year == category$Y2, bin_id == category$L2) |>
    dplyr::transmute(age_2 = age, p_age_2 = p_age_given_length)
  stopifnot(nrow(first) > 0, nrow(second) > 0,
    abs(sum(first$p_age_1) - 1) < 1e-10, abs(sum(second$p_age_2) - 1) < 1e-10)
  tidyr::crossing(first, second) |>
    dplyr::mutate(
      sampling_year_1 = category$Y1,
      sampling_year_2 = category$Y2,
      birth_1 = sampling_year_1 - age_1,
      birth_2 = sampling_year_2 - age_2,
      joint_age_probability = p_age_1 * p_age_2,
      mature_2_at_birth_1 = birth_1 >= birth_2 + maturity_age,
      mature_1_at_birth_2 = birth_2 >= birth_1 + maturity_age,
      future_years_2_to_1 = pmax(0, birth_1 - sampling_year_2),
      future_years_1_to_2 = pmax(0, birth_2 - sampling_year_1),
      N_at_birth_1 = exp(log_N0 + growth_rate * (birth_1 - y0)),
      N_at_birth_2 = exp(log_N0 + growth_rate * (birth_2 - y0)),
      p_2_mother_of_1 = mature_2_at_birth_1 * annual_survival^future_years_2_to_1 / N_at_birth_1,
      p_1_mother_of_2 = mature_1_at_birth_2 * annual_survival^future_years_1_to_2 / N_at_birth_2,
      contribution_2_to_1 = joint_age_probability * p_2_mother_of_1,
      contribution_1_to_2 = joint_age_probability * p_1_mother_of_2,
      contribution = contribution_2_to_1 + contribution_1_to_2
    )
}
