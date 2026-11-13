# Editor input preparation; the package adapter is in laugh-calculation.R.

#' Create illustrative age schedules for the laugh-test input editor
#' @param plus_age Final input age, representing all animals at that age or older.
#' @param abundance Total abundance, both sexes combined, including the plus group.
#' @param num_years Initial number of year-specific sampling curves (1-12).
#' @return Editor state with three population curves and a sample_scores matrix.
#' Rows beyond the current sampling-year count can be retained as hidden history.
#' Numbers follow illustrative annual survival 0.75; they are not an assessment.
#' @export
laugh_default_curves <- function(plus_age = 12, abundance = 1e5, num_years = 3) {
  .pop_integer(plus_age, "Plus-group age")
  if (plus_age < 4 || plus_age > 60) stop("Choose a plus-group age from 4 to 60.")
  if (length(abundance) != 1L || !is.finite(abundance) || abundance <= 0)
    stop("Total abundance must be positive.")
  .pop_integer(num_years,"Number of sampling curves")
  if (num_years > 12) stop("Use at most 12 sampling years.")
  a <- 0:plus_age
  n <- c(.75^(0:(plus_age-1)), .75^plus_age / (1-.75))
  maturity <- plogis((a - min(3.5, plus_age / 2)) / .8)
  maturity[1] <- 0
  scores <- exp(-((a-1.3)/1.1)^2) + .35 * exp(-((a-6)/2.5)^2)
  list(plus_age = as.integer(plus_age), curves = list(
    nata = abundance * n / sum(n),
    wata = .05 + 8 * (-expm1(-.23 * (a + .2)))^3,
    pmatata = maturity), sample_scores = matrix(rep(scores,num_years),
      nrow=num_years,byrow=TRUE,dimnames=list(paste("Year",seq_len(num_years)),
        c(as.character(0:(plus_age-1)),paste0(plus_age,"+")))))
}

#' Retain cached sampling curves and initialize genuinely new years
#' @param state Editor state with a sample_scores matrix.
#' @param num_years Number of consecutive years currently requested (1-12).
#' @param copy_year Existing year whose scores initialize previously unseen years.
#' @return State with enough rows. Reducing num_years never deletes rows; restoring
#' a year recovers its prior values. New rows are independent copies of copy_year.
#' @export
laugh_ensure_years <- function(state,num_years,copy_year=1) {
  .pop_integer(num_years,"Number of sampling years")
  if (num_years>12) stop("Use at most 12 sampling years.")
  .pop_integer(copy_year,"Template year")
  scores <- state$sample_scores
  stopifnot(is.matrix(scores),nrow(scores)>=copy_year,ncol(scores)==state$plus_age+1)
  if (num_years>nrow(scores)) {
    added <- matrix(rep(scores[copy_year,],num_years-nrow(scores)),
                    nrow=num_years-nrow(scores),byrow=TRUE)
    state$sample_scores <- rbind(scores,added)
    rownames(state$sample_scores) <- paste("Year",seq_len(num_years))
  }
  state
}

#' Stable display colors for consecutive sampling years
#' @return Twelve year colors, independent of which years are visible or selected.
#' @keywords internal
.laugh_year_colors <- function() c("#3c6ea8","#bd7330","#237e83","#8065a6",
  "#ba4d65","#748529","#458e9b","#8f6346","#616cb2","#9c698f","#7a8b96","#b08434")

#' Describe the four raw editor schedules and their biological meanings
#' @return Named list with labels, colors, units, explanations, and consequences.
#' The sampling editor stores allocation scores; only their normalized values
#' are exported as samp_pata. This distinction is deliberate and visible.
#' @keywords internal
.laugh_curve_info <- function() {
  list(
    nata = list(title = "Population numbers at age", color = "#237e83", unit = "Individuals (both sexes)",
      definition = "How many animals of each age are in the population, with females and males combined. These are population numbers, not numbers sampled. The last entry is the total for every age in the plus group.",
      consequence = "Absolute scale matters: more animals generally mean fewer kin matches for a fixed sampling effort. Draw a shape, then use Set total abundance to scale that shape to a chosen total.",
      package = "Mark's function smooths these numbers to describe a steady-state population and infer survival. Editor values remain your supplied inputs; the calculation uses a smoothed copy, shown with its results. The example decline is illustrative, not a stock assessment."),
    wata = list(title = "Biological weight at age", color = "#8065a6", unit = "Body weight (e.g. kg)",
      definition = "The typical body weight of an animal at each age, in consistent units such as kilograms. This is biological body weight, not a reproductive sampling weight or a sampling probability.",
      consequence = "The function assumes reproductive output is proportional to body weight times proportion mature. Changing relative weights across ages matters. Multiplying every age's weight by the same constant changes units, but leaves kin probabilities unchanged.",
      package = "For ages in the expanded plus group, the function keeps the final supplied weight constant. This input does not specify individual variation in weight at a given age."),
    pmatata = list(title = "Proportion mature at age", color = "#bd7330", unit = "Proportion mature (0 to 1)",
      definition = "The fraction of animals of each age that are reproductively mature: 0 means none and 1 means all. Together, biological weight and maturity determine expected reproductive output at age.",
      consequence = "This is distinct from the juvenile comparison cutoff. An age can contribute a little reproduction and still be included in the juvenile comparison category; the cutoff also reflects how accurately offspring can be aged.",
      package = "The final supplied maturity value is held constant across the expanded plus group. This one-sex-symmetric function applies the same schedule to males and females."),
    samp_pata = list(title = "Age composition of the sample", color = "#3c6ea8", unit = "Relative allocation score",
      definition = "Choose a sampling year, then draw its relative allocation scores across ages. Each year's scores are divided by their own sum to obtain that year's samp_pata proportions. You do not need to make the drawn values add to one.",
      consequence = "This is NOT the fraction of the population caught at an age. A normalized value of 20% means 20% of that year's samples come from that age group. The annual count determines the amount of sampling; the year's curve determines its age mix. Doubling every score in a row leaves that mix unchanged.",
      package = "The editor's samp_pata matrix contains proportions, with one row per active year and one column per age category. Each row sums to one, including gap years. Before calling microscoping we multiply each row by its annual count: the package's samp_nya matrix contains sample NUMBERS, with the first column corresponding to age 0. Its final column allocates samples to the whole plus group, which the package splits uniformly over its expanded older ages.")
  )
}

#' Generate a simple starting shape on the current input age grid
#' @param state Editor state from laugh_default_curves().
#' @param key One of nata, wata, pmatata, samp_pata.
#' @param preset One of example, flat, increasing, declining, bell, clear.
#' @param year Sampling year to use when choosing a sampling-score scale.
#' @return Nonnegative numeric vector. Presets retain the current maximum scale
#' (or an illustrative scale if empty), and preserve total population abundance.
#' @export
laugh_curve_preset <- function(state, key, preset = "example",year=1) {
  stopifnot(key %in% names(.laugh_curve_info()))
  if (preset == "example") {
    defaults <- laugh_default_curves(state$plus_age,max(1,sum(state$curves$nata)))
    return(if(key=="samp_pata") unname(defaults$sample_scores[1,]) else defaults$curves[[key]])
  }
  x <- seq(0, 1, length.out = state$plus_age + 1)
  shape <- switch(preset, flat = rep(1, length(x)), increasing = x,
    declining = 1-x, bell = exp(-((x-.4)/.22)^2), clear = rep(0, length(x)),
    stop("Unknown starting shape."))
  current <- if(key=="samp_pata") state$sample_scores[year,] else state$curves[[key]]
  amplitude <- if (key == "pmatata") 1 else max(current, 1)
  result <- shape * amplitude
  if (key == "nata" && sum(result) > 0)
    result <- result * max(1, sum(state$curves$nata)) / sum(result)
  result
}

#' Parse consecutive sampling-year counts, including zero-sample gap years
#' @param text Comma, whitespace, or semicolon-separated nonnegative integers.
#' @return Numeric vector of 1-12 annual counts, preserving zero entries.
#' @export
laugh_parse_years <- function(text) {
  if (length(text) != 1 || is.na(text) || !nzchar(trimws(text)))
    stop("Enter sample counts for at least one year.")
  tokens <- strsplit(trimws(text), "[,;[:space:]]+")[[1]]
  values <- suppressWarnings(as.numeric(tokens))
  if (!length(values) || length(values) > 12 || any(!is.finite(values)) ||
      any(values < 0 | values != floor(values)))
    stop("Enter 1 to 12 nonnegative whole-number annual counts, separated by commas or spaces.")
  if (sum(values) <= 0) stop("At least one year must have a positive sample count.")
  values
}

#' Check a raw curve edit without requiring a completed sampling design
#' @param state Current editor state.
#' @param key Curve name.
#' @param values Proposed values in age order 0 through the plus-group age.
#' @return Numeric vector with valid length and bounds. All-zero curves are
#' allowed while editing; readiness validation is handled separately.
#' @export
laugh_validate_curve <- function(state, key, values) {
  if (!key %in% names(.laugh_curve_info())) stop("Unknown age curve.")
  if (!is.numeric(values) || length(values) != state$plus_age + 1 || any(!is.finite(values)) || any(values < 0))
    stop("Every age needs a finite, nonnegative value.")
  if (key == "pmatata" && any(values > 1)) stop("Maturity proportions must lie between zero and one.")
  as.numeric(values)
}

#' Prepare normalized editor inputs for the laugh-test design
#' @param state Raw editor curves.
#' @param years Annual sample counts from laugh_parse_years().
#' @param juvenile Maximum exact age classified as juvenile for comparisons.
#' @param catch_min Youngest age used for fitting the steady-state catch curve.
#' @param after_spawning Allow a parent caught in the offspring's birth year.
#' @return Named argument list, with row-normalized year-by-age samp_pata matrix
#' and explicit package
#' defaults. No package call, age smoothing, plus-group expansion, or kin-count
#' computation is performed. Convert with laugh_function_arguments() before
#' calling microscoping: its samp_nya argument requires counts, not proportions.
#' @export
laugh_prepare_inputs <- function(state, years, juvenile = 2, catch_min = 1, after_spawning = FALSE) {
  keys <- c("nata","wata","pmatata")
  curves <- lapply(keys, function(key)
    laugh_validate_curve(state, key, state$curves[[key]]))
  names(curves) <- keys
  if (sum(curves$nata) <= 0) stop("Population numbers must have a positive total.")
  if (sum(curves$nata * curves$wata * curves$pmatata) <= 0)
    stop("At least one populated age must have positive weight and maturity.")
  if (length(juvenile) != 1 || !is.finite(juvenile) || juvenile < 0 || juvenile != floor(juvenile) || juvenile >= state$plus_age)
    stop("Choose a whole-number juvenile cutoff below the plus-group age.")
  if (length(catch_min) != 1 || !is.finite(catch_min) || catch_min < 0 || catch_min != floor(catch_min) || catch_min > state$plus_age-3)
    stop("The catch-curve start must leave at least three exact-age classes to fit.")
  if (!is.numeric(years) || !length(years) || any(!is.finite(years)) || any(years < 0 | years != floor(years)) || sum(years) <= 0)
    stop("Supply nonnegative annual sample counts with a positive total.")
  if (!is.logical(after_spawning) || length(after_spawning) != 1 || is.na(after_spawning)) stop("Choose a sampling-timing option.")
  scores <- laugh_ensure_years(state,length(years))$sample_scores[seq_along(years),,drop=FALSE]
  if(any(!is.finite(scores)) || any(scores<0)) stop("Sampling scores must be finite and nonnegative.")
  empty <- which(rowSums(scores)<=0)
  if(length(empty)) stop("Draw a positive sampling score for year ",paste(empty,collapse=", ")," (including zero-sample years).")
  curves$samp_pata <- scores / rowSums(scores)
  c(curves, list(samp_ny = years, MAX_JUVE_AGE = as.integer(juvenile),
    MIN_CATCH_CURVE_AGE = as.integer(catch_min), beqy_POP = after_spawning,
    etriv = .01, extraplus = FALSE, want_gory_details = FALSE))
}

#' Preview expected sample allocation by age and sampling year
#' @param proportions Row-normalized year-by-age sample composition matrix.
#' A vector is also accepted and applied unchanged across all years.
#' @param years Nonnegative annual sample counts.
#' @return Tibble of age, year, and expected_samples. Values can be fractional:
#' these are allocations, not a stochastic or rounded sample realization.
#' @export
laugh_sample_allocation <- function(proportions, years) {
  if(!is.matrix(proportions)) proportions <- matrix(rep(proportions,length(years)),nrow=length(years),byrow=TRUE)
  stopifnot(nrow(proportions)==length(years), all(is.finite(proportions)), all(proportions >= 0), all(abs(rowSums(proportions)-1) < 1e-10),
    all(is.finite(years)), all(years >= 0))
  grid <- expand.grid(age = seq_len(ncol(proportions))-1L, year = seq_along(years))
  grid$expected_samples <- proportions[cbind(grid$year,grid$age+1L)] * years[grid$year]
  tibble::as_tibble(grid)
}
