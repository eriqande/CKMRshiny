#' Plot sibling population or sample family sizes by sex and episode
#' @param replicate Saved sibling replicate.
#' @param sampled Show sampled offspring counts, using all census adults as the
#' denominator. These family assignments are simulation truth, not adult surveys.
#' @param y_scale Linear or log10 count axis.
#' @return A ggplot with the same colors and histogram conventions as POP.
#' @export
hsp_plot_reproduction <- function(replicate, sampled = FALSE, y_scale = c("linear", "log10")) {
  y_scale <- match.arg(y_scale)
  x <- if (sampled) replicate$sampling_summary else replicate$reproduction_summary
  x$episode <- paste("Episode", x$episode)
  width <- max(1, ceiling((max(x$offspring) + 1) / 60))
  x$bin <- floor((x$offspring + .5) / width)
  bins <- stats::aggregate(n_parents ~ sex + episode + bin, x, sum)
  bins$lower <- bins$bin * width - .5; bins$upper <- bins$lower + width
  base <- if (y_scale == "log10") .5 else 0
  plot <- ggplot2::ggplot(bins) +
    ggplot2::geom_rect(ggplot2::aes(xmin = lower, xmax = upper, ymin = base, ymax = n_parents, fill = sex), color = "white") +
    ggplot2::facet_grid(episode ~ sex, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::scale_x_continuous(breaks = .pop_count_breaks) +
    ggplot2::labs(x = if (sampled) "Sampled offspring per parent" else "Total offspring per parent",
      y = "Number of parents", caption = paste("All census adults are included, including those with zero offspring.",
        if (sampled) "Family assignments use simulation truth; no adults are survey-sampled." else "Each row is a separate breeding episode.",
        if (y_scale == "log10") "Log10 bar baseline: 0.5." else "", sep = "\n")) + .pop_plot_theme()
  plot + if (y_scale == "log10") ggplot2::scale_y_log10(breaks = .pop_log_count_breaks) else
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks)
}

#' Compare persistent propensities or realized offspring counts across episodes
#' @param replicate Saved cross-cohort sibling replicate.
#' @param realized Plot realized offspring counts rather than mean-one sampling weights.
#' @param lengths Plot episode-specific lengths instead (size scenario only).
#' @return Sex-faceted scatter plot of a fixed display subset of at most 500 per sex.
#' @export
hsp_plot_persistence <- function(replicate, realized = FALSE, lengths = FALSE) {
  x <- replicate$display_adults
  x$x <- x[[if (realized) "offspring_1" else "weight_1"]]
  x$y <- x[[if (realized) "offspring_2" else "weight_2"]]
  if (lengths) {
    x$x <- x$length
    x$y <- if (is.null(x$length_2)) x$length else x$length_2
  }
  ggplot2::ggplot(x, ggplot2::aes(x, y, color = sex)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
    ggplot2::geom_point(alpha = .4, size = 1.4) +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free") +
    ggplot2::scale_color_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::labs(x = paste("Episode 1", if (lengths) "length" else if (realized) "offspring count" else "ERRO"),
      y = paste("Episode 2", if (lengths) "length" else if (realized) "offspring count" else "ERRO"),
      caption = if (lengths) "The same adults appear on both axes; display subset only.\nDashed: unchanged length. Growth cannot move a point below this line." else
        "The same adults appear on both axes; display subset only.\nDashed: equal values in both episodes. Relative ERRO can fall even when length increases.") + .pop_plot_theme()
}

#' Select saved episode-specific traits for population diagnostic views
#' @param replicate Size-based sibling replicate.
#' @param second Show episode two; legacy fixed-size archives fall back to episode one.
#' @return List of display adults, population histogram summary, and growth model.
#' @keywords internal
.hsp_episode_traits <- function(replicate, second = FALSE) {
  adults <- replicate$display_adults
  model <- replicate$settings$size_model
  summary <- replicate$size_summary
  if (second && !is.null(replicate$size_summary_2)) {
    adults$length <- adults$length_2; adults$age <- adults$age_2
    adults$reproductive_weight <- adults$weight_2
    adults$expected_reproductive_weight <- adults$weight_2
    summary <- replicate$size_summary_2
    if (isTRUE(replicate$settings$growth)) model <- lapply(model, function(p) {
      p$min_age <- p$min_age + 1; p$max_age <- p$max_age + 1; p
    })
  }
  list(adults = adults, summary = summary, model = model)
}

#' Plot trial-level sibling match counts with naive and conditional expectations
#' @param trial Saved sibling trial.
#' @param rep_id Selected replicate.
#' @return Faceted distribution; red is C/N, dotted is mean conditional expectation
#' over simulated populations, and dashed is the selected replicate's count.
#' @export
hsp_plot_pairs <- function(trial, rep_id = 1) {
  x <- trial$summary[trial$summary$method == "Moments", ]
  expected <- data.frame(sex = c("Female", "Male"),
    naive = trial$reps[[1]]$comparisons / c(trial$settings$Nf, trial$settings$Nm),
    conditional = rowMeans(vapply(trial$reps, function(r) r$effective$expected_sample_matches, numeric(2))))
  chosen <- x[x$rep_id == rep_id, ]
  width <- max(1, ceiling((diff(range(x$matches)) + 1) / 40))
  ggplot2::ggplot(x, ggplot2::aes(matches, fill = sex)) +
    ggplot2::geom_histogram(binwidth = width, boundary = -.5, color = "white") +
    ggplot2::geom_vline(data = expected, ggplot2::aes(xintercept = naive), color = "#bc423c", linewidth = 1) +
    ggplot2::geom_vline(data = expected, ggplot2::aes(xintercept = conditional), linetype = "dotted", linewidth = .9) +
    ggplot2::geom_vline(data = chosen, ggplot2::aes(xintercept = matches), linetype = "dashed", linewidth = .8) +
    ggplot2::facet_wrap(~sex, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::scale_x_continuous(breaks = .pop_count_breaks) +
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks) +
    ggplot2::labs(x = "Shared-parent matches in the offspring samples", y = "Number of replicates",
      caption = "Red: equal-reproduction expectation C/N. Dotted: mean expectation given simulated sampling weights.\nDashed: selected replicate. Full siblings contribute one match in each sex panel.") + .pop_plot_theme()
}

#' Compare population traits with offspring or sibling-parent representations
#' @param replicate Size-based sibling replicate retaining compact trait summaries.
#' @param matched Use one contribution per sibling match, otherwise per sampled offspring.
#' @param variable Compare reproductive_weight (ERRO) or length.
#' @return Normalized histograms and a dashed reproductive-weight expectation.
#' Cross-cohort sibling parents are weighted by the product of their two
#' episode-specific ERROs. Trait bins refer to episode one throughout.
#' @export
hsp_plot_traits <- function(replicate, matched = TRUE, variable = c("reproductive_weight", "length")) {
  variable <- match.arg(variable)
  x <- replicate$size_summary[replicate$size_summary$variable == variable, ]
  observed <- if (matched) "match_sum" else "sample_sum"
  expected <- if (matched) "sibling_weight_sum" else
    if ("sample_weight_sum" %in% names(x)) "sample_weight_sum" else "erro_sum"
  x$population_density <- x$observed_density <- x$expected_density <- 0
  for (sex in c("Female", "Male")) {
    rows <- which(x$sex == sex); width <- x$upper[rows] - x$lower[rows]
    x$population_density[rows] <- x$count[rows] / (sum(x$count[rows]) * width)
    if (sum(x[[observed]][rows]) > 0) x$observed_density[rows] <- x[[observed]][rows] / (sum(x[[observed]][rows]) * width)
    x$expected_density[rows] <- x[[expected]][rows] / (sum(x[[expected]][rows]) * width)
  }
  # Build an outline from bin edges, keeping exact histogram normalization.
  lines <- do.call(rbind, lapply(c("Female", "Male"), function(sex) {
    a <- x[x$sex == sex, ]
    data.frame(sex = sex, value = as.vector(rbind(a$lower, a$upper)),
      observed = rep(a$observed_density, each = 2), expected = rep(a$expected_density, each = 2))
  }))
  ggplot2::ggplot(x) +
    ggplot2::geom_rect(ggplot2::aes(xmin = lower, xmax = upper, ymin = 0, ymax = population_density, fill = sex), alpha = .2) +
    ggplot2::geom_line(data = lines, ggplot2::aes(value, observed, color = sex), linewidth = .9) +
    ggplot2::geom_line(data = lines, ggplot2::aes(value, expected), linetype = "dashed", color = "#233a3d", linewidth = .8) +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::scale_color_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::labs(x = if (variable == "length") "Parent length in episode 1" else "Parent ERRO in episode 1", y = "Probability density",
      caption = paste("Pale: all adults. Solid:", if (matched) "one parent representation per sibling match." else "one parent representation per sampled offspring.",
        if (matched) "Dashed: parents weighted by ERRO1 * ERRO2. No matches gives no positive solid density." else
          "Dashed: expected parental contributions, weighted by each episode's sample size. Traits refer to episode 1." , sep = "\n")) + .pop_plot_theme()
}
