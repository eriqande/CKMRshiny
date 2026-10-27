#' Get reproductive frequency summaries, including zero-offspring parents
#' @param replicate A saved replicate from pop_simulate_rep().
#' @param sampled If TRUE, restrict to sampled parents and sampled offspring.
#' @return A tibble with sex, offspring, and n_parents columns. No unsampled
#' individual records are reconstructed for plotting.
#' @export
pop_parent_counts <- function(replicate, sampled = FALSE) {
  if (!sampled) return(replicate$reproduction_summary)
  pop_summarize_reproduction(rbind(replicate$samples$females, replicate$samples$males),
                              replicate$samples$kids)
}

#' Common visual theme for the parent-offspring teaching app
#' @return A ggplot2 theme.
#' @keywords internal
.pop_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(hjust = 0),
                   legend.position = "none")
}

#' Choose nonnegative whole-number ticks for count axes
#' @param limits Numeric axis limits supplied by ggplot2.
#' @return Unique nonnegative integer tick positions.
#' @keywords internal
.pop_count_breaks <- function(limits) {
  unique(pmax(0, round(pretty(limits, n = 6))))
}

#' Plot realized offspring counts in a population or sample
#' @param replicate Saved replicate.
#' @param sampled Use sampled parents and sampled offspring only.
#' @return A faceted ggplot, with females and males in consistent panel positions.
#' @export
pop_plot_reproduction <- function(replicate, sampled = FALSE) {
  counts <- pop_parent_counts(replicate, sampled)
  # Large ranges use wider bins to keep extreme realizations legible.
  width <- max(1, ceiling((max(counts$offspring) + 1) / 60))
  ggplot2::ggplot(counts, ggplot2::aes(x = offspring, fill = sex, weight = n_parents)) +
    ggplot2::geom_histogram(binwidth = width, boundary = -0.5, color = "white") +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330"),
                                 limits = c("Female", "Male")) +
    ggplot2::scale_x_continuous(breaks = .pop_count_breaks, limits = c(-0.5, NA)) +
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks) +
    ggplot2::labs(x = if (sampled) "Sampled offspring per sampled parent" else
                    "Total offspring per parent", y = "Number of parents",
                  caption = "Parents with zero offspring are included.\nPanel heights may use different scales.") +
    .pop_plot_theme()
}

#' Plot the distribution of pair counts across replicates
#' @param trial Saved trial.
#' @param rep_id Replicate to highlight with a dashed line.
#' @return A ggplot with theoretical means (red) and selected replicate (dashed).
#' @export
pop_plot_pairs <- function(trial, rep_id = 1) {
  x <- trial$summary[trial$summary$method == "Moments", ]
  settings <- trial$settings
  expected <- data.frame(sex = c("Female", "Male"),
                          mean = c(settings$Sf / settings$Nf,
                                   settings$Sm / settings$Nm) * settings$Sk)
  chosen <- x[x$rep_id == rep_id, ]
  width <- max(1, ceiling((diff(range(x$matches)) + 1) / 40))
  ggplot2::ggplot(x, ggplot2::aes(x = matches, fill = sex)) +
    ggplot2::geom_histogram(binwidth = width, boundary = -0.5, color = "white") +
    ggplot2::geom_vline(data = expected, ggplot2::aes(xintercept = mean),
                         color = "#bc423c", linewidth = 1) +
    ggplot2::geom_vline(data = chosen, ggplot2::aes(xintercept = matches),
                         linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~sex, nrow = 1) +
    ggplot2::scale_x_continuous(breaks = .pop_count_breaks) +
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks) +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330"),
                                 limits = c("Female", "Male")) +
    ggplot2::labs(x = "Parent-offspring pairs in the sample", y = "Number of replicates",
                  caption = "Red: expected count = sampled parents * sampled kids / abundance.\nDashed: selected replicate.") +
    .pop_plot_theme()
}

#' Plot abundance estimates and explicitly count unbounded estimates
#'
#' Histograms contain finite estimates only. Every facet is annotated with its
#' unbounded count out of all replicates; empty finite distributions remain visible.
#' @param trial Saved trial.
#' @param rep_id Replicate highlighted with a dashed line if its estimate is finite.
#' @return A ggplot with methods in rows and sexes in columns.
#' @export
pop_plot_estimates <- function(trial, rep_id = 1) {
  x <- trial$summary
  x$method <- factor(x$method, levels = c("Moments", "Poisson pseudo-MLE"))
  labels <- unique(x[c("sex", "method", "truth")])
  labels$label <- vapply(seq_len(nrow(labels)), function(i) {
    rows <- x$sex == labels$sex[i] & x$method == labels$method[i]
    paste0(sum(!x$finite[rows]), " / ", sum(rows), " unbounded")
  }, character(1))
  finite <- x[x$finite, ]
  chosen <- finite[finite$rep_id == rep_id, ]
  ggplot2::ggplot(finite, ggplot2::aes(x = estimate, fill = sex)) +
    ggplot2::geom_histogram(bins = 35, color = "white") +
    ggplot2::geom_vline(data = labels, ggplot2::aes(xintercept = truth),
                         color = "#bc423c", linewidth = 1) +
    ggplot2::geom_vline(data = chosen, ggplot2::aes(xintercept = estimate),
                         linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = Inf, y = Inf, label = label),
                        inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3.7) +
    ggplot2::facet_grid(method ~ sex, scales = "free_x") +
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks) +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330"),
                                 limits = c("Female", "Male")) +
    ggplot2::labs(x = "Estimated adult abundance (finite estimates)", y = "Number of replicates",
                  caption = "Red: true abundance. Dashed: selected replicate, if finite.\nZero pairs give an unbounded estimate.") +
    .pop_plot_theme()
}

#' Prepare sorted replicate intervals for a single pair of sex panels
#'
#' The current moment and Poisson pseudo-MLE procedures have identical Wald
#' intervals; use the Moments rows once to avoid duplicating every replicate.
#' Sort within sex by lower endpoint, upper endpoint, then replicate ID, so tied
#' intervals have stable positions. Unavailable intervals retain NA rank.
#' @param trial Saved trial containing a replicate summary.
#' @return A tibble of interval rows with interval_rank increasing from bottom
#' to top within each sex. Zero-match replicates are retained without a rank.
#' @export
pop_interval_data <- function(trial) {
  x <- pop_wald_intervals(trial$summary[trial$summary$method == "Moments", ])
  x <- x[order(x$sex, x$lower, x$upper, x$rep_id, na.last = TRUE), ]
  x$interval_rank <- NA_integer_
  for (sex in c("Female", "Male")) {
    rows <- which(x$sex == sex & x$interval_available)
    x$interval_rank[rows] <- seq_along(rows)
  }
  x
}

#' Plot sorted horizontal Wald intervals for female and male abundance
#'
#' Thin intervals show all finite replicates. A thicker interval and point show
#' the selected replicate. The red line marks the true abundance. Zero-match
#' replicates have no interval and are counted explicitly, including when the
#' selected replicate has no interval. Both panels share abundance and rank axes.
#' @param trial Saved trial.
#' @param rep_id Selected replicate ID, independent of its sorted position.
#' @return A faceted ggplot. Negative lower limits are shown without clipping.
#' @export
pop_plot_intervals <- function(trial, rep_id = 1) {
  x <- pop_interval_data(trial)
  finite <- x[x$interval_available, ]
  chosen <- finite[finite$rep_id == rep_id, ]
  labels <- data.frame(sex = c("Female", "Male"),
                        truth = c(trial$settings$Nf, trial$settings$Nm))
  labels$label <- vapply(labels$sex, function(sex) {
    rows <- x$sex == sex
    paste0(sum(!x$interval_available[rows]), " / ", sum(rows), "\nno interval")
  }, character(1))
  missing <- x$sex[x$rep_id == rep_id & !x$interval_available]
  caption <- "Red: true abundance.\nThick line and point: selected replicate.\nNegative lower limits are retained."
  if (length(missing)) caption <- paste0(caption, "\nSelected replicate has no Wald interval:\n",
                                         paste(missing, collapse = " and "), ".")
  ggplot2::ggplot(finite, ggplot2::aes(color = sex)) +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper,
                                        y = interval_rank, yend = interval_rank),
                            linewidth = 0.3, alpha = 0.65) +
    ggplot2::geom_vline(data = labels, ggplot2::aes(xintercept = truth),
                          color = "#bc423c", linewidth = 0.8) +
    ggplot2::geom_segment(data = chosen, ggplot2::aes(x = lower, xend = upper,
                                                     y = interval_rank, yend = interval_rank),
                            linewidth = 1.5) +
    ggplot2::geom_point(data = chosen, ggplot2::aes(x = estimate, y = interval_rank),
                          size = 2) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = Inf, y = Inf, label = label),
                         inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3.2) +
    ggplot2::facet_wrap(~sex, nrow = 1) +
    ggplot2::scale_color_manual(values = c(Female = "#237e83", Male = "#bd7330"),
                                  limits = c("Female", "Male")) +
    ggplot2::scale_y_continuous(breaks = .pop_count_breaks,
                                  limits = c(0, max(c(1, finite$interval_rank)) * 1.15 + 1),
                                  expand = ggplot2::expansion(mult = c(0.02, 0))) +
    ggplot2::labs(x = "95% CI for adult abundance",
                    y = "Replicates, sorted by interval", caption = caption) +
    .pop_plot_theme()
}
