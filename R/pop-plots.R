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
#' @param y_scale Either "linear" (default) or "log10". On a log10 axis,
#' positive bin counts are drawn from a baseline of 0.5 so a count of one is
#' visible. Empty bins are omitted, while the zero-offspring bin is retained.
#' @return A faceted ggplot, with females and males in consistent panel positions.
#' @export
pop_plot_reproduction <- function(replicate, sampled = FALSE, y_scale = c("linear", "log10")) {
  y_scale <- match.arg(y_scale)
  counts <- pop_parent_counts(replicate, sampled)
  # Large ranges use wider bins to keep extreme realizations legible.
  width <- max(1, ceiling((max(counts$offspring) + 1) / 60))
  if (y_scale == "log10") {
    bins <- .pop_histogram_bins(counts, width)
    geometry <- ggplot2::geom_rect(data = bins,
      ggplot2::aes(xmin = lower, xmax = upper, ymin = 0.5, ymax = n_parents, fill = sex),
      inherit.aes = FALSE, color = "white")
    axis <- ggplot2::scale_y_log10(limits = c(0.5, NA), breaks = .pop_log_count_breaks)
  } else {
    geometry <- ggplot2::geom_histogram(binwidth = width, boundary = -0.5, color = "white")
    axis <- ggplot2::scale_y_continuous(breaks = .pop_count_breaks)
  }
  caption <- "Parents with zero offspring are included.\nPanel heights may use different scales."
  if (y_scale == "log10") caption <- paste0(caption, "\nLog10 axis: bar baseline 0.5; empty bins omitted.")
  ggplot2::ggplot(counts, ggplot2::aes(x = offspring, fill = sex, weight = n_parents)) +
    geometry +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330"),
                                 limits = c("Female", "Male")) +
    ggplot2::scale_x_continuous(breaks = .pop_count_breaks, limits = c(-0.5, NA)) +
    axis +
    ggplot2::labs(x = if (sampled) "Sampled offspring per sampled parent" else
                    "Total offspring per parent",
                  y = if (y_scale == "log10") "Number of parents (log10 scale)" else "Number of parents",
                  caption = caption) +
    .pop_plot_theme()
}

#' Aggregate reproductive frequencies into the histogram's shared bin boundaries
#'
#' Only occupied bins are returned. Bin boundaries match geom_histogram with
#' boundary -0.5 and integer binwidth, without expanding frequencies into rows.
#' @param counts Frequency tibble with offspring, sex, and n_parents.
#' @param width Positive integer binwidth used in the linear histogram.
#' @return A data frame with sex, bin, n_parents, lower, and upper columns.
#' @keywords internal
.pop_histogram_bins <- function(counts, width) {
  counts$bin <- floor((counts$offspring + 0.5) / width)
  bins <- stats::aggregate(n_parents ~ sex + bin, data = counts, FUN = sum)
  bins$lower <- bins$bin * width - 0.5
  bins$upper <- bins$lower + width
  bins
}

#' Choose whole-count powers of ten for a log10 count axis
#' @param limits Axis limits supplied by ggplot2.
#' @return Tick positions 1, 10, 100, etc. within the visible upper limit.
#' @keywords internal
.pop_log_count_breaks <- function(limits) {
  ticks <- 10^seq(0, max(0, ceiling(log10(max(limits)))))
  ticks[ticks <= max(limits)]
}

#' Evaluate a gamma reproductive sampling-weight density for preview
#'
#' Scale is one. The plotted range ends at the 99.5th percentile. For shape < 1
#' the zero singularity is excluded by starting at 0.001 times the upper range;
#' a mixed geometric/linear grid resolves the steep left edge. Values are raw
#' gamma densities, not renormalized to the displayed range. No RNG is used.
#' @param shape Gamma shape between 0.05 and 1000.
#' @return A tibble with sampling_weight and density columns.
#' @export
pop_gamma_density <- function(shape) {
  .pop_shape(shape)
  upper <- stats::qgamma(0.995, shape = shape, scale = 1)
  lower <- upper / 1000
  x <- sort(unique(c(if (shape >= 1) 0,
                     exp(seq(log(lower), log(upper), length.out = 300)),
                     seq(lower, upper, length.out = 500))))
  tibble::tibble(sampling_weight = x, density = stats::dgamma(x, shape = shape, scale = 1))
}

#' Plot a theoretical gamma density in the scenario's sex-specific color
#' @param shape Gamma shape; scale is fixed at one.
#' @param sex Either "Female" or "Male".
#' @param y_scale Linear or log10 density axis. Exact zero densities are omitted
#' on log10 axes; no artificial offset is added to the theoretical density.
#' @return A ggplot that can be drawn before any trial exists.
#' @export
pop_plot_gamma <- function(shape, sex = c("Female", "Male"), y_scale = c("linear", "log10")) {
  sex <- match.arg(sex)
  y_scale <- match.arg(y_scale)
  density <- pop_gamma_density(shape)
  if (y_scale == "log10") density <- density[density$density > 0, ]
  axis <- if (y_scale == "log10") ggplot2::scale_y_log10() else ggplot2::scale_y_continuous()
  ggplot2::ggplot(density, ggplot2::aes(x = sampling_weight, y = density)) +
    ggplot2::geom_line(color = if (sex == "Female") "#237e83" else "#bd7330", linewidth = 1) +
    axis + ggplot2::labs(x = "Reproductive sampling weight",
      y = if (y_scale == "log10") "Density (log10 scale)" else "Probability density") +
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

#' Evaluate fixed-draw model previews without touching the trial RNG
#' @param size_model Named sex-specific model settings.
#' @param n Number of numerical integration draws per sex.
#' @return Adult traits approximating the model distribution, using fixed seeds.
#' @export
pop_size_preview <- function(size_model, n = 12000) {
  .pop_with_seed(927413, pop_generate_size_population(n, n, size_model))
}

#' Plot growth or length-based reproduction with sex-aligned panels
#'
#' Corresponding population/sample views use the same axes derived from compact
#' full-population bounds. Large samples are displayed as an 80 by 80 grid with
#' point area proportional to counts; every sampled adult contributes.
#' @param adults Adult records to display (display subsample, full survey, or preview).
#' @param size_model Saved sex-specific parameters.
#' @param kind Growth (age-length) or output (length-ERRO).
#' @param summary Optional full-population bins to set shared population/sample axes.
#' @param curves_only Draw reference curves without individual points.
#' @return A faceted ggplot, females left and males right.
#' @export
pop_plot_size_relationship <- function(adults, size_model, kind = c("growth", "output"),
                                       summary = NULL, curves_only = FALSE) {
  kind <- match.arg(kind)
  curves <- do.call(rbind, lapply(c("Female", "Male"), function(s) {
    p <- size_model[[s]]
    age <- seq(p$min_age, p$max_age, length.out = 250)
    length <- p$Linf * (-expm1(-p$K * (age - p$t0)))
    a <- adults[adults$sex == s, ]
    log_norm <- log(a$expected_reproductive_weight[1]) - p$exponent * log(a$length[1])
    if (kind == "growth") data.frame(sex = s, x = age, y = length)
    else {
      upper <- if (!is.null(summary)) max(summary$upper[summary$variable == "length"]) else max(c(length, a$length))
      length <- seq(0, upper, length.out = 250)
      data.frame(sex = s, x = length,
                 y = if (p$exponent == 0) rep(exp(log_norm), length(length)) else exp(p$exponent * log(length) + log_norm))
    }
  }))
  points <- data.frame(sex = adults$sex,
    x = if (kind == "growth") adults$age else adults$length,
    y = if (kind == "growth") adults$length else adults$reproductive_weight, n = 1)
  # Fixed bounds make population and survey views directly comparable.
  xmax <- if (kind == "growth") max(vapply(size_model, `[[`, numeric(1), "max_age")) else
    if (!is.null(summary)) max(summary$upper[summary$variable == "length"]) else max(points$x, curves$x)
  ymax <- if (!is.null(summary)) max(summary$upper[summary$variable ==
      if (kind == "growth") "length" else "reproductive_weight"]) else max(points$y, curves$y)
  binned <- nrow(points) > 10000 && !curves_only
  if (binned) {
    points$xb <- pmin(79, floor(80 * points$x / xmax))
    points$yb <- pmin(79, floor(80 * points$y / ymax))
    points <- stats::aggregate(n ~ sex + xb + yb, points, sum)
    points$x <- (points$xb + 0.5) * xmax / 80
    points$y <- (points$yb + 0.5) * ymax / 80
  }
  plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = curves, ggplot2::aes(x, y, color = sex), linewidth = 0.8)
  if (!curves_only) {
    if (binned) {
      plot <- plot + ggplot2::geom_point(data = points,
        ggplot2::aes(x, y, color = sex, size = n), alpha = 0.5) +
        ggplot2::scale_size_area(max_size = 5)
    } else {
      plot <- plot + ggplot2::geom_point(data = points,
        ggplot2::aes(x, y, color = sex), size = 0.8, alpha = 0.3)
    }
  }
  plot + ggplot2::facet_wrap(~sex, nrow = 1) +
    ggplot2::scale_color_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(x = if (kind == "growth") "Age" else "Length",
      y = if (kind == "growth") "Length" else "ERRO",
      caption = paste(if (kind == "growth") "Line: growth at mean Linf and K." else
        "Line: length^exponent,\non the same relative scale as ERRO.",
        if (binned) "\nAll adults contribute; point areas\nshow counts in grid cells." else "")) + .pop_plot_theme()
}

#' Plot full-population ERRO histogram densities from compact summaries
#' @param summary Full-population trait bins from pop_summarize_size().
#' @param y_scale Linear or log10 density axis.
#' @return Sex-faceted histogram of assigned ERRO, not realized offspring counts.
#' @export
pop_plot_erro <- function(summary, y_scale = c("linear", "log10")) {
  y_scale <- match.arg(y_scale)
  x <- summary[summary$variable == "reproductive_weight", ]
  x$density <- x$count / ave(x$count, x$sex, FUN = sum) / (x$upper - x$lower)
  if (y_scale == "log10") {
    x <- x[x$density > 0, ]
    baseline <- min(x$density) / 2
  } else baseline <- 0
  plot <- ggplot2::ggplot(x) +
    ggplot2::geom_rect(ggplot2::aes(xmin = lower, xmax = upper, ymin = baseline,
                                    ymax = density, fill = sex), color = "white") +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::labs(x = "Expected relative reproductive output (ERRO)", y = "Density",
      caption = if (y_scale == "log10") "Positive-density bins only; bar baseline is half the smallest density." else
        "ERRO has mean one within each sex.\nThese are reproductive sampling weights.") + .pop_plot_theme()
  if (y_scale == "log10") plot + ggplot2::scale_y_log10() else plot
}

#' Compare adults with parents represented once per identified pair
#'
#' All distributions use the saved full-population bin edges. Observed parents
#' are joined with match(), intentionally repeating adults with several matched
#' offspring. Expected pair-parent bin probabilities are sums of actual ERRO,
#' divided by total ERRO: the finite-population version of r*f(r). This is the
#' expectation averaging over reproduction and uniform surveys, rather than
#' a conditional prediction for a particular realized pedigree or adult sample.
#' @param replicate A saved size-based replicate.
#' @param variable ERRO (reproductive_weight) or length.
#' @return Histogram-bin tibble, separately normalized by sex and source. Empty
#' match distributions have zero density and are explicitly labeled by the plot.
#' @export
pop_pair_trait_data <- function(replicate, variable = c("reproductive_weight", "length")) {
  variable <- match.arg(variable)
  do.call(rbind, lapply(c("Female", "Male"), function(s) {
    bins <- replicate$size_summary[replicate$size_summary$sex == s &
                                    replicate$size_summary$variable == variable, ]
    adults <- if (s == "Female") replicate$samples$females else replicate$samples$males
    ids <- replicate$pairs$parent_id[replicate$pairs$sex == s]
    values <- list(`Population adults` = NULL, `Sampled adults` = adults[[variable]],
                   `Identified pairs` = adults[[variable]][match(ids, adults$parent_id)],
                   `Expected pair parents` = NULL)
    do.call(rbind, lapply(names(values), function(source) {
      counts <- if (source == "Population adults") bins$count else
        if (source == "Expected pair parents") bins$erro_sum else
          tabulate(findInterval(values[[source]], c(bins$lower, tail(bins$upper, 1)),
                                all.inside = TRUE), nbins = nrow(bins))
      data.frame(sex = s, source = source, lower = bins$lower, upper = bins$upper,
                 count = counts, total = sum(counts),
                 density = if (sum(counts) > 0) counts / sum(counts) / (bins$upper - bins$lower) else 0)
    }))
  }))
}

#' Plot normalized adult and identified-pair trait histograms with expectations
#' @param replicate A saved size-based replicate.
#' @param variable ERRO or length.
#' @return A printable pop_pair_traits object containing two ggplots. Each sex
#' has its own secondary axis: count = density * identified pairs * bin width.
#' With no matches that axis is omitted. Printing stacks panels on narrow devices.
#' @export
pop_plot_pair_traits <- function(replicate, variable = c("reproductive_weight", "length")) {
  variable <- match.arg(variable)
  x <- pop_pair_trait_data(replicate, variable)
  population <- x[x$source == "Population adults", ]
  outlines <- x[x$source %in% c("Sampled adults", "Identified pairs") & x$total > 0, ]
  # Explicit top and right edges produce true histogram outlines, including zeros.
  segments <- rbind(transform(outlines, x = lower, xend = upper, y = density, yend = density),
    do.call(rbind, lapply(split(outlines, interaction(outlines$sex, outlines$source), drop = TRUE), function(d) {
      transform(d, x = lower, xend = lower, y = c(0, head(density, -1)), yend = density)
    })))
  expected <- x[x$source == "Expected pair parents", ]
  labels <- data.frame(sex = c("Female", "Male"))
  labels$label <- vapply(labels$sex, function(s) {
    n <- sum(replicate$pairs$sex == s)
    if (n) paste0("Pairs: ", n) else "No identified pairs"
  }, character(1))
  panels <- lapply(c("Female", "Male"), function(s) {
    subset_sex <- function(d) d[d$sex == s, ]
    n <- sum(replicate$pairs$sex == s)
    bin_width <- population$upper[population$sex == s][1] - population$lower[population$sex == s][1]
    # Each sex needs its own transformation; a common faceted secondary axis
    # would be wrong whenever match counts or bin widths differ between sexes.
    secondary <- if (n > 0) ggplot2::sec_axis(~ . * n * bin_width,
      name = "Identified pairs per bin", breaks = .pop_count_breaks) else ggplot2::waiver()
    ggplot2::ggplot() +
    ggplot2::geom_rect(data = subset_sex(population), ggplot2::aes(xmin = lower, xmax = upper,
      ymin = 0, ymax = density, fill = sex), alpha = 0.16) +
    ggplot2::geom_segment(data = subset_sex(segments), ggplot2::aes(x, y, xend = xend, yend = yend,
      color = sex, linewidth = source, linetype = source)) +
    ggplot2::geom_line(data = subset_sex(expected), ggplot2::aes(x = (lower + upper) / 2, y = density,
      group = sex), linetype = "dashed", color = "#233a3d", linewidth = 0.8) +
    ggplot2::geom_text(data = subset_sex(labels), ggplot2::aes(x = Inf, y = Inf, label = label),
      hjust = 1.05, vjust = 1.5, size = 3.4) +
    ggplot2::facet_wrap(~sex, nrow = 1, scales = "free") +
    ggplot2::scale_fill_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::scale_color_manual(values = c(Female = "#237e83", Male = "#bd7330")) +
    ggplot2::scale_linewidth_manual(values = c(`Sampled adults` = 0.45, `Identified pairs` = 1)) +
    ggplot2::scale_linetype_manual(values = c(`Sampled adults` = "dotted", `Identified pairs` = "solid")) +
    ggplot2::scale_y_continuous(sec.axis = secondary) +
    ggplot2::labs(x = if (variable == "length") "Length" else "ERRO", y = "Probability density") +
    .pop_plot_theme() + ggplot2::theme(axis.title.y.right = ggplot2::element_text(
      color = if (s == "Female") "#237e83" else "#bd7330"))
  })
  structure(list(panels = stats::setNames(panels, c("Female", "Male")),
    caption = paste("Pale bars: population adults. Dotted: sampled adults.",
      "Solid: identified pairs (repeats). Dark dashed: expected pair parents.",
      "Right axis: counts for the solid outline only.", sep = "\n")), class = "pop_pair_traits")
}

#' Draw independently scaled pair-trait panels without an extra layout package
#' @param x A composition returned by pop_plot_pair_traits().
#' @param ... Reserved for compatibility with print().
#' @return Invisibly returns x after drawing on the current graphics device.
#' @keywords internal
print.pop_pair_traits <- function(x, ...) {
  wide <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) >= 7.8
  caption <- if (wide) x$caption else paste(
    "Pale bars: population adults.", "Dotted: sampled adults.",
    "Solid: identified pairs (repeats).", "Dark dashed: expected pair parents.",
    "Right axis: solid-outline counts only.", sep = "\n")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = if (wide) 2 else 3, ncol = if (wide) 2 else 1,
    heights = grid::unit.c(if (wide) grid::unit(1, "null") else grid::unit(c(1, 1), "null"),
                          grid::unit(if (wide) 50 else 75, "points")))))
  for (i in seq_along(x$panels)) {
    grid::pushViewport(grid::viewport(layout.pos.row = if (wide) 1 else i,
                                      layout.pos.col = if (wide) i else 1))
    grid::grid.draw(ggplot2::ggplotGrob(x$panels[[i]]))
    grid::popViewport()
  }
  grid::pushViewport(grid::viewport(layout.pos.row = if (wide) 2 else 3,
                                    layout.pos.col = if (wide) 1:2 else 1))
  grid::grid.text(caption, x = grid::unit(8, "points"), y = grid::unit(1, "npc"),
                  just = c("left", "top"), gp = grid::gpar(fontsize = 9))
  grid::popViewport(2)
  invisible(x)
}
