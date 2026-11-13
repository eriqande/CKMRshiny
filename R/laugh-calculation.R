# The editor retains year-by-age proportions. microscoping >= 1.1 accepts
# year-by-age NUMBERS through samp_nya. Keep this interface conversion explicit.

#' Convert an editor design to the installed laugh-test matrix interface
#' @param design Normalized editor design from laugh_prepare_inputs().
#' @return Named arguments for ckmr_laugh_test(), with samp_nya containing
#' expected sample numbers (row proportions times annual effort). Display names
#' are removed because offarray cannot combine named dimensions with age-zero
#' offsets. Row and column order are preserved. No package loading is required.
#' @export
laugh_function_arguments <- function(design) {
  p <- design$samp_pata; ny <- design$samp_ny
  if (!is.matrix(p) || nrow(p)!=length(ny) || ncol(p)!=length(design$nata) ||
      any(!is.finite(p)) || any(p<0) || any(abs(rowSums(p)-1)>1e-10))
    stop("The sampling proportions must be a row-normalized year-by-age matrix.")
  counts <- unname(p * ny)
  if(any(!is.finite(counts)) || any(counts<0) || sum(counts)<=0)
    stop("Allocated sample numbers must be finite, nonnegative, and have a positive total.")
  if(!is.finite(sum(design$nata)) || !is.finite(sum(design$nata*design$wata*design$pmatata)))
    stop("The population or reproductive-output total is too large to represent numerically.")
  keys <- c("nata","wata","pmatata","MAX_JUVE_AGE","MIN_CATCH_CURVE_AGE",
            "beqy_POP","etriv","extraplus")
  args <- design[keys]
  args$samp_nya <- counts
  args$want_gory_details <- FALSE
  args
}

#' Check availability of the package's year-by-age sample-count interface
#' @return Installed version as a string, or an informative error. The package
#' is never installed automatically, and no local source checkout is modified.
#' @export
laugh_engine_version <- function() {
  if(!requireNamespace("microscoping",quietly=TRUE))
    stop("Install microscoping with samp_nya support to calculate kin counts. The input editor and exports remain available.")
  if(!"samp_nya" %in% names(formals(microscoping::ckmr_laugh_test)))
    stop("This microscoping installation lacks samp_nya matrix support. Reinstall Mark's updated package and restart Shiny.")
  as.character(utils::packageVersion("microscoping"))
}

#' Check the expanded age grid before allocating the package's large arrays
#'
#' Calls the installed package's own smooth_n helper, then follows its default
#' plus-group expansion to estimate the size of its age-by-year-by-parent-age
#' arrays. This preflight does not replace its kin calculations. The package
#' repeats its own preparation when called. Keep this adapter checked against
#' package updates; its equations match microscoping 1.1.2 with extraplus FALSE.
#' @param args Result of laugh_function_arguments().
#' @param max_work_cells Maximum cells in one principal five-dimensional array;
#' default five million. This bounds interactive workloads, not population size.
#' @return List with expanded ages/population, survival at the plus group,
#' aggregate supplied plus-group abundance, and estimated work-cell count.
#' @export
laugh_preflight <- function(args,max_work_cells=5e6) {
  laugh_engine_version()
  if(isTRUE(args$extraplus)) stop("This app uses the supplied plus group (extraplus = FALSE).")
  n <- args$nata; n[n==0] <- .001
  a <- length(n)-1L
  indexed <- offarray::offarray(unname(n),first=c(AGE=0),dim=length(n))
  smooth <- utils::getFromNamespace("smooth_n","microscoping")(
    indexed[args$MIN_CATCH_CURVE_AGE:(a-1L)])
  plus <- n[a+1L]
  ps <- plus / (smooth[a]+plus)
  if(!is.finite(ps) || ps<=0 || ps>=1 || any(!is.finite(smooth)))
    stop("The fitted population does not give a finite plus-group survival below one. Check the population curve and catch-curve starting age.")
  fec_pre <- sum(args$wata[seq_len(a)]*args$pmatata[seq_len(a)]*smooth)
  fec_plus <- args$wata[a+1L]*args$pmatata[a+1L]*plus
  if(!is.finite(fec_pre+fec_plus) || fec_pre+fec_plus<=0)
    stop("The smoothed population has no finite positive reproductive output.")
  ancient <- floor(1+max(log(args$etriv*(1+fec_pre/fec_plus))/log(ps),0))
  amax <- a+ancient
  cells <- (amax+1)^3 * nrow(args$samp_nya)^2
  if(!is.finite(cells) || cells>max_work_cells)
    stop("This design expands to ",format(amax+1,scientific=FALSE,trim=TRUE),
      " age classes across ",nrow(args$samp_nya)," years, exceeding the interactive calculation budget. ",
      "Reduce sampling years or review the plus-group/population schedule. Inputs can still be exported for a separate R calculation.")
  list(ages=0:amax,expanded_nata=c(smooth,plus*(1-ps)*ps^(0:ancient)),
       input_plus_age=a,input_plus_n=plus,plus_survival=ps,work_cells=cells)
}

#' Run Mark's installed laugh test on an immutable editor snapshot
#' @param design Normalized editor design from laugh_prepare_inputs().
#' @param max_work_cells Interactive array-size budget passed to laugh_preflight().
#' @return Plain list with saved design, exact function arguments, named numeric
#' results, package version, calculation time, warnings, and preparation metadata.
#' Large package environments/call attributes are discarded. This function calls
#' the installed package, never a copied or locally patched implementation.
#' @export
laugh_run_design <- function(design,max_work_cells=5e6) {
  version <- laugh_engine_version()
  args <- laugh_function_arguments(design)
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  computed <- withCallingHandlers({
    prep <- laugh_preflight(args,max_work_cells)
    raw <- do.call(microscoping::ckmr_laugh_test,args)
    expected <- stats::setNames(as.numeric(raw),names(raw))
    required <- c("POP_offJ","POP_offA","HSP_JJ","GGP_JJ","nsamp_J","nsamp_A")
    if(!all(required %in% names(expected)) || any(!is.finite(expected[required])) || any(expected[required]<0))
      stop("microscoping returned missing or nonfinite expected counts; this calculation was not saved.")
    list(expected=expected[required],preparation=prep)
  },warning=function(w) {warnings <<- unique(c(warnings,conditionMessage(w))); invokeRestart("muffleWarning")})
  c(list(design=design,function_arguments=args,package_version=version,created_at=Sys.time(),
         elapsed=proc.time()[["elapsed"]]-start,warnings=warnings),computed)
}

#' Label and explain the package's six expected-count outputs
#' @param calculation Result of laugh_run_design().
#' @return Data frame of raw metric names, descriptive labels, numeric values,
#' and interpretations. Counts include both sexes and all sampling years.
#' @export
laugh_result_table <- function(calculation) {
  keys <- c("POP_offJ","POP_offA","HSP_JJ","GGP_JJ","nsamp_J","nsamp_A")
  data.frame(metric=keys,label=c("POPs: juvenile offspring","POPs: adult offspring",
    "HSPs: juvenile-juvenile","GGPs: juvenile-juvenile","Juvenile samples","Adult samples"),
    expected=unname(calculation$expected[keys]),meaning=c(
      "Parent-offspring matches with offspring at or below the juvenile cutoff when sampled.",
      "Parent-offspring matches with offspring above that cutoff when sampled.",
      "Juvenile half-sibling matches; same-birth-cohort comparisons are excluded.",
      "Juvenile grandparent-grandchild matches that can resemble half-sibling matches genetically.",
      "Allocated samples at ages up to and including the juvenile cutoff.",
      "Allocated samples above the juvenile cutoff."),stringsAsFactors=FALSE)
}

#' Plot the four types of expected kin counts from one saved calculation
#' @param calculation Result of laugh_run_design().
#' @return Horizontal bar plot; fractional counts are retained and labelled.
#' @export
laugh_plot_results <- function(calculation) {
  x <- laugh_result_table(calculation)[1:4,]
  x$label <- factor(x$label,levels=rev(x$label))
  x$text <- vapply(x$expected,function(v) format(signif(v,4),scientific=FALSE,trim=TRUE),character(1))
  ggplot2::ggplot(x,ggplot2::aes(expected,label,fill=metric)) +
    ggplot2::geom_col(width=.65) +
    ggplot2::geom_text(ggplot2::aes(label=text),hjust=-.15,size=4) +
    ggplot2::scale_fill_manual(values=c(POP_offJ="#237e83",POP_offA="#8065a6",HSP_JJ="#bd7330",GGP_JJ="#697d89")) +
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,.2))) +
    ggplot2::labs(x="Expected kin-pair count",y=NULL,
      caption="Expected counts can be fractional. These are rough scoping calculations, not observed pairs or fitted abundance estimates.") +
    .pop_plot_theme()
}
