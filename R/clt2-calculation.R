#' Assemble sex-specific editors into clt2's sex/year/age count interface
#' @param state List containing F/M editor states, a two-row counts matrix,
#' and shared settings with the microscoping argument names.
#' @return Plain named argument list for do.call(microscoping::clt2, ...).
#' @export
clt2_function_arguments <- function(state) {
  cfg <- state$settings
  if(length(cfg$etriv)!=1 || !is.finite(cfg$etriv) || cfg$etriv<=0 || cfg$etriv>=1) stop("etriv must lie strictly between 0 and 1.")
  if(!is.logical(cfg$extraplus) || length(cfg$extraplus)!=1 || is.na(cfg$extraplus)) stop("Choose a plus-group treatment.")
  if(sum(state$counts)<=0) stop("At least one sex/year needs a positive sample count.")
  ny <- ncol(state$counts)
  designs <- lapply(c("F","M"),function(sex) laugh_prepare_inputs(state[[sex]],
    state$counts[sex,],cfg$MAX_JUVE_AGE,cfg$MIN_CATCH_CURVE_AGE,cfg$beqy_POP,allow_zero=TRUE))
  args <- list(n_sa=rbind(F=designs[[1]]$nata,M=designs[[2]]$nata),
    w_sa=rbind(F=designs[[1]]$wata,M=designs[[2]]$wata),
    pmat_sa=rbind(F=designs[[1]]$pmatata,M=designs[[2]]$pmatata))
  counts <- array(0,c(2,ny,ncol(args$n_sa)),dimnames=list(SEX=c("F","M"),YEAR=NULL,AGE=NULL))
  for(i in 1:2) counts[i,,] <- designs[[i]]$samp_pata * state$counts[i,]
  args <- c(args,list(nsamp_sya=counts),cfg,list(want_gory_details=FALSE))
  if(.clt2_parent_sex_supported()) args$by_parent_sex <- TRUE
  args
}

#' Preflight sex-specific smoothing and the common expanded age grid
#' @param args Native clt2 argument list.
#' @param max_work_cells Budget including two focal sexes and parent sex.
#' @return Metadata for both populations and the estimated principal work size.
#' Uses microscoping's demographic helper, never its kin probabilities. The
#' expansion follows clt2 1.1.4, including its common maximum age for both sexes.
#' @export
clt2_preflight <- function(args,max_work_cells=5e6) {
  if(!requireNamespace("microscoping",quietly=TRUE) || !"nsamp_sya" %in% names(formals(microscoping::clt2))) stop("Install microscoping with clt2 sex/year/age sampling-array support.")
  a <- ncol(args$n_sa)-1L
  prep <- lapply(1:2,function(i) {
    n <- args$n_sa[i,]; n[n==0] <- .001
    indexed <- offarray::offarray(unname(n),first=c(AGE=0),dim=length(n))
    sm <- utils::getFromNamespace("smooth_n","microscoping")(indexed[args$MIN_CATCH_CURVE_AGE:(a-1L)])
    ps <- if(isTRUE(args$extraplus)) sm[a]/sm[a-1] else n[a+1]/(sm[a]+n[a+1])
    plus <- if(isTRUE(args$extraplus)) sm[a]*ps/(1-ps) else n[a+1]
    if(!is.finite(ps) || ps<=0 || ps>=1 || any(!is.finite(sm))) stop("Check the population schedules: fitted survival must be finite and below one for each sex.")
    pre <- sum(args$w_sa[i,seq_len(a)]*args$pmat_sa[i,seq_len(a)]*sm)
    fp <- args$w_sa[i,a+1]*args$pmat_sa[i,a+1]*plus
    ancient <- floor(1+max(log(args$etriv*(1+pre/fp))/log(ps),0))
    list(smooth=sm,plus=plus,input_plus_n=args$n_sa[i,a+1],input_plus_age=a,plus_survival=ps,ancient=ancient)
  })
  extra <- max(vapply(prep,`[[`,numeric(1),"ancient")); amax <- a+extra
  cells <- 8*(amax+1)^3 * dim(args$nsamp_sya)[2]^2
  if(!is.finite(cells) || cells>max_work_cells) stop("This two-sex design exceeds the interactive calculation budget (",amax+1," expanded ages). Reduce sampling years or export the inputs for a separate R calculation.")
  for(i in 1:2) prep[[i]]$expanded_nata <- c(prep[[i]]$smooth,prep[[i]]$plus*(1-prep[[i]]$plus_survival)*prep[[i]]$plus_survival^(0:extra))
  names(prep) <- c("F","M")
  list(sexes=prep,ages=0:amax,work_cells=cells)
}

#' Run the installed two-sex laugh test and retain a small immutable result
#' @param args Native clt2 argument list.
#' @param max_work_cells Interactive workload limit, passed to clt2_preflight().
#' @return Saved inputs, plain expected counts, warnings and demographic metadata.
#' @export
clt2_run_design <- function(args,max_work_cells=5e6) {
  if(.clt2_parent_sex_supported() && is.null(args$by_parent_sex)) args$by_parent_sex <- TRUE
  warnings <- character(); start <- proc.time()[["elapsed"]]
  r <- withCallingHandlers({
    p <- clt2_preflight(args,max_work_cells)
    result <- do.call(microscoping::clt2,args)
    expected <- stats::setNames(as.numeric(result),names(result))
    keys <- c("POP_offJ","POP_offA","HSP_JJ","GGP_JJ","nsamp_J","nsamp_A")
    if(!all(keys %in% names(expected)) || any(!is.finite(expected[keys])) || any(expected[keys]<0)) stop("microscoping returned invalid counts; the calculation was not saved.")
    by_sex <- attr(result,"kin_by_parent_sex")
    if(!is.null(by_sex)) {
      sex_keys <- c("POP_offJ","POP_offA","HSP_JJ")
      if(!is.matrix(by_sex) || !identical(dim(by_sex),c(2L,3L)) ||
         !identical(rownames(by_sex),c("F","M")) || !identical(colnames(by_sex),sex_keys) ||
         any(!is.finite(by_sex)) || any(by_sex<0) ||
         any(abs(colSums(by_sex)-expected[sex_keys])>1e-9*pmax(1,expected[sex_keys])))
        stop("Parent-sex counts are invalid or do not sum to the combined totals; the calculation was not saved.")
    } else if(isTRUE(args$by_parent_sex)) stop("Parent-sex reporting was requested but microscoping returned no breakdown.")
    list(expected=expected[keys],preparation=p,parent_sex_expected=by_sex)
  },warning=function(w) { warnings <<- unique(c(warnings,conditionMessage(w))); invokeRestart("muffleWarning") })
  c(list(function_arguments=args,package_version=as.character(utils::packageVersion("microscoping")),
    created_at=Sys.time(),elapsed=proc.time()[["elapsed"]]-start,warnings=warnings),r)
}

#' Check whether the installed clt2 offers optional parent-sex accounting
#' @return Logical; older installations still support the combined results.
#' @keywords internal
.clt2_parent_sex_supported <- function() requireNamespace("microscoping",quietly=TRUE) &&
  "by_parent_sex" %in% names(formals(microscoping::clt2))

#' Tidy the optional POP/HSP counts by parent sex from a saved calculation
#' @param calculation Result of clt2_run_design().
#' @return Data frame with metric, label, parent_sex, and expected count, or NULL
#' for calculations saved without the new package attribute. F/M refer to the
#' sampled parent for POPs and unobserved shared parent for HSPs.
#' @export
clt2_parent_sex_table <- function(calculation) {
  m <- calculation$parent_sex_expected
  if(is.null(m)) return(NULL)
  keys <- colnames(m)
  labels <- laugh_result_table(calculation)
  data.frame(metric=rep(keys,each=2),label=rep(labels$label[match(keys,labels$metric)],each=2),
    parent_sex=rep(rownames(m),length(keys)),expected=as.numeric(m),stringsAsFactors=FALSE)
}

#' Plot the optional parent-sex breakdown below the combined results
#' @param calculation Result of clt2_run_design() with parent_sex_expected.
#' @return Grouped horizontal bar plot for female and male parent contributions.
#' @export
clt2_plot_parent_sex <- function(calculation) {
  x <- clt2_parent_sex_table(calculation)
  if(is.null(x)) stop("This saved calculation has no parent-sex breakdown.")
  x$label <- factor(x$label,levels=rev(unique(x$label)))
  x$text <- vapply(x$expected,function(v)format(signif(v,4),scientific=FALSE,trim=TRUE),character(1))
  dodge <- ggplot2::position_dodge(width=.75,orientation="y")
  ggplot2::ggplot(x,ggplot2::aes(expected,label,fill=parent_sex,group=parent_sex))+
    ggplot2::geom_col(position=dodge,width=.65)+
    ggplot2::geom_text(ggplot2::aes(label=text),position=dodge,hjust=-.12,size=3.5)+
    ggplot2::scale_fill_manual(values=c(F="#237e83",M="#bd7330"),breaks=c("F","M"),labels=c("Female","Male"),name="Parent / shared parent")+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,.22)))+
    ggplot2::labs(x="Expected kin-pair count",y=NULL,
      caption="POPs: sampled parent's sex. HSPs: unobserved shared parent's sex.
Female and male contributions sum to the combined totals above; GGPs remain combined.")+.pop_plot_theme()
}
