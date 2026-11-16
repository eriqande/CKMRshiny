# Argument-list interchange for the two microscoping editors. No R expressions
# are evaluated: imports are data, and validation precedes all editor mutations.

#' Read a microscoping argument list into an editable design
#' @param inputs Named list suitable for do.call(ckmr_laugh_test, ...) or
#' do.call(clt2, ...). Full sample counts take precedence over sampling margins.
#' @return List with kind (laugh/clt2), editor states, annual counts, settings,
#' and notices. Inputs are copied; no package calculation or smoothing occurs.
#' Editors support ages 0 through a 2-200 plus group and 1-12 sampling years.
#' @export
ckmr_import_design <- function(inputs) {
  if (!is.list(inputs) || is.null(names(inputs)) || anyNA(names(inputs)) || any(!nzchar(names(inputs))) || anyDuplicated(names(inputs)))
    stop("Supply a named argument list with unique names.")
  if (is.list(inputs$function_arguments)) inputs <- inputs$function_arguments
  kind <- if ("n_sa" %in% names(inputs)) "clt2" else "laugh"
  pop_keys <- if(kind=="clt2") c("n_sa","w_sa","pmat_sa") else c("nata","wata","pmatata")
  if(!all(c(pop_keys,"MAX_JUVE_AGE") %in% names(inputs))) stop("Missing inputs: ",paste(setdiff(c(pop_keys,"MAX_JUVE_AGE"),names(inputs)),collapse=", "))
  numeric_values <- function(x,label) {
    if(!is.numeric(x) || any(!is.finite(x)) || any(x<0)) stop(label," must contain finite, nonnegative numbers.")
  }
  pop <- inputs[pop_keys]
  for(k in pop_keys) numeric_values(pop[[k]],k)
  sex_order <- if(kind=="clt2") rownames(pop[[1]]) else NULL
  if(is.null(sex_order)) sex_order <- c("F","M")
  sex_matrix <- function(x,label,order=sex_order) {
    if(length(dim(x))!=2 || nrow(x)!=2) stop(label," must have two sex rows.")
    labs <- rownames(x); if(is.null(labs)) labs <- order
    if(!identical(sort(labs),c("F","M"))) stop(label," sex labels must be F and M.")
    if(inherits(x,"offarray") && !all(offarray::dimseq(x)[[2]]==0:(ncol(x)-1))) stop(label," ages must start at 0.")
    array(as.numeric(x),dim=unname(dim(x)))[match(c("F","M"),labs),,drop=FALSE]
  }
  if(kind=="clt2") pop <- Map(sex_matrix,pop,pop_keys) else {
    for(k in pop_keys) {
      if(length(dim(pop[[k]]))>1) stop(k," must be a vector.")
      if(inherits(pop[[k]],"offarray") && !all(offarray::dimseq(pop[[k]])[[1]]==0:(length(pop[[k]])-1))) stop(k," ages must start at 0.")
    }
    pop <- lapply(pop,as.numeric)
  }
  sizes <- if(kind=="clt2") vapply(pop,ncol,integer(1)) else lengths(pop)
  if(length(unique(sizes))!=1 || sizes[1]<3 || sizes[1]>201) stop("Population schedules must share an age grid with a plus group from 2 to 200.")
  a <- sizes[1]-1L
  settings <- list(MAX_JUVE_AGE=2,MIN_CATCH_CURVE_AGE=0,beqy_POP=FALSE,etriv=.01,extraplus=FALSE)
  for(k in intersect(names(settings),names(inputs))) settings[[k]] <- inputs[[k]]
  for(k in c("MAX_JUVE_AGE","MIN_CATCH_CURVE_AGE")) {
    v <- settings[[k]]
    if(!is.numeric(v) || length(v)!=1 || !is.finite(v) || v<0 || v!=floor(v)) stop(k," must be a nonnegative whole number.")
    settings[[k]] <- as.integer(v)
  }
  if(settings$MAX_JUVE_AGE>=a || settings$MIN_CATCH_CURVE_AGE>a-2) stop("The juvenile cutoff must be below the plus group, and the catch-curve start must leave at least two exact ages.")
  for(k in c("beqy_POP","extraplus")) if(!is.logical(settings[[k]]) || length(settings[[k]])!=1 || is.na(settings[[k]])) stop(k," must be TRUE or FALSE.")
  if(!is.numeric(settings$etriv) || length(settings$etriv)!=1 || !is.finite(settings$etriv) || settings$etriv<=0 || settings$etriv>=1) stop("etriv must lie strictly between 0 and 1.")
  if(kind=="laugh" && settings$extraplus) stop("The one-sex editor currently supports only extraplus = FALSE; the clt2 editor also supports TRUE.")
  notices <- character()
  known <- c(pop_keys,names(settings),"want_gory_details","by_parent_sex",if(kind=="clt2") c("nsamp_sya","samp_p_sa","nsamp_y") else c("samp_nya","samp_pata","samp_ny"))
  unknown <- setdiff(names(inputs),known)
  if(length(unknown)) notices <- c(notices,paste("Ignored arguments:",paste(unknown,collapse=", ")))
  if(isTRUE(inputs$want_gory_details)) notices <- c(notices,"want_gory_details is set to FALSE; saved app results exclude package environments.")
  if(kind=="laugh") {
    counts <- inputs$samp_nya
    if(is.null(counts)) {
      p <- inputs$samp_pata; y <- inputs$samp_ny
      numeric_values(p,"samp_pata"); numeric_values(y,"samp_ny")
      if(is.null(p) || is.null(y) || !length(y)) stop("Supply samp_nya, or samp_pata and samp_ny.")
      if(is.null(dim(p))) p <- matrix(rep(p,length(y)),nrow=length(y),byrow=TRUE)
      if(length(dim(p))!=2 || ncol(p)!=a+1 || nrow(p)!=length(y) || any(rowSums(p)<=0)) stop("Sampling margins do not match the year/age grid, or contain an empty composition.")
      counts <- p/rowSums(p) * as.numeric(y)
    }
    numeric_values(counts,"samp_nya")
    if(length(dim(counts))!=2 || ncol(counts)!=a+1) stop("samp_nya must be a year-by-age matrix matching the population grid.")
    if(inherits(counts,"offarray") && (!all(offarray::dimseq(counts)[[1]]==seq_len(nrow(counts))) || !all(offarray::dimseq(counts)[[2]]==0:a))) stop("samp_nya indices must start at year 1 and age 0.")
    counts <- array(as.numeric(counts),dim=unname(dim(counts)))
    totals <- rowSums(counts)
    states <- list(laugh=.design_state(pop,counts,a))
  } else {
    counts <- inputs$nsamp_sya
    if(is.null(counts)) {
      p <- inputs$samp_p_sa; y <- inputs$nsamp_y
      if(is.null(p) || is.null(y) || !length(y)) stop("Supply nsamp_sya, or samp_p_sa and nsamp_y.")
      numeric_values(p,"samp_p_sa"); numeric_values(y,"nsamp_y")
      p <- sex_matrix(p,"samp_p_sa")
      if(ncol(p)!=a+1 || sum(p)<=0) stop("samp_p_sa must match the population grid and have positive total.")
      counts <- array(0,c(2,length(y),a+1))
      for(j in seq_along(y)) counts[,j,] <- p/sum(p)*y[j]
      dimnames(counts) <- list(c("F","M"),NULL,NULL)
    }
    numeric_values(counts,"nsamp_sya")
    if(length(dim(counts))!=3 || dim(counts)[1]!=2 || dim(counts)[3]!=a+1) stop("nsamp_sya must be sex-by-year-by-age, with two sexes and the input population age grid.")
    labs <- dimnames(counts)[[1]]; if(is.null(labs)) labs <- sex_order
    if(!identical(sort(labs),c("F","M"))) stop("nsamp_sya sex labels must be F and M.")
    if(inherits(counts,"offarray") && (!all(offarray::dimseq(counts)[[2]]==seq_len(dim(counts)[2])) || !all(offarray::dimseq(counts)[[3]]==0:a))) stop("nsamp_sya indices must start at year 1 and age 0.")
    counts <- array(as.numeric(counts),dim=unname(dim(counts)))[match(c("F","M"),labs),,,drop=FALSE]
    totals <- apply(counts,c(1,2),sum)
    states <- lapply(1:2,function(i) .design_state(lapply(pop,function(x) as.numeric(x[i,])),matrix(counts[i,,],nrow=dim(counts)[2]),a))
    names(states) <- c("F","M"); dimnames(totals) <- list(c("F","M"),NULL)
  }
  ny <- if(kind=="laugh") length(totals) else ncol(totals)
  if(ny<1 || ny>12 || any(!is.finite(totals)) || !is.finite(sum(totals)) || sum(totals)<=0) stop("Use 1-12 sampling years with a positive total sample count.")
  for(s in states) {
    if(any(s$curves$pmatata>1) || !is.finite(sum(s$curves$nata)) || !is.finite(sum(s$curves$nata*s$curves$wata*s$curves$pmatata)) || sum(s$curves$nata*s$curves$wata*s$curves$pmatata)<=0) stop("Each population needs maturity in [0,1] and positive reproductive output.")
  }
  if(any(totals==0)) notices <- c(notices,"Zero-sample years retain zero effort and an editable flat age composition.")
  list(kind=kind,states=states,counts=totals,settings=settings,notices=notices)
}

#' Convert numeric population schedules and sample counts into one editor state
#' @param pop Three vectors in numbers/weight/maturity order.
#' @param counts Year-by-age counts; zero rows receive flat allocation scores.
#' @param a Plus-group age.
#' @return Raw editor state, preserving the supplied population values.
#' @keywords internal
.design_state <- function(pop,counts,a) {
  scores <- counts
  scores[rowSums(scores)==0,] <- 1
  scores <- scores/rowSums(scores)
  rownames(scores) <- paste("Year",seq_len(nrow(scores)))
  colnames(scores) <- c(as.character(0:(a-1)),paste0(a,"+"))
  list(plus_age=as.integer(a),curves=stats::setNames(lapply(pop,as.numeric),c("nata","wata","pmatata")),sample_scores=scores)
}
