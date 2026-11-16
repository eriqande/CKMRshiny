#' Build paired female/male editors for the two-sex laugh test
#' @return Embedded Shiny UI; the design chooser supplies shared import controls.
#' @export
clt2_app_ui <- function() {
  info <- .laugh_curve_info()
  paired <- function(fun) shiny::div(class="clt2-pair",shiny::fluidRow(lapply(c("F","M"),function(sex)
    shiny::column(6,shiny::h4(if(sex=="F") "Females" else "Males"),fun(sex)))))
  editor <- function(sex,key) shiny::div(id=paste0("clt2_",sex,"_editor_",key),class="laugh-editor",
    `data-key`=key,`data-scope`=paste0("clt2_",sex))
  shiny::tagList(
    shiny::div(class="hero",shiny::div(class="eyebrow","CKMR \u00b7 sampling design"),
      shiny::h1("clt2: a two-sex laugh test"),
      shiny::p("Describe females and males separately, then explore the expected kin counts from your sampling design.")),
    shiny::div(class="section",shiny::h3("What are we describing?"),
      shiny::p("Two sexes in a steady-state, age-structured population. Numbers, biological weight, maturity, and sampling schedules can differ between females and males. The package smooths each population age curve and infers survival separately. These are rough feasibility calculations, not fitted abundance estimates."),
      shiny::p("Both sexes share an age grid and juvenile comparison cutoff. The final age category is a plus group containing all older animals. Expected kin counts are combined across sexes; the sex-specific inputs are used internally.")),
    shiny::div(class="section",shiny::h3("01 \u00b7 Population age grid and abundance"),
      shiny::numericInput("clt2_plus_age","Plus group begins at age",12,min=2,max=200,step=1),
      shiny::actionButton("clt2_apply_age","Set age grid"),
      shiny::p(class="hint","Setting the grid replaces both sexes' curves, including cached sampling years. Undo restores the previous design."),
      paired(function(sex) shiny::tagList(
        shiny::numericInput(paste0("clt2_abundance_",sex),"Set this sex's abundance to",50000,min=1),
        shiny::actionButton(paste0("clt2_scale_",sex),"Set abundance"))),
      shiny::numericInput("clt2_juvenile","Oldest age treated as juvenile for comparisons",2,min=0,max=11,step=1),
      shiny::p("MAX_JUVE_AGE is a comparison and age-precision cutoff, not a maturity threshold. It is shared by the sexes and shown as a red dashed boundary on the editors. Same-birth-cohort HSP comparisons are excluded."),
      shiny::tags$details(shiny::tags$summary("Calculation settings"),
        shiny::numericInput("clt2_catch_min","Youngest age used to fit each catch curve",1,min=0,max=10,step=1),
        shiny::p("MIN_CATCH_CURVE_AGE excludes younger ages from fitting survival. It leaves their population numbers and samples in the design. Leave at least two exact ages before the plus group."),
        shiny::checkboxInput("clt2_after_spawning","Allow a parent sampled in its offspring's birth year",FALSE),
        shiny::p("beqy_POP = TRUE treats those adults as having reproduced before sampling. Otherwise these same-year comparisons are excluded."),
        shiny::numericInput("clt2_etriv","Old-age reproductive-output tolerance (etriv)",.01,min=.000001,max=.99,step=.01),
        shiny::checkboxInput("clt2_extraplus","Extrapolate plus-group abundance from the catch curve",FALSE),
        shiny::p("Normally the supplied plus-group abundance is retained. With extraplus = TRUE it is replaced internally by an extrapolation from younger ages. Both sexes use the larger of the two required age expansions. Final weight and maturity stay constant at older ages; the editors retain the supplied values."))),
    shiny::div(class="laugh-toolbar",shiny::actionButton("clt2_undo","Undo last change"),
      shiny::actionButton("clt2_reset","Reset curves"),shiny::textOutput("clt2_status",inline=TRUE)),
    lapply(c("nata","wata","pmatata"),function(key) shiny::div(class="section laugh-curve-section",
      shiny::h3(info[[key]]$title),
      shiny::p(switch(key,nata="Numbers in each sex, including the entire plus group. Rescaling one sex changes its abundance while retaining its age composition.",
        wata="Biological body weight in consistent units. Reproductive output is proportional to body weight times maturity within each sex. Scaling every weight of one sex equally leaves that sex's relative output unchanged.",
        pmatata="Fraction mature at each age, between zero and one. Female and male maturity schedules can differ; the juvenile comparison cutoff remains shared.")),
      paired(function(sex) editor(sex,key)),
      shiny::actionButton(paste0("clt2_copy_",key),"Copy female curve to males"))),
    shiny::div(class="section",shiny::h3("02 \u00b7 Sampling years, numbers, and age composition"),
      paired(function(sex) shiny::textInput(paste0("clt2_years_",sex),"Samples in consecutive years","250, 250, 250")),
      shiny::p("Enter one expected count per year for each sex, including zeros for gap years. Both lists must have the same length (1-12 years). Fractional expected counts are allowed; one sex can be unsampled."),
      shiny::selectInput("clt2_sampling_year","Year to edit",choices=1:3,selected=1,selectize=FALSE),
      shiny::textOutput("clt2_cache_status"),
      shiny::p("Draw relative allocation scores for this year. Each sex/year curve is normalized separately, then multiplied by that sex's sample count. These are age proportions in the sample, not fractions of the population caught. Other years appear in muted gray behind the active curve."),
      paired(function(sex) shiny::tagList(editor(sex,"samp_pata"),
        shiny::actionButton(paste0("clt2_copy_years_",sex),"Copy this curve to all active years"))),
      shiny::actionButton("clt2_copy_samp_pata","Copy this year's female curve to males"),
      shiny::p(class="hint","Removing years hides their curves without deleting them. Re-adding years restores them; new years copy the selected year. Curve editing does not run the calculation."),
      shiny::plotOutput("clt2_allocation",height="340px"),shiny::div(class="trial-scroll",shiny::tableOutput("clt2_allocation_table"))),
    shiny::div(class="section",shiny::h3("03 \u00b7 Inspect and export inputs"),
      shiny::textOutput("clt2_readiness"),
      shiny::p("nsamp_sya contains sample NUMBERS indexed by sex, consecutive sampling year, and age. CSV rows retain the sex/year/age covariates. RDS files contain ordinary argument lists and can be imported again."),
      shiny::downloadButton("clt2_download_rds","Download inputs (RDS)")," ",
      shiny::downloadButton("clt2_download_r","Download R inputs")," ",
      shiny::downloadButton("clt2_download_csv","Download schedules (CSV)"),
      shiny::tags$details(shiny::tags$summary("Inspect age schedules"),shiny::div(class="trial-scroll",shiny::tableOutput("clt2_inputs"))),
      shiny::tags$details(shiny::tags$summary("Exact argument list"),shiny::verbatimTextOutput("clt2_exact_inputs"))),
    shiny::div(class="section",shiny::h3("04 \u00b7 Expected kin counts"),
      shiny::p("Run microscoping::clt2() explicitly when ready. Each successful calculation saves its exact inputs. Changing the design leaves previous results available."),
      shiny::textOutput("clt2_engine"),shiny::actionButton("clt2_calculate","Calculate expected kin pairs",class="btn-primary"),
      shiny::textOutput("clt2_calculation_error"),
      shiny::selectInput("clt2_result_id","Saved calculation",choices=character(),selectize=FALSE),
      shiny::uiOutput("clt2_result_status"),shiny::plotOutput("clt2_result_plot",height="340px"),
      shiny::div(class="trial-scroll",shiny::tableOutput("clt2_result_table")),
      shiny::textOutput("clt2_ggp_ratio"),
      shiny::p(class="hint","Counts are aggregated across both sexes and all years. Juvenile HSP and GGP counts use cross-cohort comparisons. GGPs can resemble HSPs genetically; these counts do not establish the precision or identifiability of a full CKMR model."),
      shiny::downloadButton("clt2_result_csv","Save expected counts (CSV)")," ",
      shiny::downloadButton("clt2_result_rds","Save calculation and inputs (RDS)"),
      shiny::tags$details(shiny::tags$summary("Population smoothing and plus-group expansion"),
        shiny::textOutput("clt2_preparation_status"),shiny::plotOutput("clt2_preparation_plot",height="340px")),
      shiny::tags$details(shiny::tags$summary("Saved argument list"),shiny::verbatimTextOutput("clt2_saved_inputs"))),
    shiny::div(class="section",shiny::h3("04.1 \u00b7 Kin counts by parent sex"),
      shiny::p("For POPs, female and male refer to the sampled parent. For HSPs, they refer to the unobserved shared parent, regardless of the sampled offspring's sexes. You can therefore predict paternal HSPs even when all sampled offspring are female."),
      shiny::p("These are contributions to the same expected counts shown above, not additional pairs. GGPs remain combined: the grandparent and the intervening parent have distinct sexes to track."),
      shiny::uiOutput("clt2_parent_sex_controls")),
    shiny::div(class="section",shiny::h3("Questions to ask / Things to try"),shiny::tags$ol(
      shiny::tags$li("Start with identical schedules and equal sample counts for the sexes. Compare with the one-sex tab using their combined population and sampling totals."),
      shiny::tags$li("Keep total sampling effort fixed, but move effort between females and males. How do different abundance and maturity schedules change the expected kin counts?"),
      shiny::tags$li("Sample young animals in one year and older animals in another. Compare HSP and GGP predictions."))))
}

#' Bind paired editors and immutable clt2 calculations
#' @param input,output,session Standard Shiny server arguments.
#' @return Reactive state, calculation handles, and an atomic import function.
#' @export
clt2_app_server <- function(input,output,session) {
  defaults <- function(a=12,ny=3) list(F=laugh_default_curves(a,5e4,ny),M=laugh_default_curves(a,5e4,ny),
    counts=matrix(250,2,ny,dimnames=list(c("F","M"),NULL)),
    settings=list(MAX_JUVE_AGE=min(2L,a-1L),MIN_CATCH_CURVE_AGE=min(1L,a-2L),beqy_POP=FALSE,etriv=.01,extraplus=FALSE))
  state <- shiny::reactiveVal(defaults()); history <- shiny::reactiveVal(list()); revision <- shiny::reactiveVal(0L)
  selected_year <- shiny::reactiveVal(1L); status <- shiny::reactiveVal("Draw or import a design, then calculate when ready.")
  counts_error <- shiny::reactiveVal(""); calculations <- shiny::reactiveVal(list()); calculation_error <- shiny::reactiveVal("")
  refresh <- function() {
    s <- state(); a <- s$F$plus_age; cfg <- s$settings
    vals <- list(clt2_plus_age=a,clt2_juvenile=cfg$MAX_JUVE_AGE,clt2_catch_min=cfg$MIN_CATCH_CURVE_AGE,clt2_etriv=cfg$etriv)
    for(id in names(vals)) { shiny::freezeReactiveValue(input,id); shiny::updateNumericInput(session,id,value=vals[[id]],
      max=if(id=="clt2_juvenile") a-1 else if(id=="clt2_catch_min") a-2 else NULL) }
    for(sex in c("F","M")) {
      id <- paste0("clt2_years_",sex); shiny::freezeReactiveValue(input,id)
      shiny::updateTextInput(session,id,value=.laugh_count_text(s$counts[sex,]))
      shiny::updateNumericInput(session,paste0("clt2_abundance_",sex),value=sum(s[[sex]]$curves$nata))
    }
    for(id in c("clt2_after_spawning","clt2_extraplus")) {
      shiny::freezeReactiveValue(input,id)
      shiny::updateCheckboxInput(session,id,value=cfg[[if(id=="clt2_extraplus") "extraplus" else "beqy_POP"]])
    }
    shiny::updateSelectInput(session,"clt2_sampling_year",choices=stats::setNames(seq_len(ncol(s$counts)),paste("Year",seq_len(ncol(s$counts)))),selected=selected_year())
  }
  commit <- function(s,message,sync=FALSE) {
    if(!identical(s,state())) { history(tail(c(history(),list(state())),25)); state(s); revision(revision()+1L) }
    status(message); if(sync) refresh()
  }
  import_design <- function(x) {
    stopifnot(x$kind=="clt2")
    commit(list(F=x$states$F,M=x$states$M,counts=x$counts,settings=x$settings),"Imported design. Undo restores the previous inputs.")
    selected_year(1L); counts_error(""); calculation_error(""); refresh()
  }
  shiny::observeEvent(list(input$clt2_years_F,input$clt2_years_M),{
    if(is.null(input$clt2_years_F) || is.null(input$clt2_years_M)) return()
    tryCatch({
      f <- laugh_parse_years(input$clt2_years_F,TRUE); m <- laugh_parse_years(input$clt2_years_M,TRUE)
      if(length(f)!=length(m)) stop("Give females and males the same number of sampling years, including zeros.")
      s <- state(); next_counts <- rbind(F=f,M=m)
      if(!.laugh_same_counts(s$counts,next_counts)) s$counts <- next_counts
      for(sex in c("F","M")) s[[sex]] <- laugh_ensure_years(s[[sex]],length(f),selected_year())
      selected_year(min(selected_year(),length(f))); counts_error("")
      commit(s,"Sampling totals updated. Cached year curves are retained.")
      shiny::updateSelectInput(session,"clt2_sampling_year",choices=stats::setNames(seq_along(f),paste("Year",seq_along(f))),selected=selected_year())
    },error=function(e) counts_error(conditionMessage(e)))
  },priority=50)
  ids <- c(MAX_JUVE_AGE="clt2_juvenile",MIN_CATCH_CURVE_AGE="clt2_catch_min",beqy_POP="clt2_after_spawning",etriv="clt2_etriv",extraplus="clt2_extraplus")
  for(k in names(ids)) local({ key <- k; id <- ids[[k]]
    shiny::observeEvent(input[[id]],{s <- state(); v <- input[[id]]
      if(key %in% c("MAX_JUVE_AGE","MIN_CATCH_CURVE_AGE") && length(v)==1 && is.finite(v) && v==floor(v)) v <- as.integer(v)
      s$settings[[key]] <- v; commit(s,"Calculation setting updated.")})
  })
  shiny::observeEvent(input$clt2_sampling_year,{
    y <- as.integer(input$clt2_sampling_year)
    if(length(y)==1 && !is.na(y) && y>=1 && y<=ncol(state()$counts)) { selected_year(y); revision(revision()+1L) }
  })
  info <- .laugh_curve_info(); colors <- c(F="#237e83",M="#bd7330")
  send <- function() {
    s <- state()
    for(sex in c("F","M")) session$sendCustomMessage("laugh_editor_state",list(scope=paste0("clt2_",sex),
      plus_age=s[[sex]]$plus_age,juvenile=s$settings$MAX_JUVE_AGE,revision=revision(),
      curves=lapply(names(info),function(key) {
        x <- info[[key]]; x$color <- colors[[sex]]; x$unit <- switch(key,nata="Individuals (this sex)",wata="Biological body weight",pmatata="Proportion mature",samp_pata="Relative allocation score")
        metadata <- if(key=="samp_pata") list(year=selected_year(),series=lapply(seq_len(ncol(s$counts)),function(y)
          list(year=y,label=paste("Year",y),color=colors[[sex]],values=unname(s[[sex]]$sample_scores[y,])))) else list()
        c(list(key=key,values=if(key=="samp_pata") unname(s[[sex]]$sample_scores[selected_year(),]) else s[[sex]]$curves[[key]]),metadata,x)
      })))
  }
  shiny::observe({state();revision();selected_year();send()})
  shiny::observeEvent(input$laugh_editor_ready,send())
  for(ss in c("F","M")) local({sex <- ss
    edit <- function(x,action=FALSE) tryCatch({
      if(is.null(x$revision) || x$revision!=revision()) { send(); return() }
      s <- state(); key <- x$key
      if(!key %in% names(info)) stop("Unknown curve.")
      y <- if(key=="samp_pata") as.integer(x$year) else 1L
      if(length(y)!=1 || is.na(y) || y<1 || y>ncol(s$counts)) stop("Choose an active sampling year.")
      current <- if(key=="samp_pata") s[[sex]]$sample_scores[y,] else s[[sex]]$curves[[key]]
      values <- if(!action) as.numeric(unlist(x$values)) else if(x$action=="preset") laugh_curve_preset(s[[sex]],key,x$preset,y) else {
        if(length(x$factor)!=1 || !x$factor %in% c(.5,2)) stop("Unknown scale action.")
        v <- current*x$factor; if(key=="pmatata") v <- pmin(v,1); v
      }
      values <- laugh_validate_curve(s[[sex]],key,values)
      if(key=="samp_pata") s[[sex]]$sample_scores[y,] <- values else s[[sex]]$curves[[key]] <- values
      commit(s,paste(if(sex=="F")"Female" else "Male",info[[key]]$title,"updated."))
    },error=function(e) { status(conditionMessage(e));send() })
    shiny::observeEvent(input[[paste0("clt2_",sex,"_curve_edit")]],edit(input[[paste0("clt2_",sex,"_curve_edit")]]))
    shiny::observeEvent(input[[paste0("clt2_",sex,"_curve_action")]],edit(input[[paste0("clt2_",sex,"_curve_action")]],TRUE))
    shiny::observeEvent(input[[paste0("clt2_scale_",sex)]],{
      tryCatch({ s <- state(); v <- input[[paste0("clt2_abundance_",sex)]]
        if(length(v)!=1 || !is.finite(v) || v<=0 || sum(s[[sex]]$curves$nata)<=0) stop("Enter a positive abundance and a nonempty curve.")
        s[[sex]]$curves$nata <- s[[sex]]$curves$nata*v/sum(s[[sex]]$curves$nata); commit(s,"Abundance rescaled.")
      },error=function(e) status(conditionMessage(e)))
    })
    shiny::observeEvent(input[[paste0("clt2_copy_years_",sex)]],{
      s <- state(); n <- ncol(s$counts); s[[sex]]$sample_scores[seq_len(n),] <- matrix(rep(s[[sex]]$sample_scores[selected_year(),],n),nrow=n,byrow=TRUE)
      commit(s,"Active years copied. Hidden curves retained; Undo is available.")
    })
  })
  for(kk in names(info)) local({ key <- kk
    shiny::observeEvent(input[[paste0("clt2_copy_",key)]],{
      s <- state()
      if(key=="samp_pata") s$M$sample_scores[selected_year(),] <- s$F$sample_scores[selected_year(),] else s$M$curves[[key]] <- s$F$curves[[key]]
      commit(s,"Female curve copied to males.")
    })
  })
  shiny::observeEvent(input$clt2_apply_age,{
    tryCatch({s <- defaults(input$clt2_plus_age,nrow(state()$F$sample_scores)); s$counts <- state()$counts
      commit(s,"Both age grids reset. Undo restores the previous design.",TRUE)
    },error=function(e) status(conditionMessage(e)))
  })
  shiny::observeEvent(input$clt2_reset,{s <- defaults(state()$F$plus_age,nrow(state()$F$sample_scores)); s$counts <- state()$counts; commit(s,"Curves reset; Undo is available.",TRUE)})
  shiny::observeEvent(input$clt2_undo,{
    h <- history(); if(!length(h)) { status("No earlier change to undo.");return() }
    state(h[[length(h)]]); history(head(h,-1)); selected_year(min(selected_year(),ncol(state()$counts))); counts_error("")
    revision(revision()+1L); status("Previous design restored."); refresh()
  })
  prepared <- shiny::reactive(tryCatch({if(nzchar(counts_error())) stop(counts_error());clt2_function_arguments(state())},error=identity))
  arguments <- shiny::reactive({x <- prepared();shiny::validate(shiny::need(!inherits(x,"error"),if(inherits(x,"error"))conditionMessage(x) else ""));x})
  table_data <- shiny::reactive({
    x <- arguments(); grid <- expand.grid(sex=c("F","M"),year=seq_len(dim(x$nsamp_sya)[2]),age=0:(ncol(x$n_sa)-1))
    idx <- cbind(match(grid$sex,c("F","M")),grid$age+1)
    grid$population <- x$n_sa[idx];grid$weight <- x$w_sa[idx];grid$maturity <- x$pmat_sa[idx];grid$relative_output <- grid$weight*grid$maturity
    grid$age_label <- ifelse(grid$age==ncol(x$n_sa)-1,paste0(grid$age,"+"),as.character(grid$age))
    grid$allocated_samples <- as.numeric(x$nsamp_sya);grid
  })
  output$clt2_status <- shiny::renderText(status())
  output$clt2_cache_status <- shiny::renderText(paste("Editing year",selected_year(),"|",ncol(state()$counts),"active years |",nrow(state()$F$sample_scores)-ncol(state()$counts),"hidden curves per sex"))
  output$clt2_readiness <- shiny::renderText({x <- prepared();if(inherits(x,"error")) paste("Inputs need attention:",conditionMessage(x)) else "Inputs ready. Calculate explicitly to save expected kin counts."})
  output$clt2_allocation_table <- shiny::renderTable(table_data()[,c("sex","year","age_label","allocated_samples")],digits=3)
  output$clt2_inputs <- shiny::renderTable(table_data(),digits=4)
  output$clt2_allocation <- shiny::renderPlot({
    x <- table_data(); a <- state()$F$plus_age; labs <- c(as.character(0:(a-1)),paste0(a,"+"));x$age <- factor(labs[x$age+1],levels=labs)
    ggplot2::ggplot(x,ggplot2::aes(age,factor(year),fill=allocated_samples))+ggplot2::geom_tile(color="white")+
      ggplot2::facet_wrap(~sex,nrow=1,labeller=ggplot2::as_labeller(c(F="Females",M="Males")))+
      ggplot2::scale_fill_gradient(low="#edf3f8",high="#3c6ea8",name="Allocated\nsamples")+
      ggplot2::labs(x="Age category",y="Sampling year",caption="Expected counts, including fractional allocations and zero-sample years.")+.pop_plot_theme()
  })
  output$clt2_exact_inputs <- shiny::renderPrint(arguments())
  output$clt2_download_rds <- shiny::downloadHandler(filename=function()"clt2-inputs.rds",content=function(file)saveRDS(arguments(),file))
  output$clt2_download_r <- shiny::downloadHandler(filename=function()"clt2-inputs.R",content=function(file)writeLines(c("inputs <-",capture.output(dput(arguments())),"# result <- do.call(microscoping::clt2, inputs)"),file))
  output$clt2_download_csv <- shiny::downloadHandler(filename=function()"clt2-age-inputs.csv",content=function(file)utils::write.csv(table_data(),file,row.names=FALSE))
  shiny::observeEvent(input$clt2_calculate,{
    tryCatch({ args <- arguments(); r <- shiny::withProgress(message="Calculating two-sex kin counts",value=.1,clt2_run_design(args))
      saved <- calculations();id <- length(saved)+1L;r$calculation_id <- id;r$editor_state <- state();saved[[as.character(id)]] <- r
      calculations(saved);calculation_error("")
      shiny::updateSelectInput(session,"clt2_result_id",choices=stats::setNames(names(saved),vapply(saved,function(z)paste("Calculation",z$calculation_id,"|",dim(z$function_arguments$nsamp_sya)[2],"years |",signif(sum(z$function_arguments$nsamp_sya),6),"samples"),character(1))),selected=as.character(id))
    },error=function(e)calculation_error(paste("Calculation not completed:",conditionMessage(e),"Previous saved results are unchanged.")))
  })
  selected_calculation <- shiny::reactive({saved <- calculations();shiny::req(length(saved)>0);id <- input$clt2_result_id;if(length(id)!=1 || is.na(id) || !id %in% names(saved)) id <- tail(names(saved),1);saved[[id]]})
  output$clt2_engine <- shiny::renderText(tryCatch(paste("Installed microscoping",laugh_engine_version(),"| native clt2 sex/year/age interface"),error=function(e)conditionMessage(e)))
  output$clt2_calculation_error <- shiny::renderText(calculation_error())
  output$clt2_result_status <- shiny::renderUI({r <- selected_calculation();x <- prepared();same <- !inherits(x,"error") && identical(x[setdiff(names(x),"by_parent_sex")],r$function_arguments[setdiff(names(r$function_arguments),"by_parent_sex")])
    shiny::tagList(shiny::p(if(same)"Saved results match the current inputs." else "Current inputs differ from this saved calculation; calculate again to save new results."),
      shiny::p(paste("microscoping",r$package_version,"|",format(r$created_at))),
      if(length(r$warnings)) shiny::tags$ul(lapply(r$warnings,shiny::tags$li)))})
  output$clt2_result_plot <- shiny::renderPlot(laugh_plot_results(selected_calculation()))
  output$clt2_result_table <- shiny::renderTable(laugh_result_table(selected_calculation())[,c("label","expected","meaning")],digits=5)
  output$clt2_ggp_ratio <- shiny::renderText({x <- selected_calculation()$expected;if(x["HSP_JJ"]>0) paste("Expected GGP/HSP ratio:",signif(100*x["GGP_JJ"]/x["HSP_JJ"],4),"%") else "GGP/HSP ratio is undefined: no juvenile HSPs predicted."})
  output$clt2_preparation_status <- shiny::renderText({p <- selected_calculation()$preparation;paste("Expanded ages: 0 to",max(p$ages),"| Female plus-group survival:",signif(p$sexes$F$plus_survival,5),"| Male:",signif(p$sexes$M$plus_survival,5))})
  output$clt2_preparation_plot <- shiny::renderPlot({
    r <- selected_calculation();p <- r$preparation;args <- r$function_arguments
    model <- do.call(rbind,lapply(c("F","M"),function(sex)data.frame(sex=sex,age=p$ages,n=p$sexes[[sex]]$expanded_nata)))
    raw <- do.call(rbind,lapply(c("F","M"),function(sex)data.frame(sex=sex,age=0:(ncol(args$n_sa)-2),n=args$n_sa[sex,-ncol(args$n_sa)])))
    ggplot2::ggplot(model,ggplot2::aes(age,n,color=sex))+ggplot2::geom_line()+ggplot2::geom_point(data=raw,shape=1)+
      ggplot2::facet_wrap(~sex,scales="free_y",labeller=ggplot2::as_labeller(c(F="Females",M="Males")))+
      ggplot2::scale_color_manual(values=colors,guide="none")+ggplot2::geom_vline(xintercept=ncol(args$n_sa)-1.5,linetype="dashed",color="grey50")+
      ggplot2::labs(x="Age",y="Population numbers",caption="Open points: supplied exact-age numbers. Lines: smoothed, expanded population. Supplied plus-group totals are kept in the saved argument list.")+.pop_plot_theme()
  })
  output$clt2_parent_sex_controls <- shiny::renderUI({
    r <- selected_calculation()
    if(is.null(r$parent_sex_expected)) return(shiny::p(class="hint",
      "This saved calculation has combined counts only. Install the microscoping fork with by_parent_sex support, restart Shiny, and calculate again to add the breakdown. Existing saved results remain available."))
    shiny::tagList(shiny::plotOutput("clt2_parent_sex_plot",height="340px"),
      shiny::div(class="trial-scroll",shiny::tableOutput("clt2_parent_sex_table")),
      shiny::downloadButton("clt2_parent_sex_csv","Save parent-sex counts (CSV)"))
  })
  output$clt2_parent_sex_plot <- shiny::renderPlot(clt2_plot_parent_sex(selected_calculation()))
  output$clt2_parent_sex_table <- shiny::renderTable({
    x <- clt2_parent_sex_table(selected_calculation());shiny::req(!is.null(x))
    data.frame(Output=x$label,Parent=ifelse(x$parent_sex=="F","Female","Male"),Expected=x$expected)
  },digits=5,rownames=FALSE)
  output$clt2_parent_sex_csv <- shiny::downloadHandler(filename=function()"clt2-parent-sex-counts.csv",
    content=function(file) {
      r <- selected_calculation();x <- clt2_parent_sex_table(r);shiny::req(!is.null(x))
      x$calculation_id <- r$calculation_id;x$microscoping_version <- r$package_version
      utils::write.csv(x,file,row.names=FALSE)
    })
  output$clt2_saved_inputs <- shiny::renderPrint(selected_calculation()$function_arguments)
  output$clt2_result_csv <- shiny::downloadHandler(filename=function()"clt2-expected-counts.csv",content=function(file)utils::write.csv(laugh_result_table(selected_calculation()),file,row.names=FALSE))
  output$clt2_result_rds <- shiny::downloadHandler(filename=function()"clt2-calculation.rds",content=function(file)saveRDS(selected_calculation(),file))
  invisible(list(state=state,history=history,revision=revision,selected_year=selected_year,prepared=prepared,arguments=arguments,
    calculations=calculations,selected_calculation=selected_calculation,calculation_error=calculation_error,import_design=import_design))
}
