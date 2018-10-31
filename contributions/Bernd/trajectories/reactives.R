drawTrajectoryHeatmap <- function(x, time, progression_group = NULL, modules = NULL, 
          show_labels_row = FALSE, show_labels_col = FALSE, scale_features = TRUE, 
          ...) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  if (!is.vector(time) || !is.numeric(time)) 
    stop(sQuote("time"), " must be a numeric vector")
  if (nrow(x) != length(time)) 
    stop(sQuote("time"), " must have one value for each row in ", 
         sQuote("x"))
  if ((!is.null(progression_group) && !is.vector(progression_group) && 
       !is.factor(progression_group)) || (!is.null(progression_group) && 
                                          length(progression_group) != nrow(x))) 
    stop(sQuote("progression_group"), " must be a vector or a factor of length nrow(x)")
  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }
  col_ann <- data.frame(row.names = rownames(x), Time = time)
  x_part <- x[order(time), , drop = FALSE]
  if (scale_features) {
    x_part <- scale_quantile(x_part)
  }
  x_part <- t(x_part)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  ann_col <- list(Time = RColorBrewer::brewer.pal(5, "RdGy"))
  if (!is.null(progression_group)) {
    if (!is.factor(progression_group)) 
      progression_group <- factor(progression_group)
    col_ann$Progression <- progression_group
    num_progressions <- length(levels(progression_group))
    progression_cols <- if (num_progressions <= 9) {
      RColorBrewer::brewer.pal(num_progressions, "Set1")
    }
    else {
      gg_color_hue(num_progressions)
    }
    ann_col$Progression <- stats::setNames(progression_cols, 
                                           levels(progression_group))
  }
  labels_row <- if (!show_labels_row) 
    rep("", nrow(x_part))
  else NULL
  labels_col <- if (!show_labels_col) 
    rep("", ncol(x_part))
  else NULL
  if (!is.null(modules)) {
    x_part <- x_part[modules$feature, ]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  }
  else {
    gaps_row <- NULL
    cluster_rows <- T
  }
  # return(list(data = x_part, cluster_cols = F, cluster_rows = cluster_rows, 
  #        annotation_col = col_ann, annotation_colors = ann_col, 
  #        gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col, 
  #        ...))
  # pheatmap::pheatmap(x_part, cluster_cols = F, cluster_rows = cluster_rows,
  #                    annotation_col = col_ann, annotation_colors = ann_col,
  #                    gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col,
  #                    ...)
  list(mat = x_part, cluster_cols = F, cluster_rows = cluster_rows,
                     annotation_col = col_ann, annotation_colors = ann_col,
                     gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col
       )
}

scorpiusSpace <- reactive({
  projections = projections()
  doCalc = input$scorpiusCalc
  dimX = input$dimScorpiusX
  dimY = input$dimScorpiusY
  
  if (!doCalc | is.null(projections)){
    if(DEBUG)cat(file=stderr(), paste("scorpiusSpace:NULL\n"))
    return(NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scorpiusSpace.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scorpiusSpace.RData")
  
  space = projections[,c(dimX, dimY)]
  return(space)
})

scorpiusTrajectory <- reactive({
  space = scorpiusSpace()
  doCalc = input$scorpiusCalc
  if (!doCalc | is.null(space) ){
    if(DEBUG)cat(file=stderr(), paste("scorpiusTrajectory:NULL\n"))
    return(NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scorpiusTrajectory.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scorpiusTrajectory.RData")
  traj <- infer_trajectory(space)
  return(traj)
})

scorpiusExpSel <- reactive({
  gbm_log = gbm_log()
  traj <- scorpiusTrajectory()
  doCalc = input$scorpiusCalc
  
  if (!doCalc | is.null(gbm_log) | is.null(traj)){
    if(DEBUG)cat(file=stderr(), paste("scorpiusExpSel:NULL\n"))
    return(NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scorpiusExpSel.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scorpiusExpSel.RData")
  
  expression = as.matrix(exprs(gbm_log))
  gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  gene_sel <- gimp[1:50,]
  expr_sel <- t(expression)[,gene_sel$gene]
  return(expr_sel)
})

scorpiusModules <- reactive({
  gbm_log = gbm_log()
  projections = projections()
  space <- scorpiusSpace()
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  
  # scorpiusModules = scorpiusModules()
  upI <- updateScorpiusInput() # needed to update input 
  dimX = input$dimScorpiusX
  dimY = input$dimScorpiusY
  # dimCol = input$dimScorpiusCol
  doCalc = input$scorpiusCalc
  
  if (!doCalc | is.null(projections) | is.null(gbm_log)){
    if(DEBUG)cat(file=stderr(), paste("scorpiusModules:NULL\n"))
    return(NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scorpiusModules.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scorpiusModules.RData")
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(gbm_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = T)
  modules = as.data.frame(modules)
  fd = fData(gbm_log)
  modules$symbol = fd[modules$feature,"symbol"]
  rownames(modules) = modules$symbol
  return(modules)
})


updateScorpiusInput <- reactive({
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)){
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "dimScorpiusX",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[1]
  )
  
  updateSelectInput(session, "dimScorpiusY",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[2]
  )
  updateSelectInput(session, "dimScorpiusCol",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[2]
  )
})

# here we define reactive values/variables
# e.g.
# DummyReactive = reactive({
#   # some debugging messages
#   if(DEBUG)cat(file=stderr(), "pca\n")
#   # call dependancies (reactives)
#   gbm_log = gbm_log()
#   # check if they are available
#   if(is.null(gbm_log)){
#     if(DEBUG)cat(file=stderr(), "pca:NULL\n")
#     return(NULL)
#   }
#   if(!is.null(getDefaultReactiveDomain())){
#     showNotification("loading", id="DummyFunc", duration = NULL)
#   }
#   if(DEBUGSAVE) 
#     save(file = "~/scShinyHubDebug/DummyReactive.RData", list = c(ls(),ls(envir = globalenv())))
#   # load(file='~/scShinyHubDebug/DummyReactive.RData')
#   
#   # actual calculation
#   retVal = DummyFunc(gbm_log)
#   
#   if(retVal == 0 & !is.null(getDefaultReactiveDomain()) ){
#     showNotification("Dummy is 0", type="warning", duration=NULL) # has to be removed by use, no removeNotification is following.
#     return(NULL)
#   }
#   if(DEBUG)cat(file=stderr(), "inputData: done\n")
#   if(!is.null(getDefaultReactiveDomain())){
#     removeNotification( id="DummyFunc")
#   }
#   return(retVal)
# })
# 
# # declare heavy calculations
# myHeavyCalculations = list(c("running DummyReactive", "DummyReactive"))


# this will never be executed as it won't be called...

# imageDummyPrecompute <- reactive({
#   if(DEBUG)cat(file=stderr(), "imageDummyPrecompute\n")
#   
#   gbm = gbm()
#   # check if they are available
#   if(is.null(gbm)){
#     if(DEBUG)cat(file=stderr(), "imageDummyPrecompute:NULL\n")
#     return(NULL)
#   }
#   width  <- session$clientData$output_plot_width
#   height <- session$clientData$output_plot_height
#   
#   
#   if(!is.null(getDefaultReactiveDomain())){
#     showNotification("loading", id="DummyFunc", duration = NULL)
#   }
#   # For high-res displays, this will be greater than 1
#   pixelratio <- session$clientData$pixelratio
#   if(is.null(pixelratio)) pixelratio = 1
#   width  <- session$clientData$output_plot_width
#   height <- session$clientData$output_plot_height
#   if(is.null(width)){width=96*7} # 7x7 inch output
#   if(is.null(height)){height=96*7}
#   
#   myPNGwidth <- width/96
#   myPNGheight <- height/96
#   
#   outfile <- paste0(tempdir(),'/dummy.png')
#   if(DEBUG)cat(file=stderr(), paste("output file: ", outfile, "\n"))
#   if(DEBUG)cat(file=stderr(), paste("output file normalized: ", normalizePath(outfile), "\n"))
#   m = data.frame("V1"=colSums(as.matrix(exprs(gbm))))
#   p=ggplot(m, aes(V1)) + geom_bar()
#   ggsave(file=normalizePath(outfile), plot=p, width=myPNGwidth, height=myPNGheight, units="in")
#   
#   if(DEBUG)cat(file=stderr(), "imageDummyPrecompute:done\n")
#   
#   return(list(src = normalizePath(outfile),
#               contentType = 'image/png',
#               width = width,
#               height = height,
#               alt = "Dummy should be here"))
#   
#   
# })
