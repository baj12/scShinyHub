# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$scropius_trajectory_plot <- renderPlot({
  projections = projections()
  upI <- updateScorpiusInput() # needed to update input 
  dimX = input$dimScorpiusX
  dimY = input$dimScorpiusY
  dimCol = input$dimScorpiusCol
  doCalc = input$scorpiusCalc
  
  if(!doCalc |  is.null(projections)  ){
    return(NULL)
  }
  
  if(DEBUG)cat(file=stderr(), paste("scropius_trajectory_plot:\n"))
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scropius_trajectory_plot.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scropius_trajectory_plot.RData")
  space = projections[,c(dimX, dimY)]
  traj <- infer_trajectory(space)
  draw_trajectory_plot(space, progression_group = projections[,dimCol], path = traj$path)
  
})

callModule(tableSelectionServer, "scorpiusTableMod", scorpiusModules)
# selected clusters heatmap module

scorpiusHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file=stderr(), paste("scorpiusHeatmapPlot:\n"))
  upI <- updateScorpiusInput() # needed to update input 
  projections = projections()
  traj = scorpiusTrajectory()
  expr_sel = scorpiusExpSel()
  modules <- scorpiusModules()

  dimCol = input$dimScorpiusCol
  doCalc = input$scorpiusCalc

  
  if (!doCalc | is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(traj) ){
    if(DEBUG)cat(file=stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/scorpiusHeatmapPlot.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/scorpiusHeatmapPlot.RData")
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(gbm_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  pixelratio <- session$clientData$pixelratio
  if(is.null(pixelratio)) pixelratio = 1
  width  <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  if (is.null(width)) {width = 96*7} # 7x7 inch output
  if (is.null(height)) {height = 96*7}
  
  outfile <- paste0(tempdir(), '/heatmapScorpius', base::sample(1:10000, 1), '.png')
  cat(file = stderr(), paste("saving to: ", outfile, '\n'))
  
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal = drawTrajectoryHeatmap(expr_sel, traj$time, projections[,dimCol], modules,
                                   filename = normalizePath(outfile))

  retVal
  
  # pheatmap(retVal$data, Rowv = retVal$cluster_rows, Colv = retVal$cluster_cols,
  #           showticklabels = c(retVal$show_colnames,T), row_side_colors = retVal$annotation_row, 
  #           col_side_colors = retVal$annotation_col)
  
  # # if(DEBUG)cat(file=stderr(), paste("scorpiusHeatmapPlot:done\n"))
  # return(list(src = normalizePath(outfile),
  #             contentType = 'image/png',
  #             width = width,
  #             height = height,
  #             alt = "heatmap should be here"))
  # 
})

callModule(
  pHeatMapModule, 
  "scorpiusHeatmapPlotModule", 
  scorpiusHeatmapPlotReactive
)
