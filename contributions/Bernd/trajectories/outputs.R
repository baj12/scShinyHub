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

output$scorpiusHeatmapPlot <- renderImage({
  if (DEBUG) cat(file=stderr(), paste("scorpiusHeatmapPlot:\n"))
  upI <- updateScorpiusInput() # needed to update input 
  projections = projections()
  # space = scorpiusSpace()
  # traj = scorpiusTrajectory()
  # expr_sel = scorpiusExpSel()
  modules <- scorpiusModules()
  
  # dimX = input$dimScorpiusX
  # dimY = input$dimScorpiusY
  dimCol = input$dimScorpiusCol
  doCalc = input$scorpiusCalc
  
  
  if (!doCalc | is.null(modules) ){
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
  
  
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  draw_trajectory_heatmap(expr_sel, traj$time, projections[,dimCol], modules)
  
  
  # if(DEBUG)cat(file=stderr(), paste("scorpiusHeatmapPlot:done\n"))
  # return(retVal)
})