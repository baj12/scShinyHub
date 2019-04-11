defaultTimeout = 20000
defaultWait = 10
snNr = 1



app <- ShinyDriver$new("../")
app$snapshotInit("mytest3")
Sys.sleep(5)

app$uploadFile(file1 = "../scEx.Rds", wait_=FALSE, values_=FALSE) # <-- This should be the path to the file, relative to the app's tests/ directory
Sys.sleep(20)
cat(file = stderr(), paste("snapshop", snNr, "\n")); snNr=snNr+1;app$snapshot()

app$setInputs(sidebarItemExpanded = "Co-expression", wait_=FALSE, values_=FALSE)
app$setInputs(sideBarID = "coexpressionSelected", wait_=FALSE, values_=FALSE)
app$setInputs(`coE_heatmapSelectedModule-pHeatMapPlot__shinyjquiBookmarkState__resizable` = c(574, 400), allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE)
app$setInputs(`coE_heatmapSelectedModule-pHeatMapPlot_size` = c(574, 400), allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE)
app$setInputs(`coE_heatmapSelectedModule-pHeatMapPlot_is_resizing` = FALSE, allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE)
app$setInputs(`coE_selected-clusterPlot__shinyjquiBookmarkState__resizable` = c(574, 400), allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE)
app$setInputs(`coE_selected-clusterPlot_size` = c(574, 400), allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE)
app$setInputs(`coE_selected-clusterPlot_is_resizing` = FALSE, allowInputNoBinding_ = TRUE, wait_=FALSE, values_=FALSE), wait_=FALSE, values_=FALSE)
app$setInputs(sidebarItemExpanded = "GeneralQC", wait_=FALSE, values_=FALSE), wait_=FALSE, values_=FALSE)
app$setInputs(sideBarID = "gQC_sampleHist", wait_=FALSE, values_=FALSE)
app$snapshot()
cat(file = stderr(), paste("snapshop", snNr, "\n")); snNr=snNr+1;app$snapshot()
app$setInputs(sideBarID = "gQC_variancePC", wait_=FALSE, values_=FALSE)
app$snapshot()
cat(file = stderr(), paste("snapshop", snNr, "\n")); snNr=snNr+1;app$snapshot()
