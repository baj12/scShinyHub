app <- ShinyDriver$new("../")
app$snapshotInit("mytest3")

app$uploadFile(file1 = "../scEx.Rds", wait_=FALSE, values_=FALSE) # <-- This should be the path to the file, relative to the app's tests/ directory
app$setInputs(sidebarCollapsed = FALSE, wait_=FALSE, values_=FALSE)
Sys.sleep(240)
app$snapshot()
app$setInputs(goCalc = "click", wait_=FALSE, values_=FALSE)
Sys.sleep(600)
app$snapshot()
