library(shinytest)

recordTest(".")

system.time(testApp(".", "mytest.R",compareImages = FALSE))

snapshotUpdate(".", "mytest")
