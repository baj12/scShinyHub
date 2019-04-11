library(shinytest)

recordTest(".")

system.time(testApp(".", "mytest.R",compareImages = TRUE))
system.time(testApp(".", "mytest2.R",compareImages = TRUE))

snapshotUpdate(".", "mytest2")
