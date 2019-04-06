library(shinytest)

recordTest(".")

system.time(testApp(".", "mytest.R",compareImages = FALSE))
system.time(testApp(".", "mytest2.R",compareImages = FALSE))

snapshotUpdate(".", "mytest")
