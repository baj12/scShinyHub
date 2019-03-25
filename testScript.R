library(shinytest)

recordTest(".")

testApp(".", "mytest.R")

snapshotUpdate(".", "mytest")
