
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyBS)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(edgeR)
library(pheatmap)
library(threejs)
library(shinyTree)
library(shinycssloaders)

# this is where the general tabs are defined:
if(file.exists('defaultValues.R')){
  source('defaultValues.R')
}
source('tabs.R')


# general tabs
allTabs = list(       
  inputTab,
  geneSelectionTab,
  cellSelectionTab
)

# Basic menu Items
allMenus = list(
  menuItem("input", tabName = "input", icon = icon("dashboard")),
  menuItem("Cell selection", tabName = "cellSelection", icon = icon("dashboard")),
  menuItem("Gene selection", tabName = "geneSelection", icon = icon("dashboard"))
)



# parse all ui.R files under contributions to include in application
uiFiles = dir(path = "contributions", pattern = "ui.R", full.names = TRUE, recursive = TRUE)
for(fp in uiFiles){
  menuList = list()
  tabList = list()
  source(fp, local = TRUE)
  
  for (li in menuList){
    if(length(li)>0){
      # if(DEBUG)cat(file=stderr(), paste("menuList:", length(allMenus)," ", li$children, "\n"))
      allMenus[[length(allMenus) + 1 ]] = li
    }
  }
  for (li in tabList){
    if(length(li)>0){
      # if(DEBUG)cat(file=stderr(), paste(li$children[[1]], "\n"))
      allTabs[[length(allTabs) + 1]] = li
    }
  }
  
}



shinyUI(
  dashboardPage(
    dashboardHeader(title = "scShinyHub"),
    dashboardSidebar(
      sidebarMenu(
        allMenus
        
      ),
      tipify(checkboxInput("noStats", "don't display stats", FALSE),
             "check this if you are working on the cell/gene selection to avoid certain calculations"),
      
      tipify(htmlOutput('summaryStatsSideBar'),
            "<h3>Data summary</h3> <ul><li>medium UMI: shows how many genes are  expressed in log2 space of normalized data</li> </ul> ", "right"),
      downloadButton("report", "Generate report"),
      actionButton('goCalc', 'Force Calculations'),
      # bookmarkButton(id = "bookmark1"),
      tipify(downloadButton("countscsv", "Download counts.csv"),
             '<h3>download current count data as CSV file</h3>'),
      checkboxInput("DEBUGSAVE", "Save for DEBUG", FALSE),
      verbatimTextOutput("DEBUGSAVEstring")
      
    ), # dashboard side bar
    dashboardBody(
      bsAlert("alert"),
      tags$div(
        allTabs,
        class = "tab-content"
      )
    ) # dashboard body
  ) # main dashboard
)


