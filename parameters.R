# tabList = list(
#   tabItem("umiHist",
#           c(tags$h3("Histogram of UMI counts"),
#           fluidRow(column(
#             10, offset = 1,
#             plotOutput('plotUmiHist') %>% withSpinner()
#           )))
#   ))
if(DEBUG)cat(file=stderr(), paste("parameters:", length(allTabs)," ", "\n"))

normaliztionChoices = list(raw = "rawNormalization")
# parameterContributions = list()
parFiles = dir(path = "contributions", pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
for(fp in parFiles){
  myNormalizationChoices = c()
  source(fp, local = TRUE)
  save(file = "normalizationCH.RData", list=ls())
  # load(file = "normalizationCH.RData")
  if (length(myNormalizationChoices) > 0){
    for (li in 1:length(myNormalizationChoices)){
      liVal = myNormalizationChoices[[li]]
      if(length(liVal)>0){
        if(DEBUG)cat(file=stderr(), paste("normalization Choice: ", liVal, "\n"))
        oldNames = names(normaliztionChoices)
        normaliztionChoices[[length(normaliztionChoices) + 1]] = liVal
        names(normaliztionChoices) = c(oldNames, names(myNormalizationChoices)[li])
      }
    }
  }
}

# here we add content to the page on the rigth (main visualization window)
allTabs[[length(allTabs) + 1]] = list(
  tabItem("normalizations",list(
    tags$h3("Parameters for normalization to be used"),
    fluidRow(column(10,
                    radioButtons(inputId = "normalizationRadioButton",
                                 label    = "Normalization to use",
                                 choices  = normaliztionChoices,
                                 width = '100%')
                    # 10, offset = 1,
                    # plotOutput('plotUmiHist') %>% withSpinner()
    )),
    fluidRow(column(10, verbatimTextOutput("normalizationRadioButtonValue"))),
    wellPanel(
      # This outputs the dynamic UI component
      uiOutput("normalizationsParametersDynamic")
    )
  )
  ,
  tableSelectionUi("normalizationResult")
  
  )
)