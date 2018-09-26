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
parameterContributions = list()
parFiles = dir(path = "contributions", pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
for(fp in parFiles){
  myNormalizationChoices = c()
  source(fp, local = TRUE)
  
  for (li in myNormalizationChoices){
    if(length(li)>0){
      if(DEBUG)cat(file=stderr(), paste("normalization Choice: ", li, "\n"))
      normaliztionChoices[[length(normaliztionChoices) + 1]] = li
    }
  }
}

# here we add content to the page on the rigth (main visualization window)
allTabs[[length(allTabs) + 1]] = list(
  tabItem("normalizations",list(
    tags$h3("Histogram of UMI counts"),
    fluidRow(column(10,
                    radioButtons(inputId = "normalizationRadioButton",
                                label    = "Normalization to use",
                                choices  = normaliztionChoices,
                                selected = 1)
      # 10, offset = 1,
      # plotOutput('plotUmiHist') %>% withSpinner()
    )),
    fluidRow(column(3, verbatimTextOutput("normalizationRadioButtonValue"))),
    wellPanel(
      # This outputs the dynamic UI component
      uiOutput("normalizationsParametersDynamic")
    )
    )
  )
)