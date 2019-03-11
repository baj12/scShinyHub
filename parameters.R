# tabList = list(
#   tabItem("umiHist",
#           c(tags$h3("Histogram of UMI counts"),
#           fluidRow(column(
#             10, offset = 1,
#             plotOutput('plotUmiHist') %>% withSpinner()
#           )))
#   ))
if (DEBUG) cat(file = stderr(), paste("parameters:", length(allTabs), " ", "\n"))

normaliztionChoices <- list(raw = "rawNormalization")
# parameterContributions = list()
parFiles <- dir(path = "contributions", pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
save(file = "normalizationCH.RData", list = ls())
# load(file = "normalizationCH.RData")
for (fp in parFiles) {
  if (DEBUG) {
    cat(file = stderr(), paste(fp, "\n"))
  }

  myNormalizationChoices <- c()
  source(fp, local = TRUE)
  if (length(myNormalizationChoices) > 0) {
    for (li in 1:length(myNormalizationChoices)) {
      liVal <- myNormalizationChoices[[li]]
      if (length(liVal) > 0) {
        if (DEBUG) cat(file = stderr(), paste("normalization Choice: ", liVal, "\n"))
        oldNames <- names(normaliztionChoices)
        normaliztionChoices[[length(normaliztionChoices) + 1]] <- liVal
        names(normaliztionChoices) <- c(oldNames, names(myNormalizationChoices)[li])
      }
    }
  }
  if (DEBUG) {
    cat(file = stderr(), paste("end:", fp, "\n"))
    cat(file = stderr(), paste("end:", normaliztionChoices, "\n"))
  }
}

# here we add content to the page on the rigth (main visualization window)
allTabs[[length(allTabs) + 1]] <- list(
  tabItem(
    "normalizations",
    list(
      tags$h3("Parameters for normalization to be used"),
      tags$p("scShinyHub generally uses normalized data (unless stated otherwise). Here, the specific method can be set. Raw means that no normalization will be performed."),
      tags$p("A table containing the first 20 cells and the normalized values is shown."),
      fluidRow(column(
        10,
        radioButtons(
          inputId = "normalizationRadioButton",
          label = "Normalization to use",
          choices = normaliztionChoices,
          selected = "gbm_logNormalization",
          width = "100%"
        )
        # 10, offset = 1,
        # plotOutput('plotUmiHist') %>% withSpinner()
      )),
      fluidRow(column(10, verbatimTextOutput("normalizationRadioButtonValue"))),
      wellPanel(
        # This outputs the dynamic UI component
        uiOutput("normalizationsParametersDynamic")
      )
    ),
    tableSelectionUi("normalizationResult")
  )
)
if (DEBUG) {
  cat(file = stderr(), paste("end: parameters.R\n"))
}
