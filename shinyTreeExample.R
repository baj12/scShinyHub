# library(devtools) install_github('trestletech/shinyTree')
library(shiny)
library(shinyTree)
server <- shinyServer(function(input, output, session) {
  output$tree <- renderTree({
    list(`I lorem impsum` = list(
      `I.1 lorem impsum` = structure(list(`I.1.1 lorem impsum` = "1", `I.1.2 lorem impsum` = "2"), stselected = TRUE),
      `I.2 lorem impsum` = structure(list(`I.2.1 lorem impsum` = "3"), stselected = TRUE)
    ))
  })
})
ui <- shinyUI(shiny::fluidPage(h4("Shiny hierarchical checkbox"), shinyTree("tree", checkbox = TRUE)))
shinyApp(ui, server)
