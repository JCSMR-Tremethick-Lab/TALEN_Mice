#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(DT)

# Define UI for application that draws a boxplot
ui <- shinyUI(
  fluidPage(
   # Application title
   titlePanel("TALEN Mice RNA-Seq data"),
     
   # Sidebar with selection box using dynamic search terms
   sidebarLayout(
      sidebarPanel(
        selectInput("gene_name", 
                    label = h4("Select gene"),
                    choices = unique(obkt$ext_gene),
                    selectize = T,
                    selected = "H2afb2"),
        width = 2
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("boxplot"),
         fluidRow(DT::dataTableOutput("table")),
         width = 6
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
   dataset.kt <- reactive({
     obkt[obkt$ext_gene == input$gene_name,]
   })
   
   output$boxplot <- renderPlot({
    p <- ggplot(dataset.kt(), aes(factor(condition), tpm)) + geom_boxplot()
    print(p)    
   })
   
   output$table <- DT::renderDataTable(DT::datatable({
     data <- obst[obst$ext_gene == input$gene_name,]})
   )
})

# Run the application 
shinyApp(ui = ui, server = server)

