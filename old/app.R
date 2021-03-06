if(!require(networkD3)){
  install.packages('networkD3')
}

library('shiny')
library('networkD3')
source('helpers.R')

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("PseudoFuN DB Search"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      selectInput("db", 
                  label = "Choose a database",
                  choices = c("BlastDB", 
                              "CUDAlign18",
                              "CUDAlign54", 
                              "CUDAlign135",
                              "CUDAlign198"),
                  selected = "BlastDB"),
      
      textInput("gene", h3("Enter a gene"), 
                value = "PTEN"),
      
      checkboxInput("isgene", "Gene Qeury", value = TRUE),
      
      selectInput("go",
                  label = "GO Analysis",
                  choices = c("Run GO Analysis: Biological Process",
                              "Run GO Analysis: Molecular Function",
                              "Run GO Analysis: Cellular Component",
                              "Do Not Run GO Analysis"),
                  selected = "Do Not Run GO Analysis"),
      
      checkboxInput("inc0", "Include GO terms without any assigned genes", value = FALSE)
    
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      #textOutput(outputId = "selected_db"),
      #textOutput(outputId = "selected_gene"),
      #textOutput(outputId = "runningGO"),
      tabsetPanel(
        tabPanel("Network 1", forceNetworkOutput("net1")),
        tabPanel("Network 2", forceNetworkOutput("net2")),
        tabPanel("Network 3", forceNetworkOutput("net3")),
        tabPanel("Network 4", forceNetworkOutput("net4")),
        tabPanel("GO Analysis", tableOutput("GOtable")),
        tabPanel("Readme", uiOutput('readme'))
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  annot <- reactive({
    get_annot()
  })
  
  dataset <- reactive({
    load_dataset(input$db)
  })

  output$selected_db <- renderText({
    paste("Database: ", input$db)
  })
  
  output$selected_gene <- renderText({
    paste("Gene: ", input$gene)
  })
  
  output$runningGO <- renderText({
    paste(input$go)
  })
  
  # Add for dynamic number of plots

  output$net1 <- renderForceNetwork({
    g <- search2network(input$gene,input$isgene,dataset(),annot(),1);
    forceNetwork(Links = g$links, Nodes=g$nodes,
                 Source = 'source', Target = 'target', NodeID = 'name',
                 Group = 'group',fontSize=16)
  })
  
  output$net2 <- renderForceNetwork({
    g <- search2network(input$gene,input$isgene,dataset(),annot(),2);
    forceNetwork(Links = g$links, Nodes=g$nodes,
                 Source = 'source', Target = 'target', NodeID = 'name',
                 Group = 'group',fontSize=16)
  })
  
  output$net3 <- renderForceNetwork({
    g <- search2network(input$gene,input$isgene,dataset(),annot(),3);
    forceNetwork(Links = g$links, Nodes=g$nodes,
                 Source = 'source', Target = 'target', NodeID = 'name',
                 Group = 'group', fontSize=16)
  })
  
  output$net4 <- renderForceNetwork({
    g <- search2network(input$gene,input$isgene,dataset(),annot(),4);
    forceNetwork(Links = g$links, Nodes=g$nodes,
                 Source = 'source', Target = 'target', NodeID = 'name',
                 Group = 'group', fontSize=16)
  })
  
  GOanalysis <- reactive({
    search2GOtbl(input$gene,input$isgene,input$go,dataset(),annot(),input$inc0)
  })
  
  output$GOtable <- renderTable(GOanalysis())
  
  output$readme <- renderUI({  
    fluidRow(
      column(width = 8, offset = 1,
             includeMarkdown("README.md")
      )
    )
  })
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

