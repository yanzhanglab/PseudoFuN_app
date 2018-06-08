library('shiny')
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
                  selected = "CUDAlign54"),
      
      textInput("gene", h3("Enter a gene"), 
                value = "ENSG00000172236"),
      
      selectInput("go",
                  label = "GO Analysis",
                  choices = c("Run GO Analysis",
                              "Do not run GO Analysis"),
                  selected = "Do not run GO Analysis"),
      
      checkboxInput("inc0", "Include 0 Significant Genes in GO Analysis", value = FALSE)
    
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      textOutput(outputId = "selected_db"),
      textOutput(outputId = "selected_gene"),
      textOutput(outputId = "runningGO"),
      plotOutput(outputId = 'network'),
      tableOutput(outputId = 'GOtable')
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
  
  output$network <- renderPlot({
    search2plotgen(input$gene,dataset(),annot())
  })
  
  GOanalysis <- reactive({
    search2GOtbl(input$gene,input$go,dataset(),annot(),input$inc0)
  })
  
  output$GOtable <- renderTable(GOanalysis())
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

