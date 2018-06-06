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
                value = "Enter Gene")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      textOutput(outputId = "selected_db"),
      textOutput(outputId = "selected_gene"),
      plotOutput(outputId = 'network')
      
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
  
  output$network <- renderPlot({
    search2plotgen(input$gene,dataset(),annot())
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

