library('shiny')

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
                  label = "Choose a variable to display",
                  choices = c("BlastDB", 
                              "CUDAlign18",
                              "CUDAlign54", 
                              "CUDAlign135",
                              "CUDAlign198"),
                  selected = "CUDAlign54"),
      
      textInput("gene", h3("Text input"), 
                value = "Enter Gene")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      textOutput(outputId = "selected_db"),
      textOutput(outputId = "selected_gene")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$selected_db <- renderText({
    paste("You have selected", input$db)
  })
  
  output$selected_gene <- renderText({
    paste("You have selected", input$gene)
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

