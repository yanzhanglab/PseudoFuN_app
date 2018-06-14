# if(!require(networkD3)){
#   install.packages('networkD3')
# }

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
                  selected = "CUDAlign54"),
      
      textInput("gene", h3("Enter a gene"),
                value = "ENST00000533288.5"),
      
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
      uiOutput('mytabs')
    )
  )
)

