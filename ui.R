# if(!require(networkD3)){
#   install.packages('networkD3')
# }
library('shiny')
library('networkD3')
source('helpers.R')
library(shinyWidgets)

# Define UI for app that draws a histogram ----


navbarPage(  
  # App title ----
  # titlePanel("PseudoFuN DB Search"),
  title=div(a(img(src="http://icons.iconarchive.com/icons/graphicloads/100-flat/256/zoom-search-2-icon.png",
                  height = 30,
                  style = "margin:0px 0px; padding-bottom: 5px"),
              "PseudoFuN DB Search", href="")
           ),
  
  tabPanel("Search Engine",
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(width = 3,
        
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
        
        prettyCheckbox(inputId = "isgene", label = "Gene Query", value = TRUE, icon = icon("check")),
        

        selectInput("go",
                    label = "GO Analysis",
                    choices = c("Run GO Analysis: Biological Process",
                                "Run GO Analysis: Molecular Function",
                                "Run GO Analysis: Cellular Component",
                                "Do Not Run GO Analysis"),
                    selected = "Do Not Run GO Analysis"),
        
        # checkboxInput("inc0", "Include GO terms without any assigned genes", value = FALSE)
        prettyCheckbox(inputId = "inc0", label = "Include GO terms without any assigned genes", value = FALSE, icon = icon("check"))
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        tabsetPanel(id = "tabs",
          tabPanel("GO Analysis", "This is the hello tab")
        )
      )
    )
  ), # end of tabPanel "Search Engine"
  
  tabPanel("Read Me"
  ),
  tabPanel("About"
  )
)# end of navbar page


