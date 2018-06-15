
library('shiny')
library('networkD3')
source('helpers.R')
library('shinyWidgets')

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
        

        awesomeRadio("go",
                    label = "GO Analysis",
                    choices = list("Do Not Run GO Analysis"=0,
                                   "Run GO Analysis: Biological Process"=1,
                                   "Run GO Analysis: Molecular Function"=2,
                                   "Run GO Analysis: Cellular Component"=3),
                              # c("Run GO Analysis: Biological Process",
                              #   "Run GO Analysis: Molecular Function",
                              #   "Run GO Analysis: Cellular Component",
                              #   "Do Not Run GO Analysis"),
                    status = "primary"),
        conditionalPanel(condition="input.go > 0",
                         helpText("Fisher classic will be run by default. KS Classic and KS elim are optional and may cost longer time (few minutes)."),
                         prettyCheckbox(inputId = "run.ks", label = "Run KS Classic", value = F, icon = icon("check")),
                         prettyCheckbox(inputId = "run.ks.elim", label = "Run KS elim", value = F, icon = icon("check")),
                         prettyCheckbox(inputId = "inc0", label = "Include GO terms without any assigned genes", value = FALSE, icon = icon("check"))
                         ),
        # checkboxInput("inc0", "Include GO terms without any assigned genes", value = FALSE)
        
        actionButton("action1", "Confirm and Run", style="color: WHITE; background-color: DODGERBLUE")
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        tabsetPanel(id = "tabs",
          tabPanel("GO Analysis",
                   h4("GO Analysis", style="color: STEELBLUE"),
                   DT::dataTableOutput("GOtable")
                   )
        )
      )
    )
  ), # end of tabPanel "Search Engine"
  
  tabPanel("Read Me",
             fluidRow(
               column(width = 8, offset = 1,
                      includeMarkdown("README.md")
               )
             )
          ),
  
  tabPanel("About"
  ),
  tags$head(tags$script(HTML("document.title = 'PseudoFuN DB Search';"))) # rename the title by JS
)# end of navbar page


