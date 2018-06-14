library('shiny')
library('networkD3')
source('helpers.R')

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  load("./annot.Rdata")
  
  dataset <- reactive({
    load_dataset(input$db)
  })
  
  GOanalysis <- reactive({
    search2GOtbl(input$gene,input$isgene,input$go,dataset(),annot,input$inc0)
  })
  
  GOtable <- reactive({
    renderTable({
      GOanalysis()
      })
    })
  
  readme <- reactive({
    renderUI({  
      fluidRow(
        column(width = 8, offset = 1,
              includeMarkdown("README.md")
        )
      )
    })})
  
  output$mytabs = renderUI({
    message('Starting mytabs renderUI')
    num_tabs <- num_networks(input$gene,input$isgene,dataset(),annot)
    message('Calculated number of networks')
    # myTabs = lapply(c(paste('Network', 1: num_tabs),'GO Analysis','README'), tabPanel)
    # do.call(tabsetPanel, myTabs)})
    myTabs = vector('list',num_tabs)
    i=1
    while(i<=num_tabs){
      helpText(sprintf("value %d", i))
      myTabs[[i]] <- tabPanel(paste0('Network ',i),forceNetworkOutput(search2network(input$gene,input$isgene,dataset(),annot,i)))
      i=i+1
    }
    # myTabs[[1]] <- tabPanel("GO Analysis", NULL)
    # myTabs[[2]] <- tabPanel("README", NULL)
    
    # save(myTabs,file = "~/Desktop/myTabs.Rdata")
    message(paste0('Number of tabs: ', length(myTabs)))
    message(class(myTabs))
    do.call(tabsetPanel, myTabs)
    message('finished mytabs renderUI')
  })
}