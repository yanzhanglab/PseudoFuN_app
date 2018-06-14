library('shiny')
library('networkD3')
source('helpers.R')

# Define server logic required to draw a histogram ----
function(input, output, session) {
  num_tabs <<- NULL
  load("./annot.Rdata")
  annot <<- annot
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
    })
  })
  
  observeEvent(input$gene,{
    print("rendering network panels...")
    t <- try(num_tabs <<- num_networks(input$gene,input$isgene,dataset(),annot))
    if("try-error" %in% class(t)) {
      removeModal()
      print("Error occured")
      smartModal(error=T, title = "Error occured", content = "Error occured. Try to enter a valid gene.")
      return()
    }
    else{
      for (i in 1:num_tabs){
        print(i)
        appendTab(inputId = "tabs",
                  tab = tabPanel(paste0('Network ',i),
                                 h2(paste0('Network ',i), style="color: STEELBLUE; font-size: 22px"),
                                 forceNetworkOutput(paste0('net',i)))
        )
        output[[paste0('net',i)]] <- renderForceNetwork({
          smartModal(error=F, title = "Processing", content = "Initializing Force Directed Networks ...")
          g <- search2network(input$gene,input$isgene,dataset(),annot,i);
          removeModal()
          forceNetwork(Links = g$links, Nodes=g$nodes,
                       Source = 'source', Target = 'target', NodeID = 'name',
                       Group = 'group')
        })
      }
    }
  })
  

}