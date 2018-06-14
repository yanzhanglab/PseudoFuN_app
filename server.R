library('shiny')
library('networkD3')
source('helpers.R')
library(DT)
# Define server logic required to draw a histogram ----
function(input, output, session) {
  num_tabs <<- 0
  tabs.list <<- NULL
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
  
  observeEvent(input$action1, {
    if (length(tabs.list) > 0){
      for (i in 1:length(tabs.list)){
        print("remove tabs")
        print(tabs.list[i])
        removeTab(inputId = "tabs", target = tabs.list[i])
      }
      tabs.list <<- NULL
    }
    for (i in 1:num_tabs)
    print("rendering network panels...")
    smartModal(error=F, title = "Processing", content = "We are processing your request ...")
    t <- try(num_tabs <<- num_networks(input$gene,input$isgene,dataset(),annot))
    if("try-error" %in% class(t)) {
      removeModal()
      print("Error occured")
      smartModal(error=T, title = "Error occured", content = "Error occured. Try to enter a valid gene.")
      return()
    }
    else{
      GOanalysis <- search2GOtbl(input$gene,input$isgene,input$go,dataset(),annot,input$inc0,
                                 input$run.ks, input$run.ks.elim)
      output$GOtable <- DT::renderDataTable({
        GOanalysis
      },selection="none",options=list(searching=F, ordering=F))#,extensions = 'Responsive'
      
      for (i in 1:num_tabs){
        print(i)
        tabs.list <<- c(tabs.list, paste0('Network ',i))
        appendTab(inputId = "tabs",
                  tab = tabPanel(paste0('Network ',i),
                                 h2(paste0('Network ',i), style="color: STEELBLUE; font-size: 22px"),
                                 forceNetworkOutput(paste0('net',i)))
        )
        removeModal()
        output[[paste0('net',i)]] <- renderForceNetwork({
          smartModal(error=F, title = "Processing", content = "Initializing Force Directed Networks ...")
          t2 <- try(g <- search2network(input$gene,input$isgene,dataset(),annot,i))
          if("try-error" %in% class(t2)) {
            removeModal()
            print("Error occured")
            smartModal(error=T, title = "Error occured", content = "Searching Network Failed. Try to enter a valid gene.")
            return()
          }
          removeModal()
          forceNetwork(Links = g$links, Nodes=g$nodes,
                       Source = 'source', Target = 'target', NodeID = 'name',
                       Group = 'group', fontSize = 16)
        })
      }
    }
  })
  

}