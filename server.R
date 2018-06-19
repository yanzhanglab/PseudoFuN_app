library('shiny')
library('networkD3')
source('helpers.R')
source('circos.R')
library(DT)
library(circlize)

# Define server logic required to draw a histogram ----
function(input, output, session) {
  num_tabs <<- 0
  g <<- NULL
  tabs.list <<- NULL
  load("./data/annot.Rdata")
  annot <<- annot
  
  observeEvent(input$action1, {
    dataset <<- load_dataset(input$db)
    if (length(tabs.list) > 0){
      for (i in 1:length(tabs.list)){
        print("remove tabs")
        print(tabs.list[i])
        removeTab(inputId = "tabs", target = tabs.list[i])
      }
      tabs.list <<- NULL
    }
    print("rendering network panels...")
    smartModal(error=F, title = "Processing", content = "We are processing your request ...")
    t <- try(num_tabs <<- num_networks(input$gene,input$isgene,dataset,annot))
    if("try-error" %in% class(t)) {
      removeModal()
      print("Error occured")
      print(dataset)
      smartModal(error=T, title = "Error occured", content = "Error occured. Try to enter a valid gene.")
      return()
    }
    else{
      if (input$go > 0){
        removeModal()
        smartModal(error=F, title = "Processing", content = "We are processing your GO analysis (may take a few minutes)")
        GOanalysis <- search2GOtbl(input$gene,input$isgene,input$go,dataset,annot,input$inc0,
                                   input$run.ks, input$run.ks.elim)
        removeModal()
        output$GOtable <- DT::renderDataTable({
          GOanalysis
        },selection="none",options=list(searching=F, ordering=F))#,extensions = 'Responsive'
        output$download_go <- downloadHandler(
          filename = function() {
             name = "GO_result.csv"
          },
          content = function(file) {
             write.table(GOanalysis, file = file, append = FALSE, quote = TRUE, sep = ',',
                         eol = "\r\n", na = "NA", dec = ".", row.names = F,
                         col.names = T, qmethod = c("escape", "double"),
                         fileEncoding = "")
        })
        session$sendCustomMessage("download_go","-")
      }
      for (i in 1:num_tabs){
        # print(i)
        tabs.list <<- c(tabs.list, paste0('Network ',i))
        appendTab(inputId = "tabs",
                  tab = tabPanel(paste0('Network ',i),
                                 h2(paste0('Network ',i),
                                 style="color: STEELBLUE; font-size: 22px"),
                                 forceNetworkOutput(paste0('net',i)),
                                 actionButton(paste0('circos',i), "Circos Plot", style="color: WHITE; background-color: #FFC300"),
                                 br(),br()
                  )
        )
      }
      removeModal()
      Map(function(i) {
        print(paste0('net',i))
        output[[paste0('net',i)]] <- renderForceNetwork({
          print("render force map")
          smartModal(error=F, title = "Processing", content = "Initializing Force Directed Networks ...")
          t2 <- try(g[[i]] <<- search2network(input$gene,input$isgene,dataset,annot,i))
          if("try-error" %in% class(t2)) {
            removeModal()
            print("Error occured")
            smartModal(error=T, title = "Error occured", content = "Searching Network Failed. Try to enter a valid gene.")
            return()
          }
          removeModal()
          targetposition = match(input$gene, sub(".*: ", "", g[[i]]$nodes$name))
          nodesize = rep(1,length(g[[i]]$nodes$name))
          nodesize[targetposition] = 50
          g[[i]]$nodes$size = nodesize
          forceNetwork(Links = g[[i]]$links, Nodes=g[[i]]$nodes,
                       Source = 'source', Target = 'target', NodeID = 'name',
                       Nodesize = 'size',
                       Group = 'group', fontSize = 16,  fontFamily = 'sans')
        })
      },
      1:num_tabs)
      
      
      Map(function(i) {
        print(num_tabs)
        observeEvent(input[[paste0('circos',i)]],{
          print(paste0('circos',i))
          print(g[[i]]$nodes$name)
        })
      },
      1:num_tabs)
      
    }
  })
  
  
  observeEvent(input$button_circos, {
    load("./data/UCSC_hg19_refGene_20180330.Rdata") # varname: hg19
    load("./data/UCSC_hg38_refGene_20180330.Rdata") # varname: hg38
    genes_str = sub(".*: ", "", g$nodes$name)
    # genes_str <- c("LOC102725121", "FAM138A", "RIMS2", "LINC01128", "MMP23A", "ULK4P1")
    hg19 <- data.frame(cbind(rownames(hg19), hg19, hg19[6]-hg19[5]))
    hg38 <- data.frame(cbind(rownames(hg38), hg38, hg38[6]-hg38[5]))
    colnames(hg38) = c("id","","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID","","","","length")
    colnames(hg19) = c("id","","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID","","","","length")
    hg19.ring <- hg19[!grepl("_", hg19$chrom),] # remove undefined chromosome
    hg38.ring <- hg38[!grepl("_", hg38$chrom),]
    hg19.ring <- hg19.ring[!grepl("chrM", hg19.ring$chrom),]
    hg38.ring <- hg38.ring[!grepl("chrM", hg38.ring$chrom),]
    hg19.matched <- hg19.ring[match(genes_str, hg19.ring$alignID, nomatch = 0), ]
    hg38.matched <- hg38.ring[match(genes_str, hg38.ring$alignID, nomatch = 0), ]
    hg19.ring.lengthsum <- aggregate(hg19.ring["length"],hg19.ring["chrom"],sum)
    hg38.ring.lengthsum <- aggregate(hg38.ring["length"],hg38.ring["chrom"],sum)
    output$circos_plot_ui_hg38 <- renderUI({
      plotOutput("circos_plot_component_hg38", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_ui_hg19 <- renderUI({
      plotOutput("circos_plot_component_hg19", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_component_hg38 <- renderPlot({
      factors_count = as.data.frame(hg38.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genome (GRCh38/hg38)",
                       T,T)
    })
    
    output$circos_plot_component_hg19 <- renderPlot({
      factors_count = as.data.frame(hg19.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                       input$circos_param_genelink,
                       input$circos_param_genesymbol)
      })
      removeModal()
    })
  

}
