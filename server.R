library('shiny')
library('networkD3')
source('helpers.R')
source('circos.R')
library(DT)
library(circlize)
library(ggplot2)
library(reshape2)


# Define server logic required to draw a histogram ----
function(input, output, session) {
  num_tabs <<- NULL
  g <<- NULL
  adjmat <<- list()
  g.circos <<- NULL
  tabs.list <<- NULL
  active_net <<- NULL # the current active network, integer value.
  load("data/annot.Rdata")
  load("data/UCSC_hg19_refGene_20180330.Rdata") # varname: hg19
  load("data/UCSC_hg38_refGene_20180330.Rdata") # varname: hg38
  annot <<- annot
  hg19 <<- hg19
  hg38 <<- hg38
  current.gene <<- NULL
  current.db <<- NULL
  current.cancer <<- NULL
  
  observeEvent(input$action1, {
    dataset <<- load_dataset(input$db)
    current.gene <<- toupper(input$gene)
    current.db <<- input$db
    output$Circos_plot_gene_db_name <- renderUI({
      h2(sprintf("Gene: %s; Database: %s.", current.gene, current.db),
         style="color: grey; font-size: 18px; font-family: Courier")
    })
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
    t <- try(num_tabs <<- num_networks(input$gene,dataset,annot))
    if("try-error" %in% class(t)) {
      removeModal()
      print("Error occured")
      smartModal(error=T, title = "Error occured", content = "Error occured. Try to enter a valid gene.")
      return()
    }
    if(num_tabs == 0) {
      removeModal()
      print("Error occured")
      sendSweetAlert(session, title = "Error occured", sprintf("%s does not return any network in %s. Please try another Gene/Database.",input$gene,input$db), type = "error",
                     btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    else{
      if (input$go > 0){
        removeModal()
        smartModal(error=F, title = "Processing", content = "We are processing your GO analysis (may take a few minutes)")
        GOanalysis <- search2GOtbl(input$gene,input$go,dataset,annot,input$inc0,
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
                                 h2("Please enable pop-ups from this website to view gene/pseudogene link out resources.",style="color: STEELBLUE; font-size: 14px"),
                                 forceNetworkOutput(paste0('net',i)),
                                 actionButton(paste0('circos_',current.gene,'_',current.db,'_',i),
                                              "Circos Plot", style="color: WHITE; background-color: #FFC300"),
                                 actionButton(paste0('TCGA_expression_',current.gene,'_',current.db,'_',i),
                                              "TCGA Expression", style="color: WHITE; background-color: #ff4ce4"),
                                 downloadButton(paste0('download_net',i), 'Download Adjacency Matrix'),
                                 helpText("The elements in adjacency matrix indicating the similarity between each genes."),
                                 br(),br()
                  ),
                  session = getDefaultReactiveDomain()
        )
      }
      removeModal()
      
      Map(function(i) {
        print(paste0('net',i))
        output[[paste0('net',i)]] <- renderForceNetwork({
          print("render force map")
          smartModal(error=F, title = "Processing", content = "Initializing Force Directed Networks ...")
          t2 <- try(g[[i]] <<- search2network(current.gene,dataset,annot,i))
          if("try-error" %in% class(t2)) {
            removeModal()
            print("Error occured")
            smartModal(error=T, title = "Error occured", content = "Searching Network Failed. Try to enter a valid gene.")
            return()
          }
          t2 <- try(adjmat[[i]] <<- search2adjmat(current.gene,dataset,annot,i))
          if("try-error" %in% class(t2)) {
            removeModal()
            print("Error occured")
            smartModal(error=T, title = "Error occured", content = "Searching Network Failed. Try to enter a valid gene.")
            return()
          }
          removeModal()
          mapped_genes <- map_genes(current.gene,annot);
          message(current.gene)
          targetposition = match(mapped_genes, substr(sub(".*: ", "", g[[i]]$nodes$name),1,15))
          nodesize = rep(1,length(g[[i]]$nodes$name))
          nodesize[targetposition] = 50
          g[[i]]$nodes$size = nodesize
          MyClickScript <- ''
          MyClickScript <- 'var split = d.name.split(": ");
          var genename = split[2].substr(0,15);
          if (genename.length > 0){
          window.open("https://www.genecards.org/Search/Keyword?queryString="+genename, "_blank");
          window.open("https://useast.ensembl.org/Homo_sapiens/Gene/TranscriptComparison?db=core;g="+genename, "_blank");
          }
          else{
          alert(d.name + " doesn\'t contain any gene symbol!");
          }'
          forceNetwork(Links = g[[i]]$links, Nodes=g[[i]]$nodes,
                       Source = 'source', Target = 'target', NodeID = 'name',
                       Nodesize = 'size', opacity = 1,
                       Group = 'group', fontSize = 16,  fontFamily = 'sans',
                       clickAction = MyClickScript)
        })
      },
      1:num_tabs)
      
      Map(function(i) {
        output[[paste0('download_net',i)]] <- downloadHandler(
          filename = function() {
            name = sprintf("%s_net_%d_adjmat.csv", input$gene, i)
          },
          content = function(file) {
            write.table(adjmat[[i]], file = file, append = FALSE, quote = TRUE, sep = ',',
                        eol = "\r\n", na = "NA", dec = ".", row.names = T,
                        col.names = NA, qmethod = c("escape", "double"),
                        fileEncoding = "")
          })
      },
      1:num_tabs)
      
      
      Map(function(i) {
        observeEvent(input[[paste0('TCGA_expression_',current.gene,'_',current.db,'_',i)]],{
          active_net <<- i
          message(sprintf("active net: %d", active_net))
          current.cancer <<- input$TCGA_cancer
          output$TCGA_Expression_gene_db_cancer_name <- renderUI({
            h2(sprintf("Gene: %s; Database: %s; Cancer: %s; Network: %d.", current.gene, current.db, current.cancer, active_net),
               style="color: grey; font-size: 18px; font-family: Courier")
          })
          smartModal(error=F, title = "Calculating", content = "Calculating TCGA Expression ...")
          print(sprintf("Gene: %s; Database: %s; Cancer: %s.", current.gene, current.db, current.cancer))
          expr_analysis(g[[i]], input$TCGA_cancer, session, input$miRNA_pred_targ)
          if (dim(Etcga)[1]==1){
            output$normal_heatmap = NULL
            output$normal_heatmap <- renderPlot({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("Correlation matrix cannot be calculated.\n",
                                           "Only 1 gene/pseudogene expression was found."), 
                   cex = 1.6, col = "black")
            })
            output$tumor_heatmap = NULL
            output$tumor_heatmap <- renderPlot({
              par(mar = c(0,0,0,0))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, paste("Correlation matrix cannot be calculated.\n",
                                           "Only 1 gene/pseudogene expression was found."), 
                                          cex = 1.6, col = "black")
            })
          }else{
            if (is.null(Cn)){
              output$normal_heatmap = NULL
              output$normal_heatmap <- renderPlot({
                par(mar = c(0,0,0,0))
                plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                text(x = 0.5, y = 0.5, paste("All disease type are tumor.\n",
                                             "No normal sample available."), 
                                              cex = 1.6, col = "black")
              })
            }
            else{
              output$normal_heatmap <- renderPlot({
                order = heatmap(Cn)$rowInd
                Cn.reorder = Cn[,order]
                Cn.reorder = Cn.reorder[order,]
                diag(Cn.reorder) = NA
                Cn.melt = melt(Cn.reorder)
                ggplot(Cn.melt, aes(Var1, Var2, value))+
                  geom_tile(aes(fill = value), color = "white") +
                  scale_fill_gradient2(low = "blue", mid = "white", high = "yellow") +
                  ylab("") +
                  xlab("") +
                  theme(legend.title = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        plot.title = element_text(size=16),
                        axis.title=element_text(size=14,face="bold"),
                        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                        axis.text.y = element_text(size=12)) +
                  labs(fill = "Expression level")
              })
            }
            output$tumor_heatmap <- renderPlot({
              order = heatmap(Ct)$rowInd
              Ct.reorder = Ct[,order]
              Ct.reorder = Ct.reorder[order,]
              diag(Ct.reorder) = NA
              Ct.melt = melt(Ct.reorder)
              ggplot(Ct.melt, aes(Var1, Var2, value))+
                geom_tile(aes(fill = value), color = "white") +
                scale_fill_gradient2(low = "blue", mid = "white", high = "yellow") +
                ylab("") +
                xlab("") +
                theme(legend.title = element_text(size = 12),
                      legend.text = element_text(size = 10),
                      plot.title = element_text(size=16),
                      axis.title=element_text(size=14,face="bold"),
                      axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                      axis.text.y = element_text(size=12)) +
                labs(fill = "Expression level")
            })
          }
          output$pseudo_boxplot <- renderPlot({fig_expr_box})
          output$correlation_plot <- renderPlot({fig_miR_scatter})
          #DGE table
          dgetable <- readRDS(sprintf("./data/Pseudogene_rds_DGE/TCGA_%s_%s.rds", "BRCA", current.db))
          dgetable.cutoff <<- dgetable[dgetable$FDR < 10^(-input$DGE_cutoff_value),]
          dgetable.cutoff <<- data.frame(apply(dgetable.cutoff, 2, function(x) unlist(x)))
          output$DGEtable <- DT::renderDataTable({dgetable.cutoff}, selection="none", escape = F,
                                                 options=list(searching=T, pageLength = 20))
          session$sendCustomMessage("myCallbackHandler", "tab_TCGA_Expression")
        })
      },
      1:num_tabs)
      
      
      Map(function(i) {
        observeEvent(input[[paste0('circos_',current.gene,'_',current.db,'_',i)]],{
          smartModal(error=F, title = "Calculating", content = "Generating Circos Plot ...")
          g.circos <<- g[[i]]
          print(paste0('circos',i))
          print(g.circos$nodes$name)
          genes_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 2)
          genestype_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 1)
          pesudogenes_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 3)
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
          print("Plot HG38 Circos")
          output$circos_plot_component_hg38 <- renderPlot({
            factors_count = as.data.frame(hg38.ring.lengthsum)
            factors = factor(factors_count[,1], levels = factors_count[,1])
            xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
            rownames(xlim) = factors_count[,1]
            BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
            circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genome (GRCh38/hg38)",
                             input$circos_param_genelink,
                             input$circos_param_genesymbol,
                             input$font.scale,
                             input$link.width,
                             input$color.picker)
          })
          print("Plot HG19 Circos")
          output$circos_plot_component_hg19 <- renderPlot({
            factors_count = as.data.frame(hg19.ring.lengthsum)
            factors = factor(factors_count[,1], levels = factors_count[,1])
            xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
            rownames(xlim) = factors_count[,1]
            BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
            circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                             input$circos_param_genelink,
                             input$circos_param_genesymbol,
                             input$font.scale,
                             input$link.width,
                             input$color.picker)
          })
          print("Plot Finished.")
          removeModal()
          session$sendCustomMessage("myCallbackHandler", "tab_circos")
        })
      },
      1:num_tabs)
    }
  })
  
  
  output$download_TCGA_exp <- downloadHandler(
    filename = 'Gene_Mir_expr_results.zip',
    content = function(fname) {
      separator = ','
      fs <- c('Etcga.csv', 'miR_gene_cor.csv', sprintf('DGEtable_cutoff=1e-%d.csv', input$DGE_cutoff_value))
      write.table(Etcga, file = fs[1], sep = separator, col.names = NA)
      write.table(miR_gene_cor, file = fs[2], sep = separator, col.names = NA)
      write.table(dgetable.cutoff, file = fs[3], sep = separator, col.names = NA)
      zip(zipfile=fname, files=fs)
      if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
    },
    contentType = "application/zip"
  )
  
  observeEvent(input$action3,{
    smartModal(error=F, title = "Updating", content = "Updating TCGA Expression ...")
    current.gene <<- input$gene
    current.db <<- input$db
    current.cancer <<- input$TCGA_cancer
    output$TCGA_Expression_gene_db_cancer_name <- renderUI({
      h2(sprintf("Gene: %s; Database: %s; Cancer: %s; Network: %d.", current.gene, current.db, current.cancer, active_net),
         style="color: grey; font-size: 18px; font-family: Courier")
    })
    expr_analysis(g[[active_net]], current.cancer, session, input$miRNA_pred_targ)
    if (dim(Etcga)[1]==1){
      output$normal_heatmap = NULL
      output$normal_heatmap <- renderPlot({
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Correlation matrix cannot be calculated.\n",
                                     "Only 1 gene/pseudogene expression was found."), 
             cex = 1.6, col = "black")
      })
      output$tumor_heatmap = NULL
      output$tumor_heatmap <- renderPlot({
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Correlation matrix cannot be calculated.\n",
                                     "Only 1 gene/pseudogene expression was found."), 
             cex = 1.6, col = "black")
      })
    }else{
      if (is.null(Cn)){
        output$normal_heatmap = NULL
        output$normal_heatmap <- renderPlot({
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          text(x = 0.5, y = 0.5, paste("All disease type are tumor.\n",
                                       "No normal sample available."), 
               cex = 1.6, col = "black")
        })
      }else{
        output$normal_heatmap <- renderPlot({
          order = heatmap(Cn)$rowInd
          Cn.reorder = Cn[,order]
          Cn.reorder = Cn.reorder[order,]
          diag(Cn.reorder) = NA
          Cn.melt = melt(Cn.reorder)
          ggplot(Cn.melt, aes(Var1, Var2, value))+
            geom_tile(aes(fill = value), color = "white") +
            scale_fill_gradient2(low = "blue", mid = "white", high = "yellow") +
            ylab("") +
            xlab("") +
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  plot.title = element_text(size=16),
                  axis.title=element_text(size=14,face="bold"),
                  axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                  axis.text.y = element_text(size=12)) +
            labs(fill = "Expression level")
        })
      }
      output$tumor_heatmap <- renderPlot({
        order = heatmap(Ct)$rowInd
        Ct.reorder = Ct[,order]
        Ct.reorder = Ct.reorder[order,]
        diag(Ct.reorder) = NA
        Ct.melt = melt(Ct.reorder)
        ggplot(Ct.melt, aes(Var1, Var2, value))+
          geom_tile(aes(fill = value), color = "white") +
          scale_fill_gradient2(low = "blue", mid = "white", high = "yellow") +
          ylab("") +
          xlab("") +
          theme(legend.title = element_text(size = 12),
                legend.text = element_text(size = 10),
                plot.title = element_text(size=16),
                axis.title=element_text(size=14,face="bold"),
                axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                axis.text.y = element_text(size=12)) +
          labs(fill = "Expression level")
      })
    }
    output$pseudo_boxplot <- renderPlot({fig_expr_box})
    output$correlation_plot <- renderPlot({fig_miR_scatter})
    #DGE table
    dgetable <- readRDS(sprintf("./data/Pseudogene_rds_DGE/TCGA_%s_%s.rds", "BRCA", current.db))
    dgetable.cutoff <- dgetable[dgetable$FDR < 10^(-input$DGE_cutoff_value),]
    dgetable.cutoff <<- data.frame(apply(dgetable.cutoff, 2, function(x) unlist(x)))
    output$DGEtable <- DT::renderDataTable({dgetable.cutoff}, selection="none", escape = F,
                                           options=list(searching=T, pageLength = 20))
  })
  
  
  observeEvent(input$action2,{
    print(g.circos$nodes$name)
    genes_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 2)
    genestype_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 1)
    pesudogenes_str = sapply(strsplit(as.character(g.circos$nodes$name), ": "), "[[", 3)
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
    print("Plot HG38 Circos")
    output$circos_plot_component_hg38 <- renderPlot({
      factors_count = as.data.frame(hg38.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genome (GRCh38/hg38)",
                       input$circos_param_genelink,
                       input$circos_param_genesymbol,
                       input$font.scale,
                       input$link.width,
                       input$color.picker)
    })
    print("Plot HG19 Circos")
    output$circos_plot_component_hg19 <- renderPlot({
      factors_count = as.data.frame(hg19.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                       input$circos_param_genelink,
                       input$circos_param_genesymbol,
                       input$font.scale,
                       input$link.width,
                       input$color.picker)
    })
    print("Plot Finished.")
    removeModal()
  })
  
  

}
