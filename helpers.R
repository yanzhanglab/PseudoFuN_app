if(!require(topGO)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("topGO")
}
if(!require(org.Hs.eg.db)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}
if(!require(BiocGenerics)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("BiocGenerics")
}
if(!require(biomaRt)){
  install.packages('biomaRt')
}
if(!require(igraph)){
  install.packages('igraph')
}
if(!require(networkD3)){
  install.packages('networkD3')
}

library('biomaRt')
library('igraph')
library('networkD3')
library('topGO')
library('org.Hs.eg.db')
library(reshape2)
library(ggplot2)
library(ggrepel)


letters_only <- function(x) !grepl("[^A-Za-z]", x)

numbers_only <- function(x) !grepl("\\D", x)

load_dataset <- function(dataset){
  message('Loading dataset')
  if(dataset=='BlastDB'){
    load('data/pgAmats_BlastDB.RData')
    return(dataset)
  }else if(dataset=='CUDAlign18'){
    load('data/pgAmats_CUDAlign18.RData')
    return(dataset)
  }else if(dataset=='CUDAlign54'){
    load('data/pgAmats_CUDAlign54.RData')
    return(dataset)
  }else if(dataset=='CUDAlign135'){
    load('data/pgAmats_CUDAlign135.RData')
    return(dataset)
  }else if(dataset=='CUDAlign198'){
    load('data/pgAmats_CUDAlign198.RData')
    return(dataset)
  }else{
    message('Error in dataset name: fun-load_dataset')
    return(NA)
  }
}

# get_annot <- function(){
#   message('Retrieving GENCODE annotation')
#   listMarts(host='dec2016.archive.ensembl.org')
#   ensembl74 <- useMart(host='dec2016.archive.ensembl.org', 
#                        biomart='ENSEMBL_MART_ENSEMBL', 
#                        dataset='hsapiens_gene_ensembl')
#   annot <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'start_position', 'end_position', 'band','gene_biotype'),
#                  mart=ensembl74)
#   return(annot)
# }

map_genes <- function(gene,annot){
  message('Mapping genes back to gencode annotation')
  # Formatting Ensembl gene name
  if(substr(gene,1,4) == 'ENSG'){
    tmp = substr(gene,1,regexpr("\\.",gene)-1)
    if (tmp!=""){
      if (nchar(tmp) == 15){
        return(tmp);
      }else{
        return(NA)
      }
    }else{
      if (nchar(gene) == 15){
        return(gene);
      }else{
        return(NA)
      }
    }
  # Formatting Ensembl transcript name
  }else if(substr(gene,1,4) == 'ENST'){
    tmp = substr(gene,1,regexpr("\\.",gene)-1);
    if (tmp!=""){
      if (nchar(tmp) == 15){
        return(tmp);
      }else{
        return(NA)
      }
    }else{
      if (nchar(gene) == 15){
        return(gene);
      }else{
        return(NA)
      }
    }
  # Formatting entrez gene ID
  }else if (numbers_only(gene)){
    tmp = annot[which(annot[,"entrezgene"]==gene),c("gene_biotype","ensembl_gene_id","ensembl_transcript_id")];
    tmp1 = tmp[,"ensembl_gene_id"]
    tmp1[grep("pseudogene",tmp[,"gene_biotype"])] = tmp[grep("pseudogene",tmp[,"gene_biotype"]),"ensembl_transcript_id"]
    tmp1 = unique(tmp1);
    if (length(tmp1) == 0){
      return(NA)
    }else{
      return(tmp1)
    }
  # Formatting Hugo gene symbols
  }else{
    tmp = annot[which(annot[,"hgnc_symbol"]==gene),c("gene_biotype","ensembl_gene_id","ensembl_transcript_id")];
    tmp1 = tmp[,"ensembl_gene_id"]
    tmp1[grep("pseudogene",tmp[,"gene_biotype"])] = tmp[grep("pseudogene",tmp[,"gene_biotype"]),"ensembl_transcript_id"]
    tmp1 = unique(tmp1);
    if (length(tmp1) == 0){
      return(NA)
    }else{
      return(tmp1)
    }
  }
}

in_pgAmat <- function(pgAmat,gene){
  if(substr(gene,1,4) == 'ENSG'){
    if(gene %in% names(pgAmat)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else if(substr(gene,1,4) == 'ENST'){
    if(length(grep(gene,names(pgAmat)))>0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    message("Error in query gene")
    return(NA)
  }
}

find_pgAmats <- function(genes,dataset){
  message('Finding pgAmats containing gene')
  tmp = list()
  if(all(is.na(genes)==TRUE)){
    return(NA)
  }else{
    for(gene in genes){
      tmp[[gene]] = names(dataset)[which(unlist(lapply(dataset,in_pgAmat, gene=gene)))]
    }
  }
  tmp = unlist(tmp)
  tmp = unique(tmp)
  return(tmp)
}

int_graph <- function(pgAmat){
  message('Generating and plotting network figure')
  names <- colnames(pgAmat);
  i=0;
  for(name in names){
    i=i+1;
    if(name %in% annot[,"ensembl_gene_id"]){
      message('Gene with gene symbol')
      tmp = annot[which(annot[,"ensembl_gene_id"]==name),"hgnc_symbol"]
      if(!is.na(tmp)){
        names[i] = paste0("Gene: ", tmp, ": ", name)
      }else{
        names[i] = paste0("Gene: ",name)
      }
    }else if(substr(name,1,regexpr("\\.",name)-1) %in% annot[,"ensembl_transcript_id"]){
      message('Pseudogene with gene symbol')
      tmp = annot[which(annot[,"ensembl_transcript_id"]==substr(name,1,regexpr("\\.",name)-1)),"hgnc_symbol"]
      if(!is.na(tmp)){
        names[i] = paste0("Pseudogene: ", tmp, ": ", name)
      }else{
        names[i] = paste0("Pseudogene: ",name)
      }
    }else if(substr(name,1,4)=="ENST"){
      message('Pseudogene without gene symbol')
      names[i] = paste0("Pseudogene: ",name)
    }else{
      message('Gene without gene symbol')
      names[i] = paste0("Gene: ",name)
    }
  }
  row.names(pgAmat) = names;
  colnames(pgAmat) = names;
  g = graph_from_adjacency_matrix(log2(pgAmat),mode="undirected",weighted=TRUE);
  E(g)$weight = max(E(g)$weight)-E(g)$weight;
  mstg = mst(g,weights=E(g)$weight);
  E(g)$weight = max(E(g)$weight)+E(g)$weight;
  E(mstg)$weight = max(E(mstg)$weight)+E(mstg)$weight
  wc <- cluster_walktrap(mstg)
  members <- membership(wc)
  nd3g = igraph_to_networkD3(mstg,group=members)
  return(nd3g)
  # plotted D3 graph now returns D3 graph
  #forceNetwork(Links = nd3g$links, Nodes=nd3g$nodes,
  #             Source = 'source', Target = 'target', NodeID = 'name',
  #             Group = 'group')
  #plot.igraph(mstg,vertex.label=V(mstg)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(mstg)$weight)
}

search2network <- function(gene,dataset,annot,idx){
  genes <- map_genes(gene,annot);
  pgAmats <- find_pgAmats(genes,dataset);
  if(idx<=length(pgAmats)){
    return(int_graph(as.matrix(dataset[[pgAmats[idx]]])))
  }else{
    return(NA) 
  }
}

search2adjmat <- function(gene,dataset,annot,idx){
  genes <- map_genes(gene,annot);
  pgAmats <- find_pgAmats(genes,dataset);
  if(idx<=length(pgAmats)){
    return(as.matrix(dataset[[pgAmats[idx]]]))
  }else{
    return(NA) 
  }
}

num_networks <- function(gene,dataset,annot){
  genes <- map_genes(gene,annot)
  pgAmats <- find_pgAmats(genes,dataset)
  return(length(pgAmats))
}

expr_analysis <-function(g,tcga_cancer_type){
  exprG = readRDS(paste0("data/TCGA_rds_gene/",tcga_cancer_type,".rds"))
  exprP = readRDS(paste0("data/dreamBase_rds_pseudo/TCGA_",tcga_cancer_type,".rds"))
  miR_cor = readRDS(paste0("data/miR-gene_rds/",tcga_cancer_type,"-TP.rds"))
  
  nodes <- g$nodes$name
  node.mat <- matrix(unlist(strsplit(as.character(nodes),': ')),ncol=3,byrow=TRUE)
  Eg = exprG[unique(node.mat[node.mat[,1]=="Gene",2],na.rm=TRUE),]
  grnames = row.names(Eg)[row.names(Eg)!="NA"]
  if(length(grnames)>0)
    Eg = Eg[grnames,]
  Eg = apply(Eg,2,as.numeric)
  if(length(grnames)<2){Eg = t(as.matrix(Eg))}
  row.names(Eg) = grnames
  Ep = exprP[unique(node.mat[node.mat[,1]=="Pseudogene",2],na.rm=TRUE),]
  prnames = row.names(Ep)[row.names(Ep)!="NA"]
  if(length(prnames)>0){
    Ep = Ep[prnames,]
    Ep = apply(Ep,2,as.numeric)
    if(length(prnames)<2){Ep = t(as.matrix(Ep))}
    row.names(Ep) = prnames
  }
  cnames = intersect(colnames(Eg),colnames(Ep))
  if(length(grnames)>0){
    Eg = Eg[,cnames]
    if(length(grnames)<2){Eg = t(as.matrix(Eg))}
    row.names(Eg) = grnames
  }
  if(length(prnames)>0){
    Ep = Ep[,cnames]
    if(length(prnames)<2){Ep = t(as.matrix(Ep))}
    row.names(Ep) = prnames
  }
  if(length(prnames)>0 & length(grnames)>0){
    Etcga <<- rbind(Eg,Ep)
  }else if(length(grnames)>0){
    Etcga <<- Eg
  }else if(length(prnames)>0){
    Etcga <<- Ep
  }else{
    Etcga <<- NULL
  }
  if(!is.null(Etcga)){
    Ctumor = as.numeric(substr(cnames,14,16))<10;
    Et = Etcga[,Ctumor]
    En = Etcga[,!Ctumor]
    Ct <<- cor(log2(t(Et)+1))
    Cn <<- cor(log2(t(En)+1))
    Ct[is.na(Ct)] <<- 0;
    Cn[is.na(Cn)] <<- 0;
    E.df = melt(Etcga)
    E.df$Var1 = as.character(E.df$Var1)
    Disease = ifelse(as.numeric(substr(E.df$Var2,14,16))<10,"tumor","normal")
    genes = unique(E.df$Var1)
    if (length(unique(Disease)) > 1){
      smartModal(error=F, title = "Warning", content = "There's only one group (tumor) in current TCGA dataset. Normal heatmap may not show up.")
      for(gene in genes){
        tmp = t.test(E.df[E.df$Var1==gene,"value"]~Disease[E.df$Var1==gene])
        E.df[E.df$Var1==gene,"Var1"] = paste0(gene,"\n(pval=",formatC(tmp$p.value, format = "e", digits = 2),")")
      }
    }
    fig_expr_box <<- ggplot(aes(y = log2(value+1),x = as.factor(Var1), fill = Disease), data = E.df) +
                            geom_boxplot() +labs(x="",y="Norm Counts (log2)") +
                            theme_bw() +
                            theme(legend.title = element_text(size = 12),
                                  legend.text = element_text(size = 10),
                                  axis.title=element_text(size=14,face="bold"),
                                  axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                                  axis.text.y = element_text(size=12))
                            
  }else{
    Ct <<- NULL;
    Cn <<- NULL;
    fig_expr_box <<- NULL;
  }
  
  miR_gene_cor <<- miR_cor[miR_cor$Gene %in% genes,]
  if(dim(miR_gene_cor)[1]>0){
    fig_miR_scatter <<- ggplot(aes(y=Corr,x=as.factor(Gene),label=Mir),data=miR_gene_cor) +
                              geom_point(size = 3) + geom_text_repel() +
                              labs(x="", y="Correlation") +
                              scale_y_continuous(trans = "reverse") +
                              theme_bw() +
                              theme(legend.title = element_text(size = 12),
                                    legend.text = element_text(size = 10),
                                    axis.title=element_text(size=14,face="bold"),
                                    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
                                    axis.text.y = element_text(size=12))
  }else{
    miR_gene_cor <<- NULL;
  }
}

search2GOtbl <- function(gene,go,dataset,annot,inc0,
                         run.ks, run.ks.elim){
  if(go != "Do Not Run GO Analysis"){
    if(go=="Run GO Analysis: Biological Process"){
      ontol="BP"; top_nodes=10000
    } else if(go =="Run GO Analysis: Molecular Function"){
      ontol="MF"; top_nodes=2500
    } else{ontol="CC"; top_nodes=1000}
    message('Generating GO table')
    genes <- map_genes(gene,annot);
    pgAmats <- find_pgAmats(genes,dataset);
    gene_set = c()
    for(pgAmat in pgAmats){
      gene_set <- rbind(names(dataset[[pgAmat]]))
    }
    genes_all <- factor(as.integer(annot[,'ensembl_gene_id'] %in% gene_set))
    names(genes_all) <- annot[,'ensembl_gene_id'];
    #genes_all = genes_all[!is.na(names(genes_all))]
    #geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
    #idx = c(which(genes_all==1),sample(which(genes_all==0), 10000, replace = FALSE))
    #genes_all = genes_all[idx]
    GOdata <- new("topGOdata",ontology = ontol,
                allGenes = genes_all,
                geneSel=function(p) p == 1,
                description ="inNetwork",
                annot=annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    if(run.ks){resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")}
    if(run.ks.elim){resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")}
    
    if(run.ks & run.ks.elim){
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       classicKS = resultKS, elimKS = resultKS.elim,
                       orderBy = "elimKS", ranksOf = "elimKS", topNodes = top_nodes)
    } else if(run.ks & !run.ks.elim){
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS,
                         orderBy = "classicKS", ranksOf = "classicKS", topNodes = top_nodes)
    } else if(!run.ks & run.ks.elim){
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "elimKS", topNodes = top_nodes)
    } else {
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = top_nodes)
    }
    
    if(inc0){
      return(allRes)
    }else{
        return(allRes[allRes[,'Significant']>0,])
    }
  }else{
    return("No GO Analysis")
  }
}

# smartModal is for showing processing windows. To close a smartModal, a "removeModal()" should followed.
smartModal <- function(error=c(T,F), title = "Title", content = "Content"){
  if(error){
    showModal(modalDialog(
      title = title, footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p(content),
          style = "margin: auto; text-align: center"
      )
    ))
  }
  else{
    showModal(modalDialog(
      title = title, footer = NULL,
      div(class = "busy",
          p(content),
          img(src="https://cdn.dribbble.com/users/503653/screenshots/3143656/fluid-loader.gif"),
          style = "margin: auto; text-align: center"
      )
    ))
  }
}

