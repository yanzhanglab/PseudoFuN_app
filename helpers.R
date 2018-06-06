if(FALSE){
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("clusterProfiler")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("DOSE")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("RDAVIDWebService")
  source("https://bioconductor.org/biocLite.R")
  biocLite("topGO")
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}


library('biomaRt')
library('igraph')
library('topGO')
#library('RDAVIDWebService')
#library('DOSE')
#library('clusterProfiler')

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

get_annot <- function(){
  message('Retrieving GENCODE annotation')
  listMarts(host='dec2016.archive.ensembl.org')
  ensembl74 <- useMart(host='dec2016.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
  annot <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'start_position', 'end_position', 'band'),
                 mart=ensembl74)
  return(annot)
}

map_gene <- function(gene, annot){
  message('Mapping genes back to gencode annotation')
  # Formatting Ensembl gene name
  if(substr(gene,1,4) == 'ENSG'){
    tmp = substr(gene,1,regexpr("\\.",gene))
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
    tmp = substr(gene,1,regexpr("\\.",gene));
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
    #annot <- get_annot();
    tmp = annot[which(annot[,"entrezgene"]==gene),"ensembl_gene_id"];
    if (length(tmp) == 0){
      return(NA)
    }else{
      return(tmp)
    }
  # Formatting Hugo gene symbols
  }else{
    #annot <- get_annot();
    tmp = annot[which(annot[,"hgnc_symbol"]==gene),"ensembl_gene_id"];
    if (length(tmp) == 0){
      return(NA)
    }else{
      return(tmp)
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
    if(length(grep(gene,names(dataset[['pgAmat16.csv']])))>0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    message("Error in query gene")
    return(NA)
  }
}

find_pgAmat <- function(genes,dataset){
  message('Finding pgAmat containing gene')
  tmp = c()
  if(is.na(genes)){
    return(NA)
  }else if(genes == ""){
    return(NA)
  }else{
    for(gene in genes){
      tmp[gene] = names(dataset)[which(unlist(lapply(dataset,in_pgAmat, gene=gene)))]
    }
  }
  tmp = tmp[which(tmp!="")];
  return(tmp)
}

int_graph <- function(pgAmat){
  message('Generating and plotting network figure')
  g = graph_from_adjacency_matrix(log2(pgAmat),mode="undirected",weighted=TRUE);
  E(g)$weight = max(E(g)$weight)-E(g)$weight;
  mstg = mst(g,weights=E(g)$weight);
  E(g)$weight = max(E(g)$weight)+E(g)$weight;
  E(mstg)$weight = max(E(mstg)$weight)+E(mstg)$weight
  plot.igraph(mstg,vertex.label=V(mstg)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(mstg)$weight)
}

search2plotgen <- function(gene,dataset,annot){
  genes <- map_gene(gene,annot);
  pgAmats <- find_pgAmat(genes,dataset);
  for(pgAmat in pgAmats){
    int_graph(as.matrix(dataset[[pgAmats]]))
  }
}

search2GOtbl <- function(gene,go,dataset,annot){
  if(go == "Run GO Analysis"){
  message('Generating GO table')
  library('org.Hs.eg.db')
  genes <- map_gene(gene,annot);
  pgAmats <- find_pgAmat(genes,dataset);
  gene_set <- names(dataset[[pgAmats]])
  message(gene_set)
  genes_all <- ifelse(annot[,'ensembl_gene_id'] %in% gene_set,1,0)
  names(genes_all) <- annot[,'ensembl_gene_id'];
  idx = c(which(genes_all==1),sample(which(genes_all==0), 20000, replace = FALSE))
  genes_all = genes_all[idx]
  GOdata <- new("topGOdata",ontology = "BP",
                allGenes = genes_all,
                geneSel=function(p) p == 1,
                description ="inNetwork",
                annot=annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,
                     orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
  return(allRes)
  }else{
    return("No GO Analysis")
  }
}



