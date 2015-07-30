add_biomart_gc <- function(df){
  # df has a column called "gene_ID" that contains ensembl ids
  # returns the same df with a column "gc_perc"
  check_gene_ID_col(df)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gc <- getBM(attributes=c("ensembl_gene_id","percentage_gc_content"),filters="ensembl_gene_id", values=df$gene_ID, mart=ensembl)
  gc <- gc[!duplicated(gc$ensembl_gene_id),]
  colnames(gc) <- c("gene_ID", "gc_perc")
  gene <- character(nrow(df))
  
  ## I use a for loop because I want to preserve the order of df. apply is no faster.
  for(i in 1:nrow(df)){
    temp <- gc$gc_perc[which(gc$gene_ID == df$gene_ID[i])]
    if(length(temp) == 0){
      temp <- NA
    }
    if(length(temp) > 1){
      temp <- paste(temp, collapse = ",")
      warning(paste0("multiple hits: ", temp))
    }
    gene[i] <- temp
  }
  df$gc_perc <- gene
  return(df)
}

add_biomart_gene <- function(df){
  # df has a column called "gene_ID" which is an ensembl ID
  # returns df with a new column : gene_biomart
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  annotations <- getGene(id = df$gene_ID, type = "ensembl_gene_id", mart = mart)
  annotations <- annotations[,c(1,2)]
  colnames(annotations) <- c("gene_ID","gene_biomart")
  gene <- character(nrow(df))
  for(i in 1:nrow(df)){
    temp <- annotations$gene_biomart[which(annotations$gene_ID == df$gene_ID[i])]
    if(length(temp) > 1){
      temp <- paste(temp, collapse = ",")
      warning(paste0("multiple hits: ", temp))
    }
    if(length(temp) == 0){
      gene[i] <- NA
    } else {
      gene[i] <- temp
    }
  }
  df$gene_biomart <- gene
  return(df)
}


## pull gene names from biomaRt
add_gene_names <- function(df) {
  # df has a column called "gene_ID" that contains ensembl ids
  # returns the same df with a column "gene" containing gene name
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  annotations <- getGene(id = df$gene_ID, type = "ensembl_gene_id", mart = mart)
  annotations <- annotations[,c(1,2)]
  colnames(annotations) <- c("gene_ID","gene")
  temp <- merge(df, annotations, by = "gene_ID", all.x = TRUE)
  return(temp)
}
