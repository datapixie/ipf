#!/usr/bin/env Rscript

## chmod 755 this file to make it executable
## use as: ./runDEXSeq.R &> runDEXSeq.out

## requires exon.count files in folder data/, which are
##    named: reverse.<sample>.noM.mq10.namesorted.exon.count for 
##    <sample> = Norm1 - Norm6, Norm8, IPF1 - IPF8
## creates dexseq_testedDEU.RData, dexseq_results.txt

## These results are calculated on exon count data where: 
## 1. Norm7 has been removed
## 2. chrM has been filtered out
## 3. counts from reads with mapping quality < 10 have been filtered out
## 4. genes overlapping any other gene have been discarded

library(DEXSeq)
library(multicore)
library(DESeq)

print(Sys.time())

OUTPUT_FILENAME_BASE <- "data/dexseq"

## model formulae
formulaDisp <- count ~ sample + (condition + sex + dem.group) * exon
formula0 <- count ~ sample + (sex + dem.group) * exon + condition
formula1 <- count ~ sample + (sex + dem.group) * exon + condition * I(exon == exonID)

## exon count files: these have had chr M removed, and are filtered by mapping quality >= 10
files <- c("reverse.Norm1.noM.mq10.namesorted.exon.count","reverse.Norm2.noM.mq10.namesorted.exon.count","reverse.Norm3.noM.mq10.namesorted.exon.count",
           "reverse.Norm4.noM.mq10.namesorted.exon.count","reverse.Norm5.noM.mq10.namesorted.exon.count","reverse.Norm6.noM.mq10.namesorted.exon.count",
           "reverse.Norm8.noM.mq10.namesorted.exon.count",
           "reverse.IPF1.noM.mq10.namesorted.exon.count","reverse.IPF2.noM.mq10.namesorted.exon.count","reverse.IPF3.noM.mq10.namesorted.exon.count",
           "reverse.IPF4.noM.mq10.namesorted.exon.count","reverse.IPF5.noM.mq10.namesorted.exon.count","reverse.IPF6.noM.mq10.namesorted.exon.count",
           "reverse.IPF7.noM.mq10.namesorted.exon.count","reverse.IPF8.noM.mq10.namesorted.exon.count")
files <- vapply(files, function(x) paste0("data/", x), character(1))

samples <- c("Norm1","Norm2","Norm3",
             "Norm4","Norm5","Norm6",
             "Norm8",
             "IPF1","IPF2","IPF3",
             "IPF4","IPF5","IPF6",
             "IPF7","IPF8")
## covariates
sex <- as.factor(c("M", "M", "F", "F", "M", "F", "M",
                  "F", "F", "F", "M", "M", "M", "F", "M"))
dem.group <- as.factor(c(0,0,0,0,1,0,1,
                         2,1,0,0,0,0,2,0))
condition <- c(rep("Normal", 7), rep("IPF", 8))

IPF.design <- data.frame(row.names = samples, condition = condition, 
                                                  sex = sex, dem.group = dem.group)

## read counts
exon.counts <- read.HTSeqCounts(countfiles = files, design = IPF.design, flattenedfile = "/home/tnance/annotation/gencode.v14.annotation.dexseq.noOverlap.gff")
design(exon.counts)

############################
## Run DEXSeq

## estimate size factors
exon.counts <- estimateSizeFactors(exon.counts)
sizeFactors(exon.counts)

## estimate dispersions
cat("estimate dispersions...\n")
exon.counts <- estimateDispersions(exon.counts, formula = formulaDisp,
                                   nCores = 6, minCount = 10, maxExon = 70, quiet = TRUE, file = paste0(OUTPUT_FILENAME_BASE,"_estimateDispersions.out"))
cat("estimated dispersions...\n")
print(Sys.time())
exon.counts <- fitDispersionFunction(exon.counts)


## test for differential exon usage
exon.counts <- testForDEU(exon.counts, formula0 = formula0, formula1 = formula1, nCores = 6, quiet = TRUE, file = paste0(OUTPUT_FILENAME_BASE,"_testForDEU"))
exon.counts <- estimatelog2FoldChanges(exon.counts)
res <- DEUresultTable(exon.counts)
## save data here
save.image(file = paste0(OUTPUT_FILENAME_BASE, "_testedDEU.RData"))

print(Sys.time())

###########################
## helper functions to add gene symbols from biomart
library(biomaRt)

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

get_not_na <- function(x,y){
  if(!is.na(x)) return(x)
  return(y)
}


###########################
## add gene symbols to DEX results
res <- res
info <- exon.counts@featureData@data
## add exon location
res$chr_pos_str <- mapply(function(x,y,z) paste0(x, ":", y, "-", z), info$chr, info$start, info$end)
## add exon bin length
res$bin_len <- mapply(function(s,e) e - s + 1, info$start, info$end)
colnames(res)[1] <- "long_gene_ID"
res$long_gene_ID <- as.character(res$long_gene_ID)
res$exonID <- as.character(res$exonID)
res$gene_ID <- vapply(res$long_gene_ID, function(x)  unlist(strsplit(x, ".", fixed = TRUE))[1], character(1))
## add gene symbol
res <- add_biomart_gene(res)
## arrange by increasing pval
res <- arrange(res, pvalue)

write.table(res, "dexseq_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print(Sys.time())



