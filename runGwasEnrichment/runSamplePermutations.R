####################################
## Second half: sample permutations 

## This script calculates a sample-permutation based empirical p-value for the enrichment of GWAS hits in DE genes.

## file structure:
## requires ../runDESeq/data/DESeq_012_norm7_fix_pval.RData and data/GWAS_genes.RData
## output to STDOUT

## Output from this script is given in data/perms.numOverlaps.100.txt
## The script was rerun with set.seed(1863), with output in data/perms.2.numOverlaps.100.txt, for a total of 200 permutations.



#!/usr/bin/env Rscript                                                                                                                                                                                                      

myread <- function(file, ...){
  return(read.table(file, header = TRUE, as.is = TRUE, sep = "\t", quote = "", comment = "", ...))
}
mywrite <- function(df, file, ...){
  write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, ...)
}

load("../runDESeq/data/DESeq_012_norm7_fix_pval.RData")
load("data/GWAS_genes.RData")

## renamed things for this script
res <- results.4.012
GWAS_genes <- ensembl_ids

library(DESeq)
set.seed(1862)

k <- integer(100)

for(i in 1:100){
  IPF.perm.design <- IPF.design
  IPF.perm.design$condition <- sample(IPF.perm.design$condition, 15, replace = FALSE)
  
  cds.perm <- newCountDataSet(genes, IPF.perm.design)
  cds.perm <- estimateSizeFactors(cds.perm)
  cds.perm <- estimateDispersions(cds.perm, method = "pooled-CR", modelFormula = form1)
  # filter out the lower 40th percentile genes (in terms of total expression across samples)                                                                                    
  cds.perm <- cds.perm[use,]
  fitFilt1  = fitNbinomGLMs(cds.perm, form1)
  fitFilt0  = fitNbinomGLMs(cds.perm, form0)
  pvals.perm <- pchisq(fitFilt0$deviance - fitFilt1$deviance,
                       attr(fitFilt0, "df.residual") - attr(fitFilt1, "df.residual"),
                       lower.tail = FALSE)
  padj.perm  <- p.adjust(pvals.perm, method="BH" )
  results.perm <- data.frame(gene_ID = gene_IDs[use], fit1_conv = fitFilt1$converged, fit0_conv = fitFilt0$converged,
                             pval = pvals.perm, padj = padj.perm, stringsAsFactors = FALSE)
  sig.perm <- subset(results.perm, padj < 0.05 & fit1_conv & fit0_conv)
  
  k[i] <- length(which(GWAS_genes$gene_ID %in% sig.perm$gene_ID))
  cat(paste0(i , "\t", k[i],"\n"), file = "perms.2.numOverlaps.100.txt", append = TRUE)
}

