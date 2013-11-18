#!/usr/bin/env Rscript
## chmod 755 this file to make it executable
## use as: ./runDESeq.R &> runDESeq.out

## file structure: from working directory, plots/, data/
## requires data/genes.txt
## creates DESeq_012_norm7_fix_pval.RData in working directory
## creates plots/deseq_012.pdf and plots/deseq_000.pdf

## This script runs DESeq to find differentially expressed genes from
## the prepared gene count file genes.txt (which still containes Norm7).
## We calculate DESeq results
## with and without covariates, each with and without a filter on 
## overall expression level. 

## Figures are output to the folder plots/ as pdf files. 
## Figures for calculations without covariates are in file deseq_000.pdf
## and figures with covariates are in deseq_012.pdf

## notation: "012" stands for things with covariates included; 
## "4" stands for things were I filtered out the lower 40% expressed 
## reads (by sums across samples)


#########
# know when we start
print(Sys.time())

# load required libraries
library(plyr)
library(DESeq)
library(biomaRt)
library(ggplot2)
library(gplots)

###################################
## put all covariates into IPF.design
samples <- c("Norm1", "Norm2", "Norm3", "Norm4", "Norm5", "Norm6", "Norm8", "IPF1", "IPF2", "IPF3", "IPF4", "IPF5", "IPF6", "IPF7", "IPF8")
conditions <- factor(c(rep("Norm", 7), rep("IPF", 8)))
sex <- c("M", "M", "F", "F", "M", "F", "M", 
         "F", "F", "F", "M", "M", "M", "F", "M")

# demography groups: 0 = European ancestry, 2 = Asian ancestry, 1 = admixed or American ancestry
dem.group <- as.factor(c(0,0,0,0,1,0,1,
                         2,1,0,0,0,0,2,0))

for(i in 1:15){
  if(dem.group[i] == 0){
    cat("European ancestry\n")
  } else if (dem.group[i] == 1){
    cat("Admixed or American ancestry\n")
  } else {
    cat("Asian ancestry\n")
  }
}

IPF.design <- data.frame(row.names = samples, condition = conditions, 
                         sex = sex, dem.group = dem.group)

###################################
## read the gene counts
genes <- read.table("data/genes.txt", header = TRUE, as.is = TRUE, sep = "\t",
                    quote = "", na = "NA", comment = "", row.names = 1)
## remove Norm 7 which is an expression outlier
genes <- genes[,-7]

###################################
## helper functions

## pull gene names from biomaRt
add.gene.names <- function(df) {
  # df has a column called "gene_ID" that contains ensembl ids
  # returns the same df with a column "gene" containing gene name
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  annotations <- getGene(id = df$gene_ID, type = "ensembl_gene_id", mart = mart)
  annotations <- annotations[,c(1,2)]
  colnames(annotations) <- c("gene_ID","gene")
  temp <- merge(df, annotations, by = "gene_ID", all.x = TRUE)
  return(temp)
}

## draw MA plots
plotDE <- function(res){
  plot(
    res$baseMean,
    res$log2FoldChange,
    log="x", pch=20, cex=.3,
    col = ifelse( res$padj < .05, "red", "black" ),
    xlab = "Mean of counts / size factors",
    ylab = "log 2 Fold Change (IPF/normal)",
    main = "STAR DESeq results (all lung samples)",
    #   xlim = c(1e-1, 1e6),
    #   ylim = c(-10, 10)
  )
  points(res$baseMean[res$padj < .05], res$log2FoldChange[res$padj < .05],
         pch = 20, cex = 1,
         col = "red")
}

###################################
## run DESeq with no covariates
cds.null <- newCountDataSet(genes, conditions)
cds.null <- estimateSizeFactors(cds.null)
cds.null <- estimateDispersions(cds.null)
results.null <- nbinomTest(cds.null, "Norm", "IPF")

## get sig genes from model without covariates at 5% FDR,
## get their gene names
sig.null <- subset(results.null, padj < .05)
sig.null <- arrange(sig.null, pval)
colnames(sig.null)[1] <- "gene_ID"
sig.null <- add.gene.names(sig.null)

## Plot dispersion estimates and make MA plots for model without 
## covariates
pdf(file = "plots/deseq_000.pdf", paper="USr")
plotDispEsts(cds.null)
plotDE(results.null)
plotMA(results.null)
par(mfrow=c(2,2))
hist(results.null$pval, breaks = 500)
hist(results.null$padj, breaks = 500)
hist(results.null$padj[results.null$padj<.9], breaks = 500)
par(mfrow = c(1,1))

###################################
## filter out 40% lowest-expressed genes (by sum of counts across samples)

rs <- rowSums(counts(cds.null))
# plot showing that lowly expressed genes do not have significant p-values
plot(rank(rs)/length(rs), -log10(results.null$pval), pch=16, cex=0.45,
     xlab = "total count across all samples (percentile)", ylab = "-log10(pval)",
     main = "pvals vs total read count")

theta <- 0.4
use <- (rs > quantile(rs, probs=theta))
# show how many genes are filtered
table(use)
textplot(table(use))

## remove the filtered genes and recalculate differential expression
cds.4.null <- cds.null[use,]
results.4.null <- nbinomTest(cds.4.null, "Norm", "IPF")
sig.4.null <- subset(results.4.null, padj < .05)
sig.4.null <- arrange(sig.4.null, pval)
colnames(sig.4.null)[1] <- "gene_ID"
sig.4.null <- add.gene.names(sig.4.null)

## MA plot for filtered results
plotDE(results.4.null)
plotMA(results.4.null)

## histogram showing that filtered genes occur evenly across p-values
hist(results.4.null$pval, breaks = 500)
hist(results.4.null$padj, breaks = 500)
hist(results.4.null$padj[results.4.null$padj<.9], breaks = 500)
h1 <- hist(results.null$pval[!use], breaks=50, plot=FALSE)
h2 <- hist(results.null$pval[use], breaks=50, plot=FALSE)
colori <- c( `do not pass` ="khaki",  `pass filter` ="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori,
        space = 0, main = "", ylab="frequency", xlab = "pval")
legend("topleft", fill=rev(colori), legend=rev(names(colori)), cex = .8)
## finished with plots for no covariates
garbage <- dev.off()

###################################
## Do DESeq analysis with covariates: demographic group and sex

# model formulae
form0 <- count ~ dem.group + sex
form1 <- count ~ condition + dem.group + sex
cds <- newCountDataSet(genes, IPF.design)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method = "pooled-CR", modelFormula = form1)

# visualize the dispersion estimates
pdf(file = paste0("plots/deseq_012.pdf"), paper = "USr")
plotDispEsts(cds)

fit1 <- fitNbinomGLMs(cds, form1)
fit0 <- fitNbinomGLMs(cds, form0)
#    pval = nbinomGLMTest(fit1, fit0)

## Avoid getting p-values of zero
pval <- pchisq(fit0$deviance - fit1$deviance,
               attr(fit0, "df.residual") - attr(fit1, "df.residual"),
               lower.tail = FALSE)
padj <- p.adjust(pval, method="BH")
gene_IDs <- rownames(genes)
results.012 <- data.frame(gene_ID = gene_IDs, fit1_conv = fit1$converged, fit0_conv = fit0$converged,
                          pval = pval, padj = padj, stringsAsFactors = FALSE)
results.012$baseMean <- results.null$baseMean
results.012$log2FoldChange <- results.null$log2FoldChange

sig.012 <- subset(results.012, padj < .05)
colnames(sig.012)[1] <- "gene_ID"
sig.012 <- add.gene.names(sig.012)
sig.012 <- arrange(sig.012, padj)

## visualize results
plotDE(results.012)
plotMA(results.012)
par(mfrow=c(2,2))
hist(results.012$pval, breaks = 500)
hist(results.012$padj, breaks = 500)
hist(results.012$padj[results.012$padj<.9], breaks = 500)
par(mfrow=c(1,1))


###################################
## filter out 40% lowest-expressed genes (by sum of counts across samples) 
## and use covariates

rs <- rowSums(counts(cds))
plot(rank(rs)/length(rs), -log10(results.012$pval), pch=16, cex=0.45,
     xlab = "total count across all samples (percentile)", ylab = "-log10(pval)",
     main = "pvals vs total read count")
theta <- 0.4
use <- (rs > quantile(rs, probs=theta))

cds.4.012 <- cds[use,]
fitFilt1  = fitNbinomGLMs(cds.4.012, form1)
fitFilt0  = fitNbinomGLMs(cds.4.012, form0)
#   pvalsFilt = nbinomGLMTest(fitFilt1, fitFilt0)
pvalsFilt <- pchisq(fitFilt0$deviance - fitFilt1$deviance,
                    attr(fitFilt0, "df.residual") - attr(fitFilt1, "df.residual"),
                    lower.tail = FALSE)
padjFilt  = p.adjust(pvalsFilt, method="BH" )
results.4.012 <- data.frame(gene_ID = gene_IDs[use], fit1_conv = fitFilt1$converged, fit0_conv = fitFilt0$converged,
                            pval = pvalsFilt, padj = padjFilt, stringsAsFactors = FALSE)

results.4.012$baseMean <- results.4.null$baseMean
results.4.012$log2FoldChange <- results.4.null$log2FoldChange

## visualize results
plotDE(results.4.012)
plotMA(results.4.012)
par(mfrow=c(2,2))
hist(results.4.012$pval, breaks = 500)
hist(results.4.012$padj, breaks = 500)
hist(results.4.012$padj[results.4.012$padj<.9], breaks = 500)
par(mfrow=c(1,1))

## table showing numbers with filter
padjFiltForComparison = rep(+Inf, length(results.012$padj))
padjFiltForComparison[use] = results.4.012$padj
tab3 <- table( "no filtering" = results.012$padj < .05,
               "with filtering"  = padjFiltForComparison < .05 )
temp <- apply(addmargins(tab3), c(1,2), as.integer)
textplot(temp)

## histogram showing that filtered genes occur evenly across p-values
h1 <- hist(results.012$pval[!use], breaks=50, plot=FALSE)
h2 <- hist(results.012$pval[use], breaks=50, plot=FALSE)
colori <- c( `do not pass` ="khaki",  `pass filter` ="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori,
        space = 0, main = "", ylab="frequency", xlab = "pval")
legend("topleft", fill=rev(colori), legend=rev(names(colori)), cex = .8)
garbage <- dev.off()

## get significant results, add gene names
sig.4.012 <- subset(results.4.012, padj < .05)
colnames(sig.4.012)[1] <- "gene_ID"
sig.4.012 <- add.gene.names(sig.4.012)
sig.4.012 <- arrange(sig.4.012, padj)

sig.4.012 <- subset(sig.4.012, fit0_conv == TRUE & fit1_conv == TRUE)
sig.012 <- subset(sig.012, fit0_conv == TRUE & fit1_conv == TRUE)

rm(garbage)
save.image("data/DESeq_012_norm7_fix_pval.RData")

## write results sig at 1% FDR for later
source("../helperFunctions/myrw.R")
mywrite(subset(sig.4.012, padj < .01), "data/sigDESeq.txt")

print(Sys.time())
