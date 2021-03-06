## Get microarray data and find differentially expressed genes; find overlap between microarray DE genes and RNASeq DE genes.

#######################
## First run differential expression analysis on gse24206
## results from this section available as:
## load("data/GSE24206.RData")
## This R Script was adapted from one generated by geo2R: a gui on ncbi that does 
##  differential expression calculations on any GSE dataset; you can also pull the R code that it uses to call them. go here:
#   http://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE24206
# 
# + assign the samples to groups
# + click "top 250" at bottom of page
# + click the R script tab to see R code used to generate results

# note: 24206 is already rma-normalized and log2 transformed:
#   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM595407

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)


# load series and platform data from GEO (takes a while)
gset24206 <- getGEO("GSE24206", GSEMatrix =TRUE)
if (length(gset24206) > 1) idx <- grep("GPL570", attr(gset24206, "names")) else idx <- 1
gset24206 <- gset24206[[idx]]
rm(idx)
# make proper column names to match toptable 
fvarLabels(gset24206) <- make.names(fvarLabels(gset24206))

# group names for all samples
sml24206 <- c("G0","G0","G0","G0","G0","G0","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1")

# takes a gset and sml, returns ex = expression matrix, tT = results of DE analysis
# log2 transform -- does not go through on this data set; it has already been log transformed
ex24206 <- exprs(gset24206)

qx <- as.numeric(quantile(ex24206, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex24206[which(ex24206 <= 0)] <- NaN
            exprs(gset24206) <- log2(ex24206) }
# set up the data and proceed with analysis
fl <- factor(sml24206, levels = c("G0", "G1"))
gset24206$description <- fl
design <- model.matrix(~ description + 0, gset24206)
colnames(design) <- levels(fl)
fit <- lmFit(gset24206, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# tabulate the results and adjust pvalues
tT24206 <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex24206))

# load NCBI platform annotation
gpl <- annotation(gset24206)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd24206 <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT24206 <- tT24206[setdiff(colnames(tT24206), setdiff(fvarLabels(gset24206), "ID"))]
tT24206 <- merge(tT24206, ncbifd24206, by="ID")
tT24206 <- tT24206[order(tT24206$P.Value), ]  # restore correct order

tT24206 <- subset(tT24206, select=c("ID","adj.P.Val","P.Value","logFC","Gene.symbol","Gene.title"))
colnames(tT24206) <- c("probe", "padj", "pval", "log2FC", "gene", "gene_title")

# check that fold changes are going in the right direction: G1 is IPF
df <- data.frame(sample = colnames(ex24206), condition = fl, ZMAT3 = ex24206["219628_at",])
library(ggplot2)
library(grid)
ggplot(df, aes(x = condition, y = ZMAT3)) +
  geom_boxplot(outlier.size = 3) + 
  geom_point(aes(color = condition), size = 3, show_guide = FALSE, 
             position = position_jitter(width = 0.1, height = 0)) +
  theme_bw() + 
  theme(plot.title = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.5), vjust = 0.2),
        legend.title = element_text(size = rel(1.5)), 
        legend.text = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)), 
        legend.key.size = unit(1.5, "cm"),
        axis.ticks = element_blank(),
        panel.border = element_blank())
tT24206[1,]
# G1 = IPF is higher and log2FC is positive, so log2FC has proper sign

rm(cont.matrix, design, df, LogC, fit, fit2, fl, gpl, platf, qx, ncbifd24206)
save.image("data/GSE24206.RData")

#######################
## Run differential expression analysis on gse24206
## results from this section available as:
## load("data/GSE32537.RData")
## This R Script is adapted from one generated by geo2R. you can also pull the R code that it uses to call them. go here:
#   http://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE32537
# 
# + assign the samples to groups if it's not already done
# + click "top 250" at bottom of page
# + click the R script tab to see R code used to generate results


#   Differential expression analysis with limma
# takes a while

gset32537 <- getGEO("GSE32537", GSEMatrix =TRUE) 
if (length(gset32537) > 1) idx <- grep("GPL6244", attr(gset32537, "names")) else idx <- 1
gset32537 <- gset32537[[idx]]
rm(idx)
# make proper column names to match toptable 
fvarLabels(gset32537) <- make.names(fvarLabels(gset32537))

# group names for all samples
sml32537 <- c("G1","G1","X","G1","G1","G1","G1","X","G1","X","G1","X","G1","X","X","G1","G1","G1","G1","G1","X","G1","G1","G1","G1","G1","G1","G1","G1","X","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","X","X","X","G1","X","X","G1","G1","G1","G1","G1","G1","X","G1","G1","G1","G1","G1","G1","G1","G1","X","X","G1","G1","G1","X","X","X","G1","G1","G1","G1","G1","G1","G1","X","G1","X","X","G1","G1","G1","G1","X","X","G1","G1","G1","G1","G1","X","X","G1","G1","G1","G1","G1","G1","X","X","G1","G1","G1","G1","G1","G1","X","X","G1","G1","G1","G1","X","G1","X","G1","G1","G1","X","G1","X","X","G1","X","G1","G1","X","G1","G1","X","G1","G1","X","G1","G1","G1","X","G1","X","G1","X","G1","X","X","G1","X","G1","G1","X","G1","G1","G1","G1","G1","G1","X","G1","G1","G1","G1","G1","G1","X","G1","G1","G1","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0","G0");

# eliminate samples marked as "X"
sel32537 <- which(sml32537 != "X")
sml32537 <- sml32537[sel32537]
gset32537 <- gset32537[ ,sel32537]

ex32537 <- exprs(gset32537)
qx <- as.numeric(quantile(ex32537, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex32537[which(ex32537 <= 0)] <- NaN
            exprs(gset32537) <- log2(ex32537) }
# set up the data and proceed with analysis
fl <- factor(sml32537, levels = c("G0", "G1"))
gset32537$description <- fl
design <- model.matrix(~ description + 0, gset32537)
colnames(design) <- levels(fl)
fit <- lmFit(gset32537, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

summary(fit2)
# tabulate the results and adjust pvalues
tT32537 <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex32537))

# load NCBI platform annotation
gpl <- annotation(gset32537)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT32537 <- tT32537[setdiff(colnames(tT32537), setdiff(fvarLabels(gset32537), "ID"))]
tT32537 <- merge(tT32537, ncbifd, by="ID")
tT32537 <- tT32537[order(tT32537$P.Value), ]  # restore correct order

# select and rename cols
tT32537 <- subset(tT32537, select=c("ID","adj.P.Val","P.Value","logFC","Gene.symbol","Gene.title"))
colnames(tT32537) <- c("probe", "padj", "pval", "log2FC", "gene", "gene_title")

# check that fold changes are going in the right direction: G1 is IPF
subset(tT32537, gene == "COMP")
df <- data.frame(sample = colnames(ex32537), condition = fl, expression = ex32537["8035517",])
library(ggplot2)
library(grid)
ggplot(df, aes(x = condition, y = expression)) +
  geom_boxplot(outlier.size = 3) + 
  geom_point(aes(color = condition), size = 3, show_guide = FALSE, 
             position = position_jitter(width = 0.1, height = 0)) +
  theme_bw() + 
  theme(plot.title = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.5), vjust = 0.2),
        legend.title = element_text(size = rel(1.5)), 
        legend.text = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)), 
        legend.key.size = unit(1.5, "cm"),
        axis.ticks = element_blank(),
        panel.border = element_blank())
subset(tT32537, gene == "COMP")
# G1 = IPF is higher and log2FC is positive, so log2FC has proper sign

rm(cont.matrix, design, df, LogC, fit, fit2, fl, gpl, platf, qx, ncbifd, sel32537)
save.image("data/GSE32537.RData")

####################################
# compare to RNASeq DE genes from DESeq
setwd("~/Dropbox/ipf/ipf/runMicroarrayComparison/")
load("data/GSE24206.RData")
load("data/GSE32537.RData")

library(plyr)
source("../helperFunctions/myrw.R")

sigDESeq <- myread("../runDESeq/data/sigDESeq.txt")
sig24206 <- subset(tT24206, padj < 0.01)
sig32537 <- subset(tT32537, padj < 0.01)

# list the genes DE in all three experiments
overlap <- subset(sigDESeq, gene %in% sig24206$gene & gene %in% sig32537$gene)
overlap <- arrange(overlap, pval)
overlap <- subset(overlap, gene != "")

# get direction of FC for each experiment for each gene in overlaps
overlap$upreg_24206 <- vapply(overlap$gene, 
                              function(g){ (sig24206$log2FC[sig24206$gene == g] > 0)[1] }, # 24206 may have multiple probes in a gene; we take direction of the first
                              logical(1))

overlap$upreg_32537 <- vapply(overlap$gene, 
                              function(g){ (sig32537$log2FC[sig32537$gene == g] > 0)[1] },
                              logical(1))
overlap$upreg_deseq <- overlap$log2FoldChange > 0

# check that all directions are the same 
all(overlap$upreg_24206 == overlap$upreg_32537)
all(overlap$upreg_24206 == overlap$upreg_deseq)

# direction of fold change is the same for all probes in genes that come up in the overlap:
check_probes_agree <- function(sig){
  temp <- sig[,c("gene", "log2FC")]
  temp <- subset(temp, gene %in% overlap$gene)
  temp$upreg <- temp$log2FC > 0
  all_same <- function(values){all(values) | all(!values)}
  temp <- ddply(temp[,c("gene", "upreg")], ~ gene, summarize, same_dir = all_same(upreg))
  cat(all(temp$same_dir))
}
check_probes_agree(sig24206)
check_probes_agree(sig32537)


mywrite(overlap, "data/overlap.txt")
#######
