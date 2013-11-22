## This first half of this analysis calculates the hypergeometric p-value associated with the overlap between DE genes and GWAS hits, 
##   and creates the QQ plots.
## The second half (runSamplePermutations.R) calculates a sample-permutation based empirical p-value for the enrichment.

## file structure: from working directory, data/
## requires ../runDESeq/data/DESeq_012_norm7_fix_pval.RData
## results from running this script are available by 
## load("data/GWAS_genes.RData")

####################
## First half: hypergeometric p-values
## Run this section twice, once for discovery SNPs and once for validation SNPs

library(biomaRt)
library(plyr)

source("../helperFunctions/getBiomart.R")

## read in the rsids of the 198 SNPs in the GWAS discovery set  
rsids <- read.table("data/discoverySNPs.txt", header = FALSE, as.is = TRUE, sep = "\t", quote = "", na = "", comment = "")
## for the validation set, use
# rsids <- read.table("data/validationSNPs.txt", header = FALSE, as.is = TRUE, sep = "\t", quote = "", na = "", comment = "")

## get genes nearest to each snp
snp_db <- useMart("snp", dataset="hsapiens_snp")
annotations <- getBM(c("refsnp_id", "allele", "minor_allele", "minor_allele_freq", "chr_name","chrom_start",                   
                       "chrom_strand","associated_gene",
                       "ensembl_gene_stable_id", "phenotype_description"),
                     filters="snp_filter",
                     values=rsids,
                     mart=snp_db)

annotations <- arrange(annotations, refsnp_id)
length(unique(annotations$refsnp_id))

ensembl_ids <- as.character(annotations$ensembl_gene_stable_id)
ensembl_ids <- unique(ensembl_ids)
ensembl_ids <- ensembl_ids[ensembl_ids != ""]

## hit genes contains genes with nearby SNPs in the GWAS discovery set
hit_genes <- data.frame(gene_ID = ensembl_ids, stringsAsFactors = FALSE)
hit_genes <- add_gene_names(hit_genes)

## load DE data
load("../runDESeq/data/DESeq_012_norm7_fix_pval.RData")

## total number of genes tested
nTested <- nrow(results.4.012)
## number of those that were significantly differentially expressed at FDR < 5%
nDEsig <- nrow(sig.4.012)
## number of DE genes in GWAS hits
nOverlap <- length(which(hit_genes$gene_ID %in% sig.4.012$gene_ID))
## number of GWAS genes tested for DE
nGwasTested <- length(which(hit_genes$gene_ID %in% results.4.012$gene_ID))
## calculate the hypergeometric p-value
phyper(q = nOverlap - 1, m = nDEsig, n = nTested - nDEsig, k = nGwasTested, lower.tail = FALSE)

save.image("data/GWAS_discovery.RData")
## for the validation set
# save.image("data/GWAS_validation.RData")

#####################
## create data for QQ plots: plot.df, plot.df.meta
## Run this section twice, once for discovery SNPs and once for validation SNPs

load("data/GWAS_discovery.RData")
## for validation set
# load("data/GWAS_validation.RData")

set.seed(1865)

num_sig <- numeric(10000)
perm_pval <- matrix(data = rep(-1, 10000*nGwasTested), nrow = 10000, ncol = nGwasTested)
perm_padj <- matrix(data = rep(-1, 10000*nGwasTested), nrow = 10000, ncol = nGwasTested)

for(i in 1:10000){
  perm_hit_genes <- results.4.012[sample(1:nrow(results.4.012), nGwasTested, replace = FALSE),]
  perm_pval[i,] <- arrange(perm_hit_genes, pval)$pval
  perm_padj[i,] <- arrange(perm_hit_genes, padj)$padj
  num_sig[i] <- length(which(perm_hit_genes$gene_ID %in% sig.4.012$gene_ID))
}
hist(num_sig, 100)

nlog_perm_pval <- -log10(perm_pval)
de.gwas <- subset(results.4.012, gene_ID %in% ensembl_ids)
plot.df <- data.frame(gwas.score = -log10(arrange(de.gwas, pval)$pval), mean.score = colMeans(nlog_perm_pval), sd = apply(nlog_perm_pval, 2, sd))

source("../helperFunctions/myrw.R")
mywrite(plot.df, "data/plot_df_discovery.txt")
## for validation set
# mywrite(plot.df, "data/plot_df_validation.txt")

#####################
## make QQ plots
library(plyr)
library(xtable)
library(ggplot2)
library(grid)
source("../helperFunctions/myrw.R")
plot.df <- myread("data/plot_df_discovery.txt")
plot.df.meta <- myread("data/plot_df_validation.txt")

g <- ggplot(plot.df, aes(x = gwas.score, y = mean.score)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = mean.score - sd, ymax = mean.score + sd, width = gwas.score/30.1), color="#E41A1C",  size = 1) +
  geom_point(size = 5, color = "#E41A1C") +
  xlab("-log10(pval) of GWAS genes") +
  ylab("-log10(pval) of random set") +
  coord_fixed() +
  xlim(0,8.7) + ylim(0,8.7) +
  theme_bw() + 
  theme(plot.title = element_text(size = rel(2)), 
        axis.title.x = element_text(size = rel(2), vjust = -0.2), axis.title.y = element_text(size = rel(2), vjust = 0.2),
        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)), 
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.key.size = unit(1.5, "cm"))
g + coord_flip()

g <- ggplot(plot.df.meta, aes(x = gwas.score, y = mean.score)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = mean.score - sd, ymax = mean.score + sd, width = gwas.score/30.1), color="#377EB8", size = 1) +
  geom_point(size = 5, color = "#377EB8") +
  xlab("-log10(pval) of GWAS genes") +
  ylab("-log10(pval) of random set") +
  coord_fixed() +
  xlim(0,8.7) + ylim(0,8.7) +
  theme_bw() + 
  theme(plot.title = element_text(size = rel(2)), 
        axis.title.x = element_text(size = rel(2), vjust = -0.2), axis.title.y = element_text(size = rel(2), vjust = 0.2),
        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)), 
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.key.size = unit(1.5, "cm"))
g + coord_flip()

temp1 <- transform(plot.df, set = "discovery")
temp2 <- transform(plot.df.meta, set = "validated")
plot.df.all <- rbind(temp1, temp2)

g <- ggplot(plot.df.all, aes(x = gwas.score, y = mean.score, color = set)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 5) +
  geom_line() + 
  xlab("-log10(pval) of GWAS genes") +
  ylab("mean -log10(pval) of 10,000 random sets") +
  coord_fixed() +
  xlim(0,7.5) + ylim(0,7.5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), guide = guide_legend(title = NULL)) +
  theme_bw() + 
  theme(plot.title = element_text(size = rel(2)), 
        axis.title.x = element_text(size = rel(2), vjust = -0.2), axis.title.y = element_text(size = rel(2), vjust = 0.2),
        legend.title = element_text(size = rel(2)), legend.text = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)), 
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.key.size = unit(1.5, "cm"))
g + coord_flip()
