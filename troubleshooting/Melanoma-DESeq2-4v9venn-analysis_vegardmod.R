### DESeq2 analysis of melanoma cell lines, comparing 4 versus 9 cell lines (Venn comparison distinct groups)

# Installing packages (if needed): 
# source("http://bioconductor.org/biocLite.R")
# biocLite("ggplot2")
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
# install.packages("pheatmap")
# install.packages("gplots")
# source("https://bioconductor.org/biocLite.R")
# biocLite("regionReport")


# Loading pachages: 
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("pheatmap")
library("ggplot2")
library("genefilter")
library("DESeq2")
library("regionReport")

#### copy values from venngroup4 into 4 cellines from venngrou9 and remove the rest of venngroup9 to get a 4 vs 4 with identical values, and no diff between the two groups.
## in order to see what deseq reports as DEG for the different comparisons.
sanitycheckmode = FALSE

# Load and prepare dataset:
# setwd("~/Desktop/Analysis")
data<-read.table("Counts4v9.txt",header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
rownames(data) <- make.unique(data[,1])
dataset <- data[,-1]
dataset <- round(dataset) # The DESeq2 package only works on raw counts, witch should be whole numbers. Not all numbers in this counts dataset is whole numbers and thus the values in the dataset need to be rounded. 
dataset <- as.matrix(dataset)

# Setting factors/metadata
cell.line <- factor(c("FemxI","FemxI", 
                      "FemxV", "FemxV",
                      "MeWo","MeWo",
                      "Skmel28","Skmel28",
                      "WM115", "WM115",
                      "WM1341","WM1341",
                      "WM1382","WM1382",
                      "WM239","WM239", 
                      "WM35", "WM35",
                      "A375", "A375",
                      "WM1366", "WM1366",
                      "WM793", "WM793",
                      "WM852","WM852"), 
                    levels = c("FemxI", 
                               "FemxV",
                               "MeWo",
                               "Skmel28",
                               "WM115",
                               "WM1341",
                               "WM1382",
                               "WM239",
                               "WM35",
                               "A375",
                               "WM1366",
                               "WM793",
                               "WM852"))

treatment <- factor(c(rep(c("DMSO","G007-LK"),13)), 
                    levels = c("DMSO","G007-LK"))
venngroup <- factor(c(rep("Venngroup9",18),rep("Venngroup4",8)),levels = c("Venngroup9","Venngroup4"))
cell.n <- factor(c(rep(1:9,each=2), rep(1:4,each=2)))
# cell.n <- factor(c(rep(9:1,each=2), rep(1:4,each=2)))
myfactors <- data.frame(venngroup, treatment, cell.line, cell.n)
rownames(myfactors) <- colnames(dataset)
myfactors


### make a dummy analysis with the exact same measurements in a 4 vs 4 comparison

if(sanitycheckmode)
{
	myfactors = myfactors[c(1:8, 19:26),]
	dataset[, 1:8] = dataset[, 19:26] # copying data
	dataset = dataset[, c(1:8, 19:26)] # subsetting so it is examctly the same 
	print(dataset[1:4,])
}
# Prepping for a DESeq2 analysis: 
# To do a venngroup:treatment interaction analysis, a comparison matrix is needed due to #Levels without samples". 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#levels-without-samples

m1 <- model.matrix(~ venngroup + venngroup:cell.n + venngroup:treatment, myfactors)
all.zero <- apply(m1, 2, function(x) all(x==0))
all.zero
idx <- which(all.zero)
m1 <- m1[,-idx]
unname(m1)



# Generating a DESeq2 object:
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = dataset,
  colData = myfactors,
  design = ~ 1)




dds <- ddsFullCountTable
as.data.frame(colData(dds)) #Checking factors


# Running a DEseq2 analysis:
dds <- DESeq(dds, full=m1)


# # PCA of DESeq2data
# rld <- rlog(dds, blind=FALSE)
# pdf(file="PCA Melanoma DESeq2 4v9 comparison.pdf") 
# plotPCA(rld, intgroup=c("treatment"))
# plotPCA(rld, intgroup=c("venngroup"))
# pcaData <- plotPCA(rld, intgroup=c("cell.n", "treatment", "venngroup"), returnData=TRUE)
# ggplot(pcaData, aes(PC1, PC2, color=cell.n, shape=treatment)) + geom_point(size=3) + coord_fixed()
# dev.off()
# 
# 
# # Get a list of the top variable genes used to plot PCA
# mat <- assay(rld)
# rv = apply(mat, 1, var)
# ntop = 500
# select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# pca = prcomp(t(mat[select,]))
# pcaweights <- as.data.frame(pca$rotation)
# write.table(pcaweights, file=paste("4v9 gruppe-treatment-comparison-PCAweights.txt", sep="."), quote=FALSE, sep="\t")
# 

# Inspecting results for the primary condition: 
resultsNames(dds) # The interaction term, "venngroupVenngroup4.treatmentG007.LK" is the difference between the effect of the treatment across the groups, 
                  # (V4_treatment - V4_ctrl) - (V9_treatment - V9_ctrl) - and this is what comes up as primary condition
res <- results(dds)
summary(res)
hist( res$pvalue, breaks=100)

# Writing file of results
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), file=paste("4v9 gruppe-treatment-comparison.txt", sep="."), quote=FALSE, sep="\t")

### vegards alternative
res2 <- results(dds, list(c("venngroupVenngroup4.treatmentG007.LK"),c("venngroupVenngroup9.treatmentG007.LK")))
summary(res2)
hist( res2$pvalue, breaks=100)
res2Ordered <- res2[order(res2$padj),]
write.table(as.data.frame(res2Ordered), file=paste("4v9_gruppe-treatment-comparison_vegards.txt", sep="."), quote=FALSE, sep="\t")


# p-value-histogram
#pdf(file="p-value.pdf") 
#hist( res$pvalue, breaks=100)
#dev.off()


## p-value analysis... (as suggested by the DESeq2 protocol)
#use <- res$baseMean > metadata(res)$filterThreshold
#h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
#h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
#colori <- c(`do not pass`="khaki", `pass`="powderblue")
#pdf(file="pval-analysis.pdf")
#barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
#text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
#legend("topright", fill=rev(colori), legend=rev(names(colori)))
#dev.off()


## Generating heatmap for top hits: 
#res <- results(dds, lfcThreshold=0.1 )
#sig = which(res$padj < 0.1) 
#mat <- assay(rld)[ sig, ]
#mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(rld)[,c("treatment","venngroup")])
#pheatmap(mat, annotation_col=df[2], cluster_cols = FALSE, cellheight=10, filename="Melanoma 4v9complex-lfc0.1-padj0.1.pdf", fontsize_row = 4)
#pheatmap(mat, annotation_col=df[2], cluster_cols = TRUE, cellheight=10, filename="Melanoma 4v9complex-lfc0.1-padj0.1clust.pdf", fontsize_row = 4)


# REGIONREPORT
report <- DESeq2Report(dds, 'DESeq2-4v9-analysis', c('venngroup', 'treatment'),
                       outdir = 'Analyse')



library(tidyverse)

sa = myfactors %>% mutate(sample=rownames(myfactors))

## must normalize before recalculating FC and plotting. DEseq probably does this internally.
cs = colSums(dataset)
normdataset = t(t(dataset)/cs) * 1000000 
colSums(normdataset)


#focusgenes = c("TRAM2")
dlong = as.data.frame(normdataset) %>%
	mutate(symbol = rownames(normdataset)) %>%
#	dplyr::filter(symbol %in% focusgenes) %>%
	gather(sample, count, -symbol) %>%
	left_join(sa)


# linechart by gene and cell line
# focusgenes = c("TRAM2", "CAV1", "BCL2L1")
# p = dlong %>% 
#	mutate(venngroup = factor(venngroup, levels=c("Venngroup4", "Venngroup9"))) %>%
#	dplyr::filter(symbol %in% focusgenes) %>%
#	ggplot( aes(x=treatment, y=count, group=cell.line)) +
#	geom_line( aes(color=cell.line), size=1) +
#	facet_wrap(~ symbol + venngroup, ncol=2, scales="free_y")
#print(p)
# ggsave("linechart_some_genes.pdf", p, width=15, height=30, units = "cm")


martinssammenlikning = as.data.frame(resOrdered) %>%
	mutate(symbol = rownames(resOrdered)) %>%
	dplyr::select( symbol, log2FoldChange) %>%
	rename( martinslog2FC = log2FoldChange)

vegardssammenlikning = as.data.frame(res2Ordered) %>%
	mutate(symbol = rownames(res2Ordered)) %>%
	dplyr::select( symbol, log2FoldChange) %>%
	rename( vegardslog2FC = log2FoldChange)

ngenes = 100
comptable = dlong %>%
	dplyr::filter(symbol %in% c( martinssammenlikning$symbol[1:ngenes],vegardssammenlikning$symbol[1:ngenes] )) %>%
	mutate(log2count = log2(count)) %>%
	dplyr::select(-sample, -count, -cell.n) %>%
	spread(treatment, log2count) %>%
	rename(G007LK = "G007-LK") %>%
	mutate(log2change = G007LK- DMSO) %>%
	group_by(venngroup, symbol) %>%
	summarize(meanchange = mean((log2change))) %>%
	ungroup() %>%
	spread(venngroup, meanchange) %>%
	mutate(betweengroupdiff = Venngroup4 - Venngroup9) %>%
	left_join(martinssammenlikning) %>%
	left_join(vegardssammenlikning)


p = comptable %>%
	gather(venngroup, recalculated, Venngroup9, Venngroup4, betweengroupdiff ) %>%
	gather(method, deseqlog2FC, martinslog2FC, vegardslog2FC) %>%
	ggplot(aes(x=deseqlog2FC, y=recalculated)) +
	geom_point(size=0.5) +
	facet_wrap( ~venngroup + method , ncol=2)

ggsave("deseqlog2FC_vs_recalculated.pdf", p, width=15, height=20, units = "cm")



### check overlap with limmaresults

limmares = read_tsv("../humananalysis/output_sleuth_human_2018-11-05/DEGtest_results_2018-11-05/venngroupcentral_vs_venngroupgreen_change_limma_results_2018-11-05.txt")
limmares = limmares %>%
	arrange(P.Value)
	
noverlapckeck = 100
sum(limmares$symbol[1:noverlapckeck] %in% rownames(resOrdered)[1:noverlapckeck] )
sum(limmares$symbol[1:noverlapckeck] %in% rownames(res2Ordered)[1:noverlapckeck] )

# difference in filtered out genes
sum(limmares$symbol %in%  rownames(res2Ordered) )
sum(rownames(res2Ordered) %in% limmares$symbol)


# Half of limma top genes not even present in the input for deseq2, i.e filtered out
sum(rownames(res2Ordered)[1:noverlapckeck] %in% limmares$symbol)
sum(limmares$symbol[1:noverlapckeck] %in%  rownames(res2Ordered) )


# limma has a left sided p-value peak also.
hist(limmares$P.Value, breaks=100)



