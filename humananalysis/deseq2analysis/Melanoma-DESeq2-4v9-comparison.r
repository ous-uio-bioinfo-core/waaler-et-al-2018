### DESeq2 analysis of melanoma cell lines, comparing 4 versus 9 cell lines (Venn comparison distinct groups)

# Installing packages (if needed): 
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("pheatmap")


# Loading pachages: 
library("DESeq2")
library("pheatmap")


# Load and prepare dataset:
setwd("~/Desktop/Analysis")
data<-read.table("Counts4v9.txt",header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
rownames(data) <- make.unique(data[,1])
dataset <- data[,-1]
dataset <- round(dataset) # ajust counts to whole numbers - kallistooutput contains some decimal values, but DESeq2 needs whole counts.

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
group <- factor(c(rep("Venngroup9",18),rep("Venngroup4",8)),levels = c("Venngroup9","Venngroup4"))
cell.n <- factor(c(rep(1:9,each=2), rep(1:4,each=2)))
myfactors <- data.frame(group, treatment, cell.line, cell.n)
rownames(myfactors) <- colnames(dataset)
myfactors


# Preparing comparison matrix for DESeq2 analysis, according to the DESeq2 protocol:
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#levels-without-samples

m1 <- model.matrix(~ group + group:cell.n + group:treatment, myfactors)
all.zero <- apply(m1, 2, function(x) all(x==0))
all.zero
idx <- which(all.zero)
m1 <- m1[,-idx]
unname(m1)


# Generating a DESeq2 object:
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = dataset,
  colData = myfactors,
  design = m1)

dds <- ddsFullCountTable
as.data.frame(colData(dds)) #Checking factors in dds


# Running a DEseq2 analysis:
dds <- DESeq(dds)



# Inspecting results for difference in treatment effect between the two groups: 
resultsNames(dds)
res <- results(dds, contrast=list(c("groupVenngroup4.treatmentG007.LK"),c("groupVenngroup9.treatmentG007.LK")))
summary(res)
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), file=paste("4v9 gruppe-treatment-comparison.txt", sep="."), quote=FALSE, sep="\t")


# Generating heatmap for top hits:
rld <- rlog(dds, blind=FALSE)
res <- results(dds, contrast=list(c("groupVenngroup4.treatmentG007.LK"),c("groupVenngroup9.treatmentG007.LK")))
sig = which(res$padj < 0.01) 
mat <- assay(rld)[ sig, ]

c.mat <- cbind(mat[,2]-mat[,1],
               mat[,4]-mat[,3],
               mat[,6]-mat[,5],
               mat[,8]-mat[,7],
               mat[,10]-mat[,9],
               mat[,12]-mat[,11],
               mat[,14]-mat[,13],
               mat[,16]-mat[,15],
               mat[,18]-mat[,17],
               mat[,20]-mat[,19],
               mat[,22]-mat[,21],
               mat[,24]-mat[,23],
               mat[,26]-mat[,25])

colnames(c.mat) <- c("FemxI", "FemxV","MeWo","Skmel28","WM115","WM1341","WM1382","WM239","WM35","A375","WM1366","WM793","WM852")
venn <- c(rep("Venngroup9",9),rep("Venngroup4",4))
df <- as.data.frame(venn)
rownames(df) <- colnames(c.mat)
pheatmap(c.mat, annotation_col = df, cluster_cols = TRUE, cellheight=10, filename="Melanoma 4v9-padj0.01change-clust.pdf", fontsize_row = 4,main="4v9-DESeq2 analysis")
