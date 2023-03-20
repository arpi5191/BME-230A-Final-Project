library(DESeq2)
metadata <- read.delim("metadata.csv", header = TRUE, sep = ',', row.names=1)

countsData <- read.csv("data_t.csv", header = TRUE, row.names = 1, sep = ',', check.names = FALSE)

countsData <- t(countsData)
head(countsData)

ncol(countsData)

nrow(metadata)

metadata$condition
unique(metadata$condition)

colData <- DataFrame(condition=factor(metadata$condition))

dds <- DESeqDataSetFromMatrix(ceiling(countsData), colData = colData, design=formula(~condition)) #gets stuck here

contrast = c("condition", "FLT", "GC")

dds <- DESeq(dds)
res <- results(dds, contrast = contrast)

resOrdered <- res[order(res$padj), ]
resOrdered

summary(resOrdered)

#out of 33067 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 169, 0.51%
#LFC < 0 (down)     : 81, 0.24%
#outliers [1]       : 0, 0%
#low counts [2]     : 10911, 33%
#(mean count < 1)

pVal <- 0.05
l2fc <- 2

significant <- resOrdered[!is.na(resOrdered$padj) & !is.nan(resOrdered$padj) & resOrdered$padj<pVal & abs(resOrdered$log2FoldChange)>=l2fc,]


selected <- rownames(significant)

selected

library(RSQLite)
packageVersion("RSQLite")
library("org.Mm.eg.db")
library("annotate")

u <- rownames(resOrdered)
u

geneIDs <- mapIds(org.Mm.eg.db, keys=row.names(significant), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")

geneIDs

#Gene Results

#condition
#ENSMUSG00000073842 ENSMUSG00000051747 ENSMUSG00000031640
#            "Mup7"              "Ttn"                 NA
#ENSMUSG00000055775 ENSMUSG00000024029 ENSMUSG00000038670
#            "Myh8"             "Tff3"           "Mybpc2"
#ENSMUSG00000030324 ENSMUSG00000058975
#             "Rho"            "Kcnc1"

# Effects of transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
library("RColorBrewer")

plotPCA(vsd, intgroup=c("condition"))

select <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:25]

pheatmap(assay(vsd)[select,],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize = 8,
         show_colnames = TRUE,
         main = "Heatmap of VST-transformed data - condition")


