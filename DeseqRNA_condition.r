metadata <- read.delim("metadata.csv", header = TRUE, sep = ',', row.names=1)

countsData <- read.csv("data_t.csv", header = TRUE, row.names = 1, sep = ',', check.names = FALSE)

countsData <- t(countsData)
head(countsData)
#countsData <- countsData[, -1]

ncol(countsData)

nrow(metadata)

metadata$condition
unique(metadata$condition)

colData <- DataFrame(condition=factor(metadata$condition))

dds <- DESeqDataSetFromMatrix(ceiling(countsData), colData = colData, design=formula(~condition)) #gets stuck here

#contrast = c("strain", "C57BL/6T", "BALB/cT")
#contrast = c("libprep", "ribodepleted", "polyA")
contrast = c("condition", "FLT", "GC")

dds <- DESeq(dds)
res <- results(dds, contrast = contrast)

resOrdered <- res[order(res$padj), ]
resOrdered

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

