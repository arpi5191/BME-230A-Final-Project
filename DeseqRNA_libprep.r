metadata <- read.delim("metadata.csv", header = TRUE, sep = ',', row.names=1)

countsData <- read.csv("data_t.csv", header = TRUE, row.names = 1, sep = ',', check.names = FALSE)

countsData <- t(countsData)
head(countsData)
#countsData <- countsData[, -1]

ncol(countsData)

nrow(metadata)

metadata$libprep
unique(metadata$libprep)

colData <- DataFrame(libprep=factor(metadata$libprep))

dds <- DESeqDataSetFromMatrix(ceiling(countsData), colData = colData, design=formula(~libprep)) #gets stuck here

#contrast = c("strain", "C57BL/6T", "BALB/cT")
contrast = c("libprep", "ribodepleted", "polyA")
#contrast = c("condition", "FLT", "GC")

dds <- DESeq(dds)
res <- results(dds, contrast = contrast)

resOrdered <- res[order(res$padj), ]
resOrdered

pVal <- 0.05
l2fc <- 16

significant <- resOrdered[!is.na(resOrdered$padj) & !is.nan(resOrdered$padj) & resOrdered$padj<pVal & abs(resOrdered$log2FoldChange)>=l2fc,]


selected <- rownames(significant)

selected
#[1] "ENSMUSG00000104222" "ENSMUSG00000105361"
#[3] "ENSMUSG00000088025" "ENSMUSG00000026535"
#[5] "ENSMUSG00000068397" "ENSMUSG00000084383"
#[7] "ENSMUSG00000036322" "ENSMUSG00000069045"

library(RSQLite)
packageVersion("RSQLite")
library("org.Mm.eg.db")
library("annotate")

u <- rownames(resOrdered)
u

geneIDs <- mapIds(org.Mm.eg.db, keys=row.names(significant), column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")

geneIDs


#Gene Results

#libprep

#ENSMUSG00000104222 ENSMUSG00000105361 ENSMUSG00000088025
#NA                 NA            "Rprl3"
#ENSMUSG00000026535 ENSMUSG00000068397 ENSMUSG00000084383
#"Ifi202b"          "Gm10240"          "Gm13370"
#ENSMUSG00000036322 ENSMUSG00000069045
#"H2-Ea"            "Ddx3y"
