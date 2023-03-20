metadata <- read.delim("metadata.csv", header = TRUE, sep = ',', row.names=1)

countsData <- read.csv("data_t.csv", header = TRUE, row.names = 1, sep = ',', check.names = FALSE)

countsData <- t(countsData)
head(countsData)
#countsData <- countsData[, -1]

ncol(countsData)

nrow(metadata)

metadata$strain
unique(metadata$strain)

colData <- DataFrame(strain=factor(metadata$strain))

dds <- DESeqDataSetFromMatrix(ceiling(countsData), colData = colData, design=formula(~strain)) #gets stuck here

#"C57BL6T" "C57BL6J" "BALBcT"

contrast = c("strain", "C57BL6T", "BALBcT")
#contrast = c("strain", "C57BL6T", "C57BL6J")
#contrast = c("strain", "BALBcT", "C57BL6J")

dds <- DESeq(dds)
res <- results(dds, contrast = contrast)

resOrdered <- res[order(res$padj), ]
resOrdered

pVal <- 0.05
l2fc <- 16

significant <- resOrdered[!is.na(resOrdered$padj) & !is.nan(resOrdered$padj) & resOrdered$padj<pVal & abs(resOrdered$log2FoldChange)>=l2fc,]
significant


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

#strain - "C57BL6T", "BALBcT"

#ENSMUSG00000105796 ENSMUSG00000089803 ENSMUSG00000020912
#                NA                 NA            "Krt12"
#ENSMUSG00000042240 ENSMUSG00000006546 ENSMUSG00000093484
#          "Crybb2"           "Cryba2"                 NA
#ENSMUSG00000024041 ENSMUSG00000068165 ENSMUSG00000098986
#           "Cryaa"                 NA                 NA
#ENSMUSG00000073658 ENSMUSG00000086691 ENSMUSG00000068457
#           "Crygb"                 NA              "Uty"
#ENSMUSG00000025952 ENSMUSG00000083405 ENSMUSG00000105704
#        "Crygc"                 NA                 NA
#ENSMUSG00000080859
#       "Rpl10-ps1"

#strain  - "C57BL6T", "C57BL6J
#ENSMUSG00000068457 ENSMUSG00000069045 ENSMUSG00000069049 ENSMUSG00000020912 ENSMUSG00000006546
#             "Uty"            "Ddx3y"          "Eif2s3y"            "Krt12"           "Cryba2"
#ENSMUSG00000042240 ENSMUSG00000073658
#          "Crybb2"            "Crygb"

#strain - "BALBcT", "C57BL6J"

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
library("RColorBrewer")

plotPCA(vsd, intgroup=c("strain"))

select <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:25]

pheatmap(assay(vsd)[select,],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize = 8,
         show_colnames = TRUE,
         main = "Heatmap of VST-transformed data - strain")
