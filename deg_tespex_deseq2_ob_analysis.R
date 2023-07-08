library (DESeq2)
library (ggplot2)
library(ggrepel)
library (dplyr)


## See Github RNA-Seq_mouse/gene_annotation.R
#system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)



## normal STAR results from the RNA-Seq IIT pipeline
a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("Gene|bam", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))


counts <- annot <- a

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

row.names (counts) <- counts$Geneid
counts <- counts[ ,grep ("IIT", colnames (counts))]

colnames (counts) <- gsub ("IIT_SHA_LINE1_KO_", "", colnames (counts))
colnames (counts) <- gsub ("_S.*", "", colnames (counts))

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "condition", "area")
samples$area <- "OB"
samples$condition <- factor (c(rep ("KO",4), rep ("WT", 4)))
row.names (samples) <- colnames (counts)
samples$sample <- row.names (samples)



## TEspex result

tesp <- read.delim ("outfile.txt", row.names=1)
tesp$transname <- gsub ("#.*", "", row.names (tesp))
tesp$transname <- transname <- gsub ("_[3|5|o].*", "", tesp$transname)

transfamily <- gsub (".*/", "", row.names (tesp))
#tesp <- tesp[ ,grep ("IIT", colnames (tesp))]
annot_trans <- unique (data.frame (cbind (transname, transfamily)))

# aggregate multiple columns
tesp <- data.frame (tesp %>% group_by (transname) %>% summarise(across(everything(), sum)))

# same as in base: tesp <- aggregate(. ~ transname, tesp, sum)
# data.frame (t (tesp[tesp$transname == "L1MdA_II", ]))

colnames (tesp) <- gsub ("IIT_SHA_LINE1_KO_", "", colnames (tesp))
colnames (tesp) <- gsub ("_S[0-9].*", "", colnames (tesp))
row.names (tesp) <- tesp[ ,1]
tesp <- tesp[ ,-1]

tesp <- tesp[ ,colnames (tesp) %in% samples$sample]

idx <- match (colnames (tesp), samples$sample)
tesp <- tesp[ , idx]
stopifnot (colnames (tesp) == samples$sample)

counts <- rbind (counts, tesp)




## DESeq2 
dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = samples, design = ~ condition)
                                 
# keep <- rowSums(counts(dds)) >= 30
keep <- rowSums(counts(dds) >= 30) >= dim (counts)[2]/2
dds <- dds[keep,]
dds

dds$condition <- relevel(dds$condition, "WT")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=list("condition_KO_vs_WT"))

#res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x= TRUE)
colnames (res)[1] <- "Geneid"
res <- resa <- res[order (res$padj), ]

idx2 <- which (is.na (res$external_gene_name))
res$external_gene_name [idx2] <- res$Geneid[is.na (res$external_gene_name)]
res$gene_name [idx2] <- res$Geneid[idx2]
res$gene_type [idx2] <- "transposon"

for (i in (1:dim (res)[1])) {
if (res$gene_type[i] == "transposon") {
print (res$gene_name[i])
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

table (res$padj < 0.05)
# 153

boxplot (res$log2FoldChange)
abline (h=0)

write.table (res, "LINE_olfactory_bulb_differential_expression_with_transposons.txt", row.names=F, quote=F, sep="\t")



# PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed () + geom_label_repel (aes(label = name))

ggsave ("PCA plot.pdf", p1)




