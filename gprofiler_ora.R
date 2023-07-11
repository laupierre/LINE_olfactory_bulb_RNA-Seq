library (gprofiler2)

deg <- readRDS ("deg.rds")


#### Positive genes
deg.pos <- na.omit (deg[deg$log2FoldChange >0 & deg$padj < 0.05, ])

genes <- as.character (deg.pos$Geneid)
genes <- gsub ("\\..*", "", genes)

geneUniverse <- as.character (deg$Geneid)
geneUniverse  <- gsub ("\\..*", "", geneUniverse)


# GO:BP vs GO
gp_up <- gost(query = genes, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO"),
             exclude_iea=TRUE, evcodes = TRUE, significant = FALSE, correction_method = "fdr", custom_bg= geneUniverse)

gp <- gp_up$result[ ,c("source", "term_id", "term_name", "p_value", "significant", "term_size", "intersection_size", "intersection")]
gp <- gp$result
gp <- gp[order (gp$p_value), ]
head(gp)



