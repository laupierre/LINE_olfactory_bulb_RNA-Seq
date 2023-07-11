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

gp <- gp_up$result[ ,c("source", "term_id", "term_name", "p_value", "significant", "intersection_size", "intersection")]
colnames (gp)[colnames (gp) == "p_value"] <- "padj"
gp <- gp[order (gp$padj), ]
gp <- gp[gp$padj < 0.05, ]
head(gp)


# Adding back gene symbols

res <- list ()

for (i in (1:dim (gp)[1])) {
print (i)
entrez <- unlist (strsplit (as.vector (gp$intersection [i]), split=","))
symbols <- AnnotationDbi::select(org.Mm.eg.db, keys= entrez, keytype="ENSEMBL", columns="SYMBOL")
res[[i]] <- unique (symbols$SYMBOL)
}

res <- lapply (res, function (x) {paste (x,collapse="," )} )
x <- t (data.frame (res))
row.names (x) <- 1:dim (gp)[1]
gsedf <- cbind (data.frame (gp), genes= x)
gsedf <- gsedf[ ,!colnames (gsedf) %in% c("intersection")]
head (gsedf)

write.table (gsedf, "GOA_gprofiler_increased_genes.txt", sep="\t", quote=F, row.names=F)







