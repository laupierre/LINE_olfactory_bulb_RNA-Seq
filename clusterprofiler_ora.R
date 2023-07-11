#library (org.Hs.eg.db)
library (org.Mm.eg.db)
library(clusterProfiler)

deg <- readRDS ("deg.rds")


#### Positive genes
deg.pos <- na.omit (deg[deg$log2FoldChange >0 & deg$padj < 0.05, ])

genes <- as.character (deg.pos$Geneid)
genes <- gsub ("\\..*", "", genes)

## both DE gene list and gene universe (with ENSEMBL) must be annotated in Entrez IDs (EG)

deGenes <- unlist(mget(genes, envir=org.Mm.egENSEMBL2EG, ifnotfound = NA))
head (deGenes)
#ENSMUSG00000034312 ENSMUSG00000025582 ENSMUSG00000028518 ENSMUSG00000021943 
#          "232227"            "18164"           "108079"            "14560" 

geneUniverse <- as.character (deg$Geneid)
geneUniverse  <- gsub ("\\..*", "", geneUniverse)
geneUniverse <- unlist(mget(geneUniverse, envir=org.Mm.egENSEMBL2EG, ifnotfound = NA))



# over-representation analysis of GO pathways (ALL or BP)

ora.go <- enrichGO (gene = deGenes, ont = "ALL",
                    OrgDb ="org.Mm.eg.db",
                    universe = geneUniverse,
                    readable=TRUE,
                    pvalueCutoff = 0.05, pAdjustMethod = "BH")

dim (ora.go)
# 1

tab.go <- as.data.frame(ora.go)
tab.go$geneID <- gsub ("/", ",", tab.go$geneID)
head (tab.go)

write.table (tab.go, "GOA_increased_genes.txt", sep="\t", quote=F, row.names=F)



# over-representation of KEGG pathways

ora.kegg <- enrichKEGG(gene = deGenes,
                       organism = 'mmu',
                       universe = geneUniverse,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")

dim (ora.kegg)
# 0

tab.kegg <- as.data.frame(ora.kegg)
head (tab.kegg) 


# Adding back gene symbols

res <- list ()

for (i in (1:dim (tab.kegg)[1])) {
print (i)
entrez <- unlist (strsplit (as.vector (tab.kegg$geneID [i]), split="/"))
symbols <- AnnotationDbi::select(org.Mm.eg.db, keys= entrez, keytype="ENTREZID", columns="SYMBOL")
res[[i]] <- unique (symbols$SYMBOL)
}

res <- lapply (res, function (x) {paste (x,collapse="," )} )
x <- t (data.frame (res))
row.names (x) <- 1:dim (tab.kegg)[1]
gsedf <- cbind (data.frame (tab.kegg), genes= x)
gsedf <- gsedf[ ,!colnames (gsedf) %in% c("geneID")]
head (gsedf)



#### Negative genes

deg.neg <- na.omit (deg[deg$log2FoldChange <0 & deg$padj < 0.05, ])
genes <- as.character (deg.neg$Geneid)
genes <- gsub ("\\..*", "", genes)

## both DE gene list and gene universe (with ENSEMBL) must be annotated in Entrez IDs (EG)

deGenes <- unlist(mget(genes, envir=org.Mm.egENSEMBL2EG, ifnotfound = NA))
head (deGenes)
#ENSMUSG00000034312 ENSMUSG00000025582 ENSMUSG00000028518 ENSMUSG00000021943 
#          "232227"            "18164"           "108079"            "14560" 

geneUniverse <- as.character (deg$Geneid)
geneUniverse  <- gsub ("\\..*", "", geneUniverse)
geneUniverse <- unlist(mget(geneUniverse, envir=org.Mm.egENSEMBL2EG, ifnotfound = NA))



# over-representation analysis of GO pathways

ora.go <- enrichGO (gene = deGenes, ont = "ALL",
                    OrgDb ="org.Mm.eg.db",
                    universe = geneUniverse,
                    readable=TRUE,
                    pvalueCutoff = 0.05, pAdjustMethod = "BH")

dim (ora.go)
# 29

tab.go <- as.data.frame(ora.go)
tab.go$geneID <- gsub ("/", ",", tab.go$geneID)
head (tab.go)

write.table (tab.go, "GOA_decreased_genes.txt", sep="\t", quote=F, row.names=F)



# over-representation of KEGG pathways

ora.kegg <- enrichKEGG(gene = deGenes,
                       organism = 'mmu',
                       universe = geneUniverse,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")

dim (ora.kegg)
# 1

tab.kegg <- as.data.frame(ora.kegg)
head (tab.kegg) 


# Adding back gene symbols

res <- list ()

for (i in (1:dim (tab.kegg)[1])) {
print (i)
entrez <- unlist (strsplit (as.vector (tab.kegg$geneID [i]), split="/"))
symbols <- AnnotationDbi::select(org.Mm.eg.db, keys= entrez, keytype="ENTREZID", columns="SYMBOL")
res[[i]] <- unique (symbols$SYMBOL)
}


res <- lapply (res, function (x) {paste (x,collapse="," )} )
x <- t (data.frame (res))
row.names (x) <- 1:dim (tab.kegg)[1]
gsedf <- cbind (data.frame (tab.kegg), genes= x)
gsedf <- gsedf[ ,!colnames (gsedf) %in% c("geneID")]
head (gsedf)


