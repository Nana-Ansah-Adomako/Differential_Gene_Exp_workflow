library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(DOSE)
library(rio)

#===================importing counts_data file 

counts_data <- import("file_path")
head(counts_data)

#==================subsetting only the required columns
counts_data <- subset(counts_data,
                  select=c(GeneID, 
                           GSM1275862, GSM1275863, 
                           GSM1275866, GSM1275867,
                           GSM1275870, GSM1275871, 
                           GSM1275874, GSM1275875 ))

colnames(counts_data)
rownames(counts_data)


#=========replace rownames with  gene ids

rownames(counts_data) <- counts_data$GeneID

#==========remove column 1(gene id column)

counts_data <- counts_data[, -1]
rownames(counts_data)


#========read in datainfo: data_info contains a metadata on cellines and samples info(treated and untreated) 
colData <- import("C:/Users/NANA ADOMAKO ANSAH/OneDrive/Documents/data_info.xlsx")

row.names(colData) <- colData$sample_id
colData <- colData[, -1]
colData


#==========================checking if all colnames of counts_data is same as rownmaes of coldata

all(colnames(counts_data) %in% rownames(colData))

#========================checking if they are in the same order

all(colnames(counts_data) == rownames(colData))



#=========================construct a DESeq dataset

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ Dexamethasone)

dds

#=========================setting a factor level: In this we use "untreated" as a ref

dds$Dexamethasone <-  relevel(dds$Dexamethasone, ref = 'untreated')



# =========================pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


#=================================Run DESeq: DESeq runs on raw read counts
dds <- DESeq(dds) 




#=============================assigning results to res
res <- results(dds, alpha = 0.001)

res
summary(res)

#=================================MA plot
plotMA(res,
  main = "MA Plot-Differentially Expressed Genes",
  colNonSig = "gray60",colSig = "darkred",
  colLine = "grey40")


rownames(res)
sum(is.na(res$padj)) #cecking the number of NA's in Padj column 

# ================Create background dataset for hypergeometric testing using all genes tested for significance in the results 

all_genes <- as.character(rownames(res))

# Extract significant results
signif_res <- res[res$padj < 0.001 & !is.na(res$padj),]

summary(signif_res)
signif_genes <- as.character(rownames(signif_res))
signif_genes


#====================GO ANALYSIS

GO_results_BP <- enrichGO(gene = signif_genes, 
                          universe = all_genes,
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          readable = TRUE)


GO_results_CC <- enrichGO(gene = signif_genes, 
                          universe = all_genes,
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENTREZID", 
                          ont = "CC",
                          pAdjustMethod = "BH",
                          readable = TRUE)

as.data.frame(GO_results_BP)
as.data.frame(GO_results_CC)


barplot(GO_results_BP, showCategory = 20)
barplot(GO_results_CC, showCategory = 20)

#==============================KEGG pathway enrichment

kegg_enrich <- enrichKEGG(gene = signif_genes, organism = 'hsa')

#Viz
dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enriched Terms")


#============Taking a closer look at the pathway enriched terms
#cnet uses entrezids only

edo <- enrichKEGG(signif_genes)

edo_x <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edo_x, showCategory = 5,
               categorySize = "pvalue",
               node_label = "category",
               color_category='darkred', 
               color_gene='steelblue')
p1

p2 <- cnetplot(edo_x, 
               showCategory = 5,
               node_label = "gene",
               color_category = 'darkred', 
               color_gene = 'steelblue',
               cex_label_category = 0.3,  
               cex_label_gene = 0.6) 
p2

p3 <- cnetplot(edo_x, showCategory = 5, circular = TRUE, 
               colorEdge = TRUE, cex_label_gene = 0.5,
               cex_label_category = NULL)
p3
