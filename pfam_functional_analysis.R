
# Analysis of the domains detected by Pfam software. 
# Venn diagrams and upset plots are constructed for each orthogroup, to analyze differences in protein domains.
# In the end of this code a Gene Ontology (GO) analysis is done for each species. However, this kind of approach is not useful for orthologues, and as expected, no differences are seen between them.

library(readr)
OG <- read_table("~/TFM/Pfam analysis/pfam_outputs/OG20.txt", 
                 col_names = FALSE)
View(OG)

categories <- unique(OG$X1)

drosophila <- elegans <- mouse <- human <- vector()

for (category in categories) {
  subset_df <- subset(OG, X1 == category)
  unique_values <- unique(subset_df$X7)
  if (category == "C25G4.2.1") { #susbtitute here elegans gene name
    elegans <- c(elegans, unique_values)
  } else if (category == "FBpp0083317") { #susbtitute here drosophila gene name
    drosophila <- c(drosophila, unique_values)
  } else if (category == "ENSMUSP00000071493.3") { #susbtitute here mouse gene name
    mouse <- c(mouse, unique_values)
  } else if (category == "ENSP00000461388.1") { #susbtitute here human gene name
    human <- c(human, unique_values)
  }
}

library(zoo)

# Find the maximum length among the vectors
max_length <- max(length(drosophila), length(elegans), length(mouse), length(human))

# Fill the vectors with NA values to match the maximum length
drosophila <- `length<-`(drosophila, max_length)
elegans <- `length<-`(elegans, max_length)
mouse <- `length<-`(mouse, max_length)
human <- `length<-`(human, max_length)

# Create the data frame with the vectors

df<- data.frame(domains_elegans=elegans, domains_drosophila=drosophila, 
                domains_mouse=mouse, domains_human=human)

elegans <- na.omit(c(df$domains_elegans))
drosophila <- na.omit(c(df$domains_drosophila))
mouse <- na.omit(c(df$domains_mouse))
human <- na.omit(c(df$domains_human))


## Venn diagram for representing shared and new/lost domains.

library(VennDiagram)
library(RColorBrewer)

myCol <- brewer.pal(4, "Spectral")

venn<-venn.diagram(
  x = list(elegans, drosophila, mouse, human),
  filename = 'OG5_shared_domains_venn.png',
  fill = myCol,
  cat.cex = 0.8,
  category.names = c("C. elegans" , "D. melanongaster " , "M. musuclus", "H. sapiens"),
  output=TRUE,
  main = "OG5 shared domains"
)

# Upset plots for representing shared and new/lost domains.

library(UpSetR)

x=list("C. elegans"= elegans,"D. melanogaster" = drosophila, 
       "M. musculus" = mouse, "H. sapiens" = human)

upsetplot<-upset(fromList(x), order.by="freq", 
                 sets.bar.color = c("lightgreen","lightblue", "darkgreen", "purple"),
                 text.scale =  2)

upsetplot


df_subset <- head(selected_orthologs, 10)


file_name <- "C:/Users/toric/OneDrive/Documentos/TFM/table.txt"
file_conn <- file(file_name)
write.table(names(df_subset), file_conn, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(df_subset, file_conn, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
close(file_conn)

table <- readLines(file_name)


cat(table, sep="\n")

# Extraction of Gene Symbols from Ensembl IDs, as these could be more informative. 
# I use biomaRt package to retrieve gene names given these ENSEMBL IDs. 
# In the case of ensembl IDs of human and mouse, the isoform is also given, which hampers the retrieval of gene identifiers. 
# For this reason the part of the name after the dot is removed. However this would be taken into account later on in specific orthologs analysis.

library(biomaRt)
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl_connect_drosophila <- useMart("ensembl", dataset ='dmelanogaster_gene_ensembl')
ensembl_connect_elegans <- useMart("ensembl", dataset ='celegans_gene_ensembl')

genes_drosophila <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                          filters = "ensembl_peptide_id",
                          values = selected_orthologs$names_drosophila,
                          mart = ensembl_connect_drosophila)
genes_elegans <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                       filters = "ensembl_peptide_id",
                       values = selected_orthologs$names_elegans,
                       mart = ensembl_connect_elegans)

genes_elegans
# I delete the isoform part of the ensembl IDs for human and mouse proteins
mouse <- c()
spl <- strsplit(x = selected_orthologs$names_mouse, split = "[.]")

for (i in 1:length(spl)) {
  mouse <- c(mouse, spl[[i]][1])
  
}

human <- c()
spl <- strsplit(x = selected_orthologs$names_human, split = "[.]")

for (i in 1:length(spl)) {
  human <- c(human, spl[[i]][1])
  
}

ensembl_connect_mouse <- useMart("ensembl", dataset ='mmusculus_gene_ensembl')
ensembl_connect_human <- useMart("ensembl", dataset ='hsapiens_gene_ensembl')

genes_mouse <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                     filters = "ensembl_peptide_id",
                     values = mouse,
                     mart = ensembl_connect_mouse)
genes_human <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"),
                     filters = "ensembl_peptide_id",
                     values = human,
                     mart = ensembl_connect_human)

# I add the gene IDs to the selected_orthologs first dataframe
genes_human <- genes_human[order(genes_human$ensembl_peptide_id),]
selected_orthologs <- selected_orthologs[order(selected_orthologs$names_human),]

selected_orthologs <- cbind(selected_orthologs, genes_id_human = genes_human$external_gene_name) # add human gene IDs 

genes_mouse <- genes_mouse[order(genes_mouse$ensembl_peptide_id),]
selected_orthologs <- selected_orthologs[order(selected_orthologs$names_mouse),]

selected_orthologs <- cbind(selected_orthologs, genes_id_mouse = genes_mouse$external_gene_name) # add mouse gene IDs 

genes_drosophila <- genes_drosophila[order(genes_drosophila$ensembl_peptide_id),]
selected_orthologs <- selected_orthologs[order(selected_orthologs$names_drosophila),]

selected_orthologs <- cbind(selected_orthologs, genes_id_drosophila = genes_drosophila$external_gene_name) # add drosophila gene IDs 

genes_elegans <- genes_elegans[order(genes_elegans$ensembl_peptide_id),]
selected_orthologs <- selected_orthologs[order(selected_orthologs$names_elegans),]

selected_orthologs <- cbind(selected_orthologs, genes_id_elegans = genes_elegans$external_gene_name) # add elegans gene IDs 

# Export of the selected orthologs 
library(data.table)

directory <- "~/TFM/Pfam analysis/pfam_outputs"

file_list <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)

dataframes <- list()

for (file in file_list) {
  dt <- fread(file, header = FALSE, sep = "\u0020", fill = TRUE, na.strings = "")
  file_name <- tools::file_path_sans_ext(basename(file))
  dt[, doc_name := file_name]
  dataframes[[file_name]] <- dt
}

combined_df <- rbindlist(dataframes)
library(writexl)
write_xlsx(combined_df, "C:/Users/toric/combined.xlsx")

# A Gene Ontology analysis if performed. However,Gene Ontology is based on sequence homology, which would probably hamper the analysis of differences in function between orthologs. 
# First databases for each species are loaded.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

BiocManager::install("org.Hs.eg.db", force=TRUE)
library(org.Hs.eg.db)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celegans.db")
library(celegans.db)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)


GO analysis for Biological Process (BP) is performed for each species.

GO_human <- enrichGO(gene          = selected_orthologs$genes_id_human,
                     OrgDb        = 'org.Hs.eg.db', 
                     keyType      = "SYMBOL",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable     = TRUE)
GO_human<- as.data.frame(GO_human)

GO_mouse <- enrichGO(gene          = selected_orthologs$genes_id_mouse,
                     OrgDb        = 'org.Mm.eg.db', 
                     keyType      = "SYMBOL",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable     = TRUE)
GO_mouse<- as.data.frame(GO_mouse)

GO_drosophila <- enrichGO(gene          = selected_orthologs$genes_id_drosophila,
                          OrgDb        = 'org.Dm.eg.db', 
                          keyType      = "SYMBOL",
                          ont          = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable     = TRUE)
GO_drosophila<- as.data.frame(GO_drosophila)

GO_elegans <- enrichGO(gene          = selected_orthologs$genes_id_elegans,
                       OrgDb        = 'celegans.db', 
                       keyType      = "SYMBOL",
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable     = TRUE)
GO_elegans<- as.data.frame(GO_elegans)

