# OrthoFinder can been installed using a Windows Subsystem for Linux (WSL) (Emms and Kelly, 2019).
# Orthogroupinference has been made with 4 species (H. sapiens, M. musculus, D. melanogaster and C. elegans).
# This R script is used to perform a descriptive analysis of the OrthoFinder results and to filter the selected orthogroups.

library(seqinr)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(ggtree)
# A descriptive analysis of the OrthoFinder results has been made. For making the phylogenetics trees the ggtree package of Bioconductor has been used, which reads trees in Newick format.
# 91.8% of the analyzed proteins have been assigned to orthogroups (230702 of 251303), which puts in evidence the good quality of the analysis.

setwd("C:/Users/toric/Documents/TFM/Results_Apr23_1/") 

orthogroups_statistics<-read.table("C:/Users/toric/Documents/TFM/Results_Apr23_1/Comparative_Geno
mics_Statistics/Statistics_Overall.tsv", header=T, stringsAsFactors = F)

# Ortogroups with a single protein per species have been chosen to facilitate de analysis. 
#These orthrogroups are presented by the OrthoFinder software in FASTA format. 

# Orthogroups with a single protein per species
# A function is built for opening the all the orthogroups files from the selected directory

abrir_archivos <- function(directorio, extension) {
  archivos <- list.files(directorio, pattern = paste0("*.", extension), 
                         full.names = TRUE)
  lista_contenidos <- list()
  for (archivo in archivos) {
    contenido <- read.fasta(archivo, as.string = TRUE, seqtype = "AA")
    lista_contenidos[[archivo]] <- contenido
  }
  return(lista_contenidos)
}
directorio <-"C:/Users/toric/Documents/TFM/Results_Apr23_1/Single_Copy_Orthologue_Sequ
ences/"
extension <- ".fa"
contenidos <- abrir_archivos(directorio, extension)

# Gene names are extracted for each orthogroup and each species from the generated list

# I choose the gene names for each specie
names_elegans <- vector()
names_drosophila <- vector()
names_mouse <- vector()
names_human <- vector()
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- names(lista_anidada_1)
  primer_valor <- lista_anidada_2[[1]]
  names_elegans <- append(names_elegans, primer_valor)
}
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- names(lista_anidada_1)
  primer_valor <- lista_anidada_2[[2]]
  names_drosophila <- append(names_drosophila, primer_valor)
}
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- names(lista_anidada_1)
  primer_valor <- lista_anidada_2[[4]]
  names_mouse <- append(names_mouse, primer_valor)
}
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- names(lista_anidada_1)
  primer_valor <- lista_anidada_2[[3]]
  names_human <- append(names_human, primer_valor)
}

# The same approach is taken to extract protein length from the FASTA sequences of the list

# The amino acid sequences in FASTA format are saved in vectors for each species.
genes_elegans <- vector()
genes_drosophila <- vector()
genes_mouse <- vector()
genes_human <- vector()
# Extract the sequences for each orthogroup for C. elegans
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- lista_anidada_1[[1]]
  primer_valor <- lista_anidada_2[[1]]
  genes_elegans <- append(genes_elegans, primer_valor)
}
# Extract the sequences for each orthogroup for D. melanogaster
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- lista_anidada_1[[2]]
  primer_valor <- lista_anidada_2[[1]]
  genes_drosophila <- append(genes_drosophila, primer_valor)
}
# Extract the sequences for each orthogroup for M. musculus
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- lista_anidada_1[[4]]
  primer_valor <- lista_anidada_2[[1]]
  genes_mouse <- append(genes_mouse, primer_valor)
}
# Extract the sequences for each orthogroup for H. sapiens
for (i in seq_along(contenidos)) {
  lista_anidada_1 <- contenidos[[i]]
  lista_anidada_2 <- lista_anidada_1[[3]]
  primer_valor <- lista_anidada_2[[1]]
  genes_human <- append(genes_human, primer_valor)
}

# Each amino acid sequence length is calculated with the nchar() function 
# These results together with gene names are saved in a table, which is then exported to a tab separated file (.tsv)
# Moreover, standard deviation (SD) is calculated, as an approximation of the difference in protein sequence length. 
# A minimum of 20 SD is established and the genes fulfilling this feature are exported to a tsv called selected_orthologs.

genes_elegans <-as.vector(genes_elegans)
genes_drosophila <-as.vector(genes_drosophila)
genes_mouse <-as.vector(genes_mouse)
genes_human <-as.vector(genes_human)
# Function for counting the character's lengths and gather the results in 
a table
almacenar_longitud_genes <- function(genes_elegans, genes_drosophila, 
                                     genes_mouse, genes_human) {
  nombres_genes <- c("genes_elegans", "genes_drosophila", "genes_mouse", 
                     "genes_human")
  resultados <- list()
  for (i in 1:length(nombres_genes)) {
    longitud <- nchar(get(nombres_genes[i]))
    resultados[[i]] <- longitud
  }
  
  # A table with the results is created
  tabla_resultados <- as.data.frame(do.call(cbind, resultados))
  colnames(tabla_resultados) <- nombres_genes
  return(tabla_resultados)
}
tabla_resultados <- almacenar_longitud_genes(genes_elegans, 
                                             genes_drosophila, genes_mouse, genes_human)
# Standard deviation for amino acid sequence length is extracted, as an 
approximation of differences between protein lengths
desv <- rowSds(as.matrix(tabla_resultados))
tabla_resultados <- cbind(tabla_resultados, desv)
View(tabla_resultados)
# Adding of gene names extracted above
tabla_resultados <- cbind(tabla_resultados, names_elegans, 
                          names_drosophila, names_mouse, names_human)

# Rows with an SD higher than 20 are selected and exported to a .tsv
selected_orthologs <- subset(tabla_resultados, desv>20)
as.numeric(selected_orthologs$genes_human)
as.numeric(selected_orthologs$genes_drosophila)
as.numeric(selected_orthologs$genes_elegans)
as.numeric(selected_orthologs$genes_mouse)

write.table(selected_orthologs, file =
              "C:/Users/toric/Documents/TFM/selected_orthologs.tsv", sep = "\t", quote 
            = FALSE, row.names = FALSE)

dim(selected_orthologs) # orthologs in which the length in very different among species


# 29 single gene per species orthogroups with very different protein length have been selected