
# Here, the steps needed for an ancestral reconstruction analysis are shown, using ape package.
# The aim of this code is to obtain the length of the ancestral protein sequnce for each node in each orthogroup selected in the orthologue_analysis.R code.
library(ape)
library(phytools)
library(seqinr)
library(treeio)
library(ggtree)

# Species tree representation
species_tree<-read.tree(text="(CE:0.350896,(DM:0.540653,(HS:0.0267447,MM:0.0702978)0.972764:0.353361)1:0.350896);")

plot(species_tree, type = "phylogram", show.tip.label = TRUE, cex = 1)



nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)


ggtree(species_tree, color="darkgreen", size=2, layout="roundrect") + 
  geom_tiplab() + geom_rootedge() + theme(plot.margin = margin(0.2, 0.2, 0.01, 0.2, "cm"))


# phylogenetic trees for orthogroups are loaded

tree <-read.tree(text="(C._elegans_dolk-1_299_aa:0.429625,(D._melanogaster_CG8311_503_aa:0.857442,(H._sapiens_DOLK_534_aa:0.027903,M._musculus_Dolk_538_aa:0.035346)n2:0.916521)n1:0.429625)n0;")

caracters <- c(299, 503, 534, 538)

ace <-ace(caracters, tree)

reconstructed_states <- ace$ace

plot(tree, type = "phylogram", show.tip.label = TRUE, cex = 1)
nodelabels(round(reconstructed_states, 2), cex = 0.8)

stop_codon_loss <- selected_orthologs[
  selected_orthologs$genes_human > selected_orthologs$genes_mouse |
    selected_orthologs$genes_mouse > selected_orthologs$genes_drosophila |
    selected_orthologs$genes_drosophila > selected_orthologs$genes_elegans,
]

# Create a new dataframe for stop codon gain or loss
stop_codon_gain <- selected_orthologs[
  !(selected_orthologs$genes_human > selected_orthologs$genes_mouse |
      selected_orthologs$genes_mouse > selected_orthologs$genes_drosophila |
      selected_orthologs$genes_drosophila > selected_orthologs$genes_elegans),
]




