# stopcodonshift_evolution

## Introduction

This repository contains a collection of R files that aim to investigate the impact of stop codon shifts on molecular evolution. The analysis focuses on orthologous genes in eukaryotic genomes, and several important steps are performed to gain insights into this phenomenon. The main objective of this project is to understand the implications of stop codon shifts in molecular evolution. By analyzing orthologous genes in eukaryotic genomes, we aim to explore the relationship between stop codon shifts and evolution, as well as to investigate the potential consequences on protein function.

## Contents

This repository includes the following files:

**ortholog_analysis.R**: This script is performs ortholog analysis of an output file supplied by OrthoFinder (Emms and Kelly, 2019). It identifies orthologous genes and extracts information on protein length for further analysis.

**ancestral_reconstruction.R**: This script focuses on ancestral reconstruction analysis. It utilizes the obtained orthologous gene data to infer ancestral states of protein length using ape package (Analyses of Phylogenetics and Evolution).

**pfam_functional_analysis.R**: Here, an analysis of Pfam outputs is conducted to assess the potential functional consequences of stop codon shifts. This analysis provides insights into the domains associated with the orthologous genes.
