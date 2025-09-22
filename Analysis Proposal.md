# EE282 Bioinformatics Project Analysis Proposal

## Author: Bobby Sims

## Introduction

The role of T regulatory cells (Tregs) in auto immune neuro inflammatory diseases is unique and complex. The pathogenesis of MS for example, is characterized by their initial failure of Tregs in distinguishing self from non-self antigens, leading to immune reactivity to myelin, the structural component of neuron sheaths that facilitate neuron signal conduction. Yet paradoxically upon onset of auto immune neuroinflammation, it is the performance of Tregs that determines outcomes of the disease, as their multipronged immune suppressive capacities are essential to recovery and remission. The gatekeepers of immune tolerance can remediate their own initial failure. 
While image analysis and immune phenotyping provide the “what” and “who” of this phenomenon, differential gene expression is at the depth of analysis that facilitates get much closer to the “how” of the autoimmune neuro inflammatory event.  To achieve greater understanding of the paradoxical idiosyncrasies of Treg behavior, it is important to identify the key distinctions in transcriptional machinery that lead to differential disease outcomes between individuals and among Treg populations within individuals.
Comparison of genes differentially expressed between normal and EAE Treg transcriptomes aligns with and enriches our lab’s goal in the comprehensive characterization of the immune regulation capacities of Tregs in the context of auto immune neuroinflammation.
My analysis will be of the RNA microarray sequencing dataset GSE164460:
Transcriptome analysis of Tregs (CD4+Foxp3+) and Tcons (CD4+Foxp3-) from spleen and CNS of naïve mice and at the peak of EAE uploaded to NCBI Gene Expression Omnibus and originally published by Pohar et al, in “Antigen receptor-engineered Tregs inhibit CNS autoimmunity in cell therapy using nonredundant immune mechanisms in mice.” 
PMID: [35579560](https://pubmed.ncbi.nlm.nih.gov/35579560/)

## Methods

I will download software packages required for R analysis of microarray data from Affymetrix gene chip including “GEOquery” and “Affy”, from Bioconductor. 
Data “cleaning” includes unpacking .tar files into .CEL format and data normalization using RMA method. 
 I will analyze genes differentially expressed between wild type and EAE Tregs in CNS tissues (brain and spinal cord) and spleen tissue at the peak of disease in R using Deseq2/ limma. Variance between replicates will be shown on dispersion plot and I will analyze the top 10 and top 5 DEGs using volcano plot and PCA. Graphical visualizations of analyses are performed using ggplot2 package in R.  

## Rationale

Many of the necessary software tools and applications are available as open-source freeware. A few of the proposed analyses (ggplot2) have been previously covered in class as part of homework 3.
