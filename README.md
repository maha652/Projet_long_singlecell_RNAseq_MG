I ) Project's Title : 

Projet_long_singlecell_RNAseq_MG

II) Project Description : 

Perform single-cell RNA-seq analysis from "Seurat Object" to markers_diff file , volcano_plot  , dotplot ,barplot,ggplot and Featureplot .

Available Steps : 

1) Analyse_differentielle_gsea_enrichissement (analyse_differentielle_gsea_enrichissement.R):


- markers_diff:csv file containing differentially expressed genes

- volcano_plot  : represents the differentially expressed genes

- dotplot : represents GSEA analysis results

- barplot : represents enrichment analysis result 

2) Proportions_cellulaires (proportions_cellulaires.R ):

- Dataframe : cell proportions dataframe

- ggplot : Cell proportions graph

3) Tracés_expressions (tracés_expressions.R):

- Featureplot : Gene expression on UMAPs


III) Runing

You should have your file "Seurat Object " resulting from sequencing of single cell suspension already transformed from a FastQ file by calling  : 

- analyse_differentielle_gsea_enrichissement.R

- proportions_cellulaires.R

- tracés_expressions.R 










