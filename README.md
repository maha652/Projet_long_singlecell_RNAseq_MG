I ) Project's Title : 

Projet_long_singlecell_RNAseq_MG

II) Project Description : 

Perform single-cell RNA-seq analysis from "Seurat Object" to differentially expressed genes {markers_diff file , volcano_plot } ,Cell_proportions { dotplot ,barplot ,ggplot} and Expressions_plots {Featureplot}.

Available Steps : 

1)Differentially expressed genes(analyse_differentielle_gsea_enrichissement.R):


- markers_diff:csv file containing differentially expressed genes

- volcano_plot  : represents the differentially expressed genes

- dotplot : represents GSEA analysis results

- barplot : represents enrichment analysis result 

2) Cell_proportions(proportions_cellulaires.R ):

- Dataframe : cell proportions dataframe

- ggplot : Cell proportions graph

3) Expressions_plots(tracés_expressions.R):

- Featureplot : Gene expression on UMAPs


III) Runing

You should have your file "Seurat Object " resulting from sequencing of single cell suspension already transformed from a FastQ file as an input and run by calling  : 

- analyse_differentielle_gsea_enrichissement.R

- proportions_cellulaires.R

- tracés_expressions.R 

IV) Added folder 

Params_pipeline_BIGR_et_scripts_bash : 

- A folder containing the parameters to be applied during the various stages of the single-cell pipeline developed by the bioinformaticians on Gustave Roussy's bioinformatics platform for the three types of analysis: 
- Individual 
- Grouped
- Integrated 

and the corresponding bash scripts to launch each parameter file.
For more details about the pipeline "single-cell analysis" , please find the wiki link :  
https://github.com/gustaveroussy/single-cell/wiki










