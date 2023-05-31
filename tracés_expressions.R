# But du script:
# Ce script a pour but de tracer l'expression de certains gènes sur les UMAP, afin de confirmer l'annotation des clusters.

# Installation des packages 
if(!require("Seurat")) install.packages("Seurat")
if(!require("ggplot2")) install.packages("ggplot2")

# Chargement des libraries
library(Seurat)
library(ggplot2)

# Paramètres
input_data_int <-".../projet_long/CB50_CTRL_MIF_Int_Seurat_NORMKEPT_pca_17_0.7.rda"
output_data <- ".../output/traces_expression/"
#genes_list <- c("CD14","FUT4","TFRC","CD34") # CD14 = monocytes; FUT4 = granulocytes; TFRC = érythroïdes; CD34 = cellules souches hématopoïétiques.
#autres listes de gènes issues de https://www.proteinatlas.org/ et https://www.cell.com/cms/10.1016/j.molcel.2020.04.008/attachment/fe0af298-d8d0-49dc-b49b-bf44fec66827/mmc1.pdf
#genes_list <- c("EPOR","KLF1","TFR2","CSF2RB","APOE","APOC1","CNRIP1","CD38","HBB") # Erythoïdes
#genes_list <- c("MPIG6B","PF4","GP9","VWF","SELP") # Megakaryocytes
#genes_list <- c("CRHBP","HLF","AVP") # HSPC souches
#genes_list <- c("ELANE","AZU1","PRTN3","MPO","CST7","CTSG","CEACAM8") # Neutrophiles
#genes_list <- c("CFD") # Granulocytes
#genes_list <- c("CD14", "CCR2","CSF1R") # Monocytes
#genes_list <- c("CPA3","MS4A2") # Basophiles
genes_list <- c("MBP" , "ECP") # Eosinophiles #ces genes n'ont pas �t� trouv�
# Chargement des données
load(input_data_ctrl)
ctrl <- sobj
load(input_data_mif)
mif <- sobj
load(input_data_int)
ctrl_mif <- sobj
rm(sobj)
gc()

# Traçage des graphes
featureplot_ctrl <- FeaturePlot(ctrl, features = genes_list)
featureplot_mif <- FeaturePlot(mif, features = genes_list)
featureplot_ctrl_mif <- FeaturePlot(ctrl_mif, features = genes_list)

# Sauvegarde des graphes
ggsave(paste0(output_data,"Featureplot_ctrl.png"), plot = featureplot_ctrl, width = 10, height = 10)
ggsave(paste0(output_data,"Featureplot_mif.png"), plot = featureplot_mif, width = 10, height = 10)
ggsave(paste0(output_data,"Featureplot_ctrl_mif.png"), plot = featureplot_ctrl_mif, width = 10, height = 10)



