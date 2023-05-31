# But du script:
# Ce script a pour but de réaliser une analyse différentielle entre les conditions des 2 échantillons pour un cluster d'intérêt

# Installation des packages 
if(!require("pacSeuratman")) install.packages("Seurat")
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano")
if(!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if(!require("msigdbr")) install.packages("msigdbr")
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if(!require("stats")) BiocManager::install("stats")
# Chargement des libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(stats)
# Paramètres
input_data <-  "....projet_long/CB50_CTRL_MIF_Int_Seurat_NORMKEPT_pca_17_0.7.rda"
output_data <- "....../projet_long/output/"
clusters_to_select <- "10"  # 10 = Monocytes; 1 et 4 = Granulocytes

# Chargement des données
load(input_data)
ctrl_mif <- sobj
rm(sobj)
gc()


###########################################
########  Analyse différentielle  #########
###########################################

# Comme les échantillons ont été normalisés séparément avec la fonction SCTransform, il faut rendre leur distribution
# comparable.
ctrl_mif <- PrepSCTFindMarkers(ctrl_mif) #uniformisation des normalisation

# Seurat propose une fonction (FindMarkers) pour comparer les clusters, mais nous souhaitons comparer les échantillons 
# sur un cluster précis.
# Donc il faut sélectionner les cellules sur ce cluster puis, dans l'objet seurat, remplacer l'information du clustering
# par l'information de l'échantillon d'origine pour pouvoir utiliser cette fonction.

# Sous-échantillonnage du cluster d'intérêt 
ctrl_mif_selected <- subset(x= ctrl_mif , idents = clusters_to_select)

# Inversement des information de clustering et de l'échantillon d'origine
head(Idents(ctrl_mif_selected))                  #clustering initial
Idents(ctrl_mif_selected) <- "orig.ident"        #inversement
head(Idents(ctrl_mif_selected))                  #vérification

# Analyse différentielle
markers_diff <- FindMarkers(ctrl_mif_selected, ident.1 = "CB50_MIF_GE" , ident.2 = "CB50_CTRL_GE", logfc.threshold = 0, test.use = "wilcox")

# Sauvegarde de la table de résultats
write.csv2(markers_diff, file = paste0(output_data, "markers_diff_clust_",paste(clusters_to_select, collapse = "_"),".csv"), row.names = TRUE)

# Volcanoplot
volcano_plot <- EnhancedVolcano(markers_diff,
                                lab = rownames(markers_diff),
                                x = "avg_log2FC",
                                y = "p_val_adj", #toujours prendre la p-valeur ajustée!
                                labSize = 3,
                                title = "Volcano Plot des gènes différentiellement exprimés",
                                subtitle ="Seuil p-values = 0.05, Seuils LogFC = -0.58 et 0.58",
                                ylab = bquote(~-Log[10]~(adjusted~P-values)),
                                axisLabSize=14,
                                pCutoff = 0.05,
                                FCcutoff = 0.58)

# Sauvegarde du volcanoplot
print(volcano_plot)
ggsave(paste0(output_data, "Volcanoplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),".png"), plot = volcano_plot, width = 10, height = 10)





############################################
#########  GSEA et Enrichissement  #########
############################################

# Parameters
ref_species <- org.Hs.eg.db
OrgDb_species <- "org.Hs.eg.db"

# Récupération des voies de signalisation et des gènes correspondants issus de la base de données MSigDB
m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, entrez_gene)
m_t2g_C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>% dplyr::select(gs_name, entrez_gene)
m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, entrez_gene)


### GSEA ###

# Identification des gènes significatifs (seuil sur la pvaleur)
All_geneSig = subset(markers_diff, markers_diff$p_val_adj < 0.05)
All_geneSig$SYMBOL = rownames(All_geneSig)

# Traduction des noms de gène en code ENTREZID (recommandé)
eg = tryCatch(bitr(All_geneSig$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb_species), error=function(err){ print(err); data.frame(SYMBOL=NA, ENTREZID=NA)})

# Remplacement des noms des gènes par le code ENTREZID
All_geneSig = merge(All_geneSig, eg, by ="SYMBOL") %>% dplyr::select(ENTREZID,avg_log2FC)
#formattage et tri par ordre décroissant de avg_log2FC pour All_geneSig
All_geneSig_gsea <- All_geneSig$avg_log2FC
names(All_geneSig_gsea) <- as.character(All_geneSig$ENTREZID)
All_geneSig_gsea <- sort(All_geneSig_gsea, decreasing = TRUE)

# GSEA
gsea_msig_H <- GSEA(All_geneSig_gsea, TERM2GENE = m_t2g_H, maxGSSize = Inf, verbose = FALSE, seed = 1234 , pvalueCutoff = 0.25)
gsea_msig_C2 <- GSEA(All_geneSig_gsea, TERM2GENE = m_t2g_C2, maxGSSize = Inf, verbose = FALSE, seed = 1234 , pvalueCutoff = 0.25)
gsea_msig_C6 <- GSEA(All_geneSig_gsea, TERM2GENE = m_t2g_C6, maxGSSize = Inf, verbose = FALSE, seed = 1234 , pvalueCutoff = 0.25)
gsea_msig_C7 <- GSEA(All_geneSig_gsea, TERM2GENE = m_t2g_C7, maxGSSize = Inf, verbose = FALSE, seed = 1234 , pvalueCutoff = 0.25)


# Remplacement des codes ENTREZID par les noms de gènes
gsea_msig_H <- setReadable(gsea_msig_H, ref_species, keyType = "ENTREZID")
gsea_msig_C2 <- setReadable(gsea_msig_C2, ref_species, keyType = "ENTREZID")
gsea_msig_C6 <- setReadable(gsea_msig_C6, ref_species, keyType = "ENTREZID")
gsea_msig_C7 <- setReadable(gsea_msig_C7, ref_species, keyType = "ENTREZID")

# Tracés des graphes
dotplot_H <- if(length(gsea_msig_H$Description)>1) dotplot(gsea_msig_H, showCategory=20) else  NULL
dotplot_C2 <- if(length(gsea_msig_C2$Description)>1) dotplot(gsea_msig_C2, showCategory=20) else  NULL
dotplot_C6 <- if(length(gsea_msig_C6$Description)>1) dotplot(gsea_msig_C6, showCategory=20) else  NULL
dotplot_C7 <- if(length(gsea_msig_C7$Description)>1) dotplot(gsea_msig_C7, showCategory=20) else  NULL

# Savegarde des graphes
if(!is.null(dotplot_H)) ggsave(paste0(output_data, "Dotplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_H.png"), plot = dotplot_H, width = 10, height = 15)
if(!is.null(dotplot_C2)) ggsave(paste0(output_data, "Dotplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C2.png"), plot = dotplot_C2, width = 10, height = 15)
if(!is.null(dotplot_C6)) ggsave(paste0(output_data, "Dotplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C6.png"), plot = dotplot_C6, width = 10, height = 15)
if(!is.null(dotplot_C7)) ggsave(paste0(output_data, "Dotplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C7.png"), plot = dotplot_C7, width = 10, height = 15)


### Enrichissement ###

#Identification des gènes différentiellement exprimés (seuil sur la pvaleur et le LogFC)
All_geneDE = subset(markers_diff, (markers_diff$p_val_adj < 0.05) & (abs(markers_diff$avg_log2FC) > 0.58))
All_geneDE$SYMBOL = rownames(All_geneDE)

# Remplacement des noms des gènes par le code ENTREZID
All_geneDE <-merge(All_geneDE, eg, by ="SYMBOL") %>% dplyr::select(ENTREZID,avg_log2FC)

# Enrichissement
enrich_msig_H <- enricher(All_geneDE$ENTREZID, TERM2GENE=m_t2g_H)
enrich_msig_C2 <- enricher(All_geneDE$ENTREZID, TERM2GENE=m_t2g_C2)
enrich_msig_C6 <- enricher(All_geneDE$ENTREZID, TERM2GENE=m_t2g_C6)
enrich_msig_C7 <- enricher(All_geneDE$ENTREZID, TERM2GENE=m_t2g_C7)

# Remplacement des codes ENTREZID par les noms de gènes
enrich_msig_H <- setReadable(enrich_msig_H, ref_species, keyType = "ENTREZID")
enrich_msig_C2 <- setReadable(enrich_msig_C2, ref_species, keyType = "ENTREZID")
enrich_msig_C6 <- setReadable(enrich_msig_C6, ref_species, keyType = "ENTREZID")
enrich_msig_C7 <- setReadable(enrich_msig_C7, ref_species, keyType = "ENTREZID")

# Tracés des graphes
barplot_H <- if(length(enrich_msig_H$Description)>1) barplot(enrich_msig_H, showCategory=20) else  NULL
barplot_C2 <- if(length(enrich_msig_C2$Description)>1) barplot(enrich_msig_C2, showCategory=20) else  NULL
barplot_C6 <- if(length(enrich_msig_C6$Description)>1) barplot(enrich_msig_C6, showCategory=20) else  NULL
barplot_C7 <- if(length(enrich_msig_C7$Description)>1) barplot(enrich_msig_C7, showCategory=20) else  NULL

# Savegarde des graphes
if(!is.null(barplot_H)) ggsave(paste0(output_data, "Barplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_H.png"), plot = barplot_H, width = 10, height = 15)
if(!is.null(barplot_C2)) ggsave(paste0(output_data, "Barplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C2.png"), plot = barplot_C2, width = 10, height = 15)
if(!is.null(barplot_C6)) ggsave(paste0(output_data, "Barplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C6.png"), plot = barplot_C6, width = 10, height = 15)
if(!is.null(barplot_C7)) ggsave(paste0(output_data, "Barplot_ctrl_mif_clust_",paste(clusters_to_select, collapse = "_"),"_C7.png"), plot = barplot_C7, width = 10, height = 15)


