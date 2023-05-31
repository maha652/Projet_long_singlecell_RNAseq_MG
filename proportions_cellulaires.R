# But du script:
# Ce script a pour but de mesurer les proportions de cellules entre les conditions

# Installation des packages 
if(!require("pacSeuratman")) install.packages("Seurat")
if(!require("ggplot2")) install.packages("ggplot2")

# Chargement des libraries
library(Seurat)
library(ggplot2)

# Chargement des données
input_data <- ".../projet_long/CB50_CTRL_MIF_Int_Seurat_NORMKEPT_pca_17_0.7.rda"
load(input_data)
ctrl_mif <- sobj
rm(sobj)
gc()

# Renommage des clusters
ctrl_mif <- RenameIdents(object = ctrl_mif,
                         '1' = "Granulocytes", 
                         '2' = "Erythrocytes",
                         '4' = "Granulocytes",
                         '3' = "Progéniteurs",
                         '0' = "Progéniteurs",
                         '5' = "Progéniteurs",
                         '11' = "Progéniteurs",
                         '7' = "Progéniteurs",
                         '6'= "Erythrocytes",
                         '9' = "Erythrocytes",
                         '12' ="Erythrocytes" ,
                         '13' = "Erythrocytes",
                         '8' = "Erythrocytes",
                         '14' = "megacaryocytes",
                         '10' = "Monocytes")

# Création du dataframe des proportions cellulaires
pct_clust_by_sample <- data.frame(apply(table(Idents(ctrl_mif),ctrl_mif$orig.ident), 2, function(x){prop.table(x)*100}))
donnees <- data.frame(
  Echantillons = c(rep("CTRL",dim(pct_clust_by_sample)[1]), rep("MIF",dim(pct_clust_by_sample)[1])),
  Proportions = c(pct_clust_by_sample$CB50_CTRL_GE, pct_clust_by_sample$CB50_MIF_GE),
  Clusters =rep(rownames(pct_clust_by_sample),2)
)

# Graphique des proportions cellulaires
ggplot(donnees, aes(x = Echantillons, y = Proportions, fill = Clusters)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportions cellulaires entre ctrl et MIF", x = "Clusters", y = "Proportions") +
  theme_minimal()

# Test de Chi2 sur toutes les proportions cellulaires 
contingence_table_ctrl_mif <- table(Idents(ctrl_mif),ctrl_mif$orig.ident) # vérification que les effectifs sont bien >5
contingence_table_ctrl_mif
test_ctrl_mif <- chisq.test(contingence_table_ctrl_mif)
test_ctrl_mif

# Sous-échantillonnage des Monocytes et des Granulocytes
ctrl_mif_mono_granulo <- subset(x= ctrl_mif , idents = c("Monocytes","Granulocytes"))

# Test de Chi2 sur les proportions cellulaires des cellules différentiées des progénoteurs myéloïdes
contingence_table_ctrl_mif_mono_granulo <- table(Idents(ctrl_mif_mono_granulo),ctrl_mif_mono_granulo$orig.ident) # vérification que les effectifs sont bien >5
contingence_table_ctrl_mif_mono_granulo
test_ctrl_mif_mono_granulo <- chisq.test(contingence_table_ctrl_mif_mono_granulo)
test_ctrl_mif_mono_granulo


# Sous-échantillonnage des Monocytes , Granulocytes et erythrocytes 
ctrl_mif_mono_granulo_erythro <- subset(x= ctrl_mif , idents = c("Monocytes","Granulocytes" , "Erythrocytes"))

# Test de Chi2 sur les proportions cellulaires des cellules différentiées des progénoteurs myéloïdes
contingence_table_ctrl_mif_mono_granulo_erythro <- table(Idents(ctrl_mif_mono_granulo_erythro),ctrl_mif_mono_granulo_erythro$orig.ident) # vérification que les effectifs sont bien >5
contingence_table_ctrl_mif_mono_granulo_erythro
test_ctrl_mif_mono_granulo_erythro <- chisq.test(contingence_table_ctrl_mif_mono_granulo_erythro)
test_ctrl_mif_mono_granulo_erythro

