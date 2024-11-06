## INSTAL.LACIÓ PAQUETS ##

# Configuració de mirror de CRAN
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Instal.lació de paquets
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
install.packages("readxl") #Readxl
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("tidyr")
install.packages("pheatmap")

## PRIMERA VISUALITZACIÓ ##

# Carreguem readxl per poder veure l'arxiu original
library(readxl)

# Definim la ruta de l'arxiu
file_path <- "C:/Users/Arnau/Desktop/PEC1_ADO/original_data/TIO2+PTYR-human-MSS+MSIvsPD.XLSX"

# Carreguem les dades com original_data
original_data <- read_excel(file_path, sheet = "originalData")

# Mostrar las primeras 5 filas
head(original_data, 5)

## CREACIÓ SUMMARIZEDEXPERIMENT I PREPROCESSING ##

# Carreguem els paquets necessaris
library(SummarizedExperiment)
library(readxl)

# Un cop vist el format original de les dades, seleccionem les que ens interessen i generem data_filtered
selected_columns <- c("SequenceModifications", 
                      "M1_1_MSS", "M1_2_MSS", "M5_1_MSS", "M5_2_MSS", 
                      "T49_1_MSS", "T49_2_MSS", 
                      "M42_1_PD", "M42_2_PD", "M43_1_PD", "M43_2_PD", 
                      "M64_1_PD", "M64_2_PD")
data_filtered <- original_data[selected_columns]

# Generem la matriu assay amb les abundàncies de fosfopèptids
assay_data <- as.matrix(data_filtered[, -1]) # Treiem la columna SequenceModifications, només volem les dades a assay
rownames(assay_data) <- data_filtered$SequenceModifications

# Definim colData manualment per fer coincidir amb assay_data
samples <- c("M1", "M1", "M5", "M5", "T49", "T49", "M42", "M42", "M43", "M43", "M64", "M64")
phenotypes <- c("MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "PD", "PD", "PD", "PD", "PD", "PD")
col_data <- DataFrame(Sample = samples, Phenotype = phenotypes)
rownames(col_data) <- colnames(assay_data)

# Afegim a row_data les metadades de SequenceModifications com informació addicional dels fosfopèptids
row_data <- DataFrame(SequenceModifications = data_filtered$SequenceModifications)
rownames(row_data) <- data_filtered$SequenceModifications

# Un cop format Assay, colData i rowData, generem el SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = assay_data),
  rowData = row_data,
  colData = col_data
)

# Visualitzem assay del nostre SummarizedExperiment per veure que les dades hi són, igual que les metadades de les columnes i les files
se
head(assay(se), 5)

## DATA VISUALIZATION ##

library(ggplot2)
library(ggrepel)
library(tidyr)
library(pheatmap)

summary(assay(se))

# Definim els colors que volem que s'utilitzin per als plots
custom_colors <- c("MSS" = "#90008B", "PD" = "#90DD77")

# Preparem les dades pel boxplot
# Utilitzant tidyr, passem Assay a un data frame per poder generar els boxplots.
boxplot_data <- as.data.frame(assay_data)

# Degut a la naturalesa de ggplot2, paquet que utilitzarem per generar els boxplots, hem de passar les dades a un format "long", ja que així facilitem el mapeig de variables estétiques (color, fill, tema...), de mateixa manera que així també podrem agrupar les dades per grups (fenotips, en aquest cas)
boxplot_data_long <- pivot_longer(boxplot_data, cols = everything(), names_to = "Sample", values_to = "Abundance")

# Creem un vector de Phenotype que coincideixi amb cada mostra generant els dos grups segons fenotip.
phenotype_map <- setNames(phenotypes, colnames(assay_data))
boxplot_data_long$Phenotype <- phenotype_map[boxplot_data_long$Sample]

# Ordenem els fenotips per a que apareguin així als boxplots
sample_order <- c(colnames(assay_data)[phenotypes == "MSS"], colnames(assay_data)[phenotypes == "PD"])
boxplot_data_long$Sample <- factor(boxplot_data_long$Sample, levels = sample_order)

# 1. Boxplot d'abundàncies sense transformació logarítmica
ggplot(boxplot_data_long, aes(x = Sample, y = Abundance, fill = Phenotype)) +
  geom_boxplot(aes(color = Phenotype), outlier.shape = 21, outlier.size = 1, outlier.stroke = 1, outlier.fill = NA) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Phosphoproteomics Experiment - Normalized Abundances", x = "Sample", y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

# 2. Boxplot logarítmic
# Creem una nova columna para les abundàncies logarítmiques. Afegim el + 1 al final per evitar problemes amb els 0s, tenint en compte que treballem en escala log10.
boxplot_data_long$LogAbundance <- log10(boxplot_data_long$Abundance + 1)

ggplot(boxplot_data_long, aes(x = Sample, y = LogAbundance, fill = Phenotype)) +
  geom_boxplot(outlier.colour = "red", outlier.fill = NA, outlier.shape = 21, outlier.size = 1) +  # Outliers vermells
  scale_fill_manual(values = custom_colors) + 
  labs(title = "Phosphoproteomics Experiment - Normalized Abundances in Log10 Scale", x = "Sample", y = "Log10(Abundance + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. PCA
# Generem el PCA amb les dades logarítmiques
pca_res <- prcomp(t(log10(assay_data + 1)), scale. = FALSE)
pca_data <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Phenotype = phenotypes)
pca_data$Sample <- colnames(assay_data)  # Cridem a assay_data per poder identificar els dots dins el PCA

# Creem el plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(title = "PCA of Phosphoproteomic Data", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, color = "grey50", linetype = "solid") +
  geom_vline(xintercept = 0, color = "grey50", linetype = "solid") +
  theme_minimal()

# 4. Heatmap + dendograma en els eixos
# Creem el log d'Assay i eliminem els NA, en cas de que n'hi hagin.
log_data <- log10(assay_data + 1)
log_data[is.na(log_data)] <- 0

# Configurem les columnes de les mostres
annotation_col <- data.frame(Phenotype = phenotypes)
rownames(annotation_col) <- colnames(log_data)

# Assignem els colors que ja veniem utilitzant per a cada fenotip
annotation_colors <- list(Phenotype = custom_colors)

# Creem el heatmap
pheatmap(log_data, 
         annotation_col = annotation_col,
         annotation_colors = list(Phenotype = custom_colors),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Heatmap of Phosphoproteomic Data",
         show_rownames = FALSE)

# Exportació de totes les dades en format text
# Directori on guardar els CSV
output_dir <- "C:/Users/Arnau/Desktop/PEC1_ADO/data/"

# 1. Dades filtrades
write.csv(data_filtered, file.path(output_dir, "data_filtered.csv"), row.names = FALSE)

# 2. assay_data
write.csv(as.data.frame(assay_data), file.path(output_dir, "assay_data.csv"), row.names = TRUE)

# 3. colData
write.csv(as.data.frame(col_data), file.path(output_dir, "col_data.csv"), row.names = TRUE)

# 4. rowData
write.csv(as.data.frame(row_data), file.path(output_dir, "row_data.csv"), row.names = TRUE)

# 5. boxplot_data_long
write.csv(boxplot_data_long, file.path(output_dir, "boxplot_data_long.csv"), row.names = FALSE)

# 6. pca_data
write.csv(pca_data, file.path(output_dir, "pca_data.csv"), row.names = FALSE)

# 7. log_data
write.csv(as.data.frame(log_data), file.path(output_dir, "log_data.csv"), row.names = TRUE)

# 8. annotation_col
write.csv(annotation_col, file.path(output_dir, "annotation_col.csv"), row.names = TRUE)