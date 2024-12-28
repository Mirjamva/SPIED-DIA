###############################################################################
# Script Name:    CellLinePanel_Global_Norm.R
# Author:         Mirjam van Bentum
# Date:           2024-12-13
#
# Description:    This script performs preprocessing, normalization, and PCA 
#                 analyses of LFQ (Label-Free Quantification) data from timsTOF 
#                 measurements. It includes:
#                 - Loading and cleaning input data.
#                 - Applying filtering criteria and normalizing intensities.
#                 - Conducting PCA and generating PCA plots for multiple 
#                   cell lines and conditions.
#                 - Have ready: the output from DIA-NN (global analysis)
#                 - ran with dia-nn settings as outlined in readme
#                 - Raw files were searched against a HpH fractionated library. 
#
# Requirements:   - R environment with required packages installed:
#                     tidyverse, data.table, patchwork, limma, ggrepel
#                 - The custom functions defined here must be sourced or 
#                    included before running.
###############################################################################


#---------------------------
# Load Required Packages
#---------------------------
library(tidyverse)
library(data.table)
library(patchwork)
library(limma)
library(ggrepel)
library(stringr)

#---------------------------
# Define Directories
#---------------------------
results_dir <- "Results"
data_dir    <- "Data"
plots_dir   <- "Plots"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(data_dir))    dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plots_dir))   dir.create(plots_dir, recursive = TRUE)


#---------------------------
# Custom Functions
#---------------------------

fun_modseq_diannMQlib2pSTY <- function(seq_p) {
  # Convert phospho modifications to pS, pT, pY; remove SILAC, Oxidation, Acetyl, and UniMod tags
  str_replace_all(seq_p, c(
    "_" = "",
    "(S|T|Y)\\(Phospho \\(STY\\)\\)" = "p\\1",
    "(S|T|Y)\\(UniMod:21\\)" = "p\\1",
    "\\(SILAC-(R|K)-(H|L)\\)" = "",
    "\\(Oxidation \\(M\\)\\)" = "",
    "\\(Acetyl \\(Protein N-term\\)\\)" = "",
    "\\(UniMod:\\d+\\)" = ""
  ))
}

fun_modseq_mq20tojpt <- function(seq_p) {
  # Convert phospho modifications (with optional spaces) to pS, pT, pY; 
  # remove Acetyl, Oxidation, underscores at ends, and various asmod/label tags.
  str_replace_all(seq_p, c(
    "(S|T|Y)\\(Phospho\\s*\\(STY\\)\\)" = "p\\1",
    "\\(Acetyl\\s*\\(Protein\\s*N-term\\)\\)" = "",
    "\\(Oxidation\\s*\\(M\\)\\)" = "",
    "^_|_$" = "",
    "(Lys8|Arg10|Pro6|Phe6|Phe10)(?:_)?(asmod|label)" = ""
  ))
}

#---------------------------
# Input Data Path
#---------------------------
path_reporttsv_HpHLFQ <- file.path(data_dir, "report.tsv")

#---------------------------
# Data Reading & Preparation
#---------------------------
if (!"reporttsv_HpHLFQ" %in% ls()) {
  reporttsv_HpHLFQ <- fread(path_reporttsv_HpHLFQ, sep = "\t")
}

setDT(reporttsv_HpHLFQ)

# Filter for runs and create new columns for sample names and groups
HpHLFQ_fil <- reporttsv_HpHLFQ[grepl("SS2", Run)][,
                                                  c("sample_simple", "samplegroup", "jptseq") := {
                                                    # Adjust run names to simplified sample names
                                                    sample_simple <- gsub("Popeye_20230220_MVB_HStdia_SS2_", "Pro2_1_", Run)
                                                    sample_simple <- gsub("Popeye_20230228_MVB_HStdia_SS2_", "Pro2_", sample_simple)
                                                    sample_simple <- gsub("Popeye_20230309_MVB_HStdia_SS2_", "Pro2_", sample_simple)
                                                    sample_simple <- gsub("_[A-Z0-9]*_[0-9]_[0-9]*$", "", sample_simple)
                                                    sample_simple <- gsub("100ng_*", "", sample_simple)
                                                    
                                                    # Assign sample groups based on naming pattern
                                                    samplegroup <- sample_simple
                                                    samplegroup <- gsub("1$|2$|3$", "_BSA_DMSO", samplegroup)
                                                    samplegroup <- gsub("4|5|6", "_BSA_MEKi", samplegroup)
                                                    samplegroup <- gsub("7|8|9", "_GF_BSA", samplegroup)
                                                    samplegroup <- gsub("10|11|12", "_GF_MEKi", samplegroup)
                                                    
                                                    # Process modified sequences
                                                    jptseq <- fun_modseq_diannMQlib2pSTY(Modified.Sequence)
                                                    .(sample_simple, samplegroup, jptseq)
                                                  }]

HpHLFQ_fil <- as.data.frame(HpHLFQ_fil)

# Save filtered data
saveRDS(HpHLFQ_fil, file = file.path(data_dir, "CellLinePanel_Global_unfiltered.rds"), compress = "gzip")

#---------------------------
# Data Normalization
#---------------------------
HpH_m_all <- HpHLFQ_fil %>%
  ungroup() %>%
  rowwise() %>%
  filter(
    Q.Value < 0.05,
    grepl("UniMod:21", Modified.Sequence),
    PTM.Site.Confidence > 0.5,
    PTM.Q.Value < 0.05
  ) %>%
  mutate(
    log_int = log10(Ms1.Area),
    sample_simple = gsub("Pro2_", "", sample_simple)
  ) %>%
  select(Precursor.Id, sample_simple, log_int) %>%
  pivot_wider(names_from = sample_simple, values_from = log_int) %>%
  column_to_rownames("Precursor.Id") %>%
  as.matrix()

HpH_norm <- normalizeCyclicLoess(HpH_m_all, weights = NULL, span = 0.7, iterations = 3, method = "fast")

saveRDS(HpH_norm, file = file.path(data_dir, "CellLinePanel_Global_norm.rds"), compress = "gzip")

#---------------------------
# PCA Analysis - All Samples
#---------------------------
HpH_t <- t(HpH_norm[complete.cases(HpH_norm), ])
n_peps <- ncol(HpH_t)

pca <- prcomp(HpH_t, scale. = FALSE)
pca_scores <- as.data.frame(pca$x)
pca$percVarExpl <- ((pca$sdev^2) / sum(pca$sdev^2)) * 100
names(pca$percVarExpl) <- colnames(pca$rotation)

pca_scores$CellLine <- substr(rownames(HpH_t), 1, 1)
pca_scores$Number <- as.numeric(substr(rownames(HpH_t), 2, nchar(rownames(HpH_t))))
labels <- c("C" = "Caco-2", "D" = "DLD-1", "H" = "HCT116")
pca_scores$CellLineLabel <- labels[pca_scores$CellLine]

shape_mapping <- scale_shape_manual(values = c("Caco-2" = 16, "DLD-1" = 17, "HCT116" = 15))
color_mapping <- scale_color_manual(values = c("Caco-2" = "deepskyblue3", "DLD-1" = "lightgoldenrod3", "HCT116" = "azure4"))

PCA_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = CellLineLabel, shape = CellLineLabel)) +
  geom_point(size = 4) +
  labs(
    x = paste("PC1 (", round(pca$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC2 (", round(pca$percVarExpl["PC2"], 2), "%)", sep = " "),
    subtitle = "All cell-lines"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_viridis_d(option = "E", end = 0.9) +
  shape_mapping +
  color_mapping

print(n_peps)
print(PCA_all)

ggsave(file.path(plots_dir, "20240308_PCA_allCellLines.pdf"), plot = PCA_all, width = 3, height = 4, units = "in", dpi = 300)

#---------------------------
# PCA Analysis - HCT116 Only
#---------------------------
Hvec <- paste0("H", 1:12)
HpH_norm_H <- HpH_norm[, Hvec]
HpH_t_H <- t(HpH_norm_H[complete.cases(HpH_norm_H), ])
n_peps <- ncol(HpH_t_H)

pca_H <- prcomp(HpH_t_H, scale. = FALSE)
pca_scores_H <- as.data.frame(pca_H$x)
pca_H$percVarExpl <- ((pca_H$sdev^2) / sum(pca_H$sdev^2)) * 100
names(pca_H$percVarExpl) <- colnames(pca_H$rotation)

pca_scores_H$Number <- as.numeric(substr(Hvec, 2, nchar(Hvec)))
pca_scores_H$ColorCategory <- cut(
  pca_scores_H$Number,
  breaks = c(0, 3, 6, 9, 12),
  labels = c('Control', 'MEKi', 'GFmix', 'MEKi+GFmix'),
  include.lowest = TRUE
)

PCA_HCT116 <- ggplot(pca_scores_H, aes(x = PC1, y = PC2, color = ColorCategory)) +
  geom_point(size = 4, shape = 15) +
  labs(
    x = paste("PC1 (", round(pca_H$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC2 (", round(pca_H$percVarExpl["PC2"], 2), "%)", sep = " "),
    subtitle = "HCT116"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom", 
    legend.title = element_blank()
  )

print("HCT116")
print(n_peps)
print(PCA_HCT116)

#---------------------------
# PCA Analysis - Caco-2
#---------------------------
Cvec <- paste0("C", 1:12)
HpH_norm_C <- HpH_norm[, Cvec]
HpH_t_C <- t(HpH_norm_C[complete.cases(HpH_norm_C), ])
n_peps <- ncol(HpH_t_C)

pca_C <- prcomp(HpH_t_C, scale. = FALSE)
pca_scores_C <- as.data.frame(pca_C$x)
pca_C$percVarExpl <- ((pca_C$sdev^2) / sum(pca_C$sdev^2)) * 100
names(pca_C$percVarExpl) <- colnames(pca_C$rotation)

pca_scores_C$Number <- as.numeric(substr(Cvec, 2, nchar(Cvec)))
pca_scores_C$ColorCategory <- cut(
  pca_scores_C$Number,
  breaks = c(0, 3, 6, 9, 12),
  labels = c('Control', 'MEKi', 'GFmix', 'MEKi+GFmix'),
  include.lowest = TRUE
)

CACO2_PCA_C6 <- ggplot(pca_scores_C, aes(x = PC1, y = PC2, color = ColorCategory)) +
  geom_point(size = 4, shape = 16) +
  labs(
    x = paste("PC1 (", round(pca_C$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC2 (", round(pca_C$percVarExpl["PC2"], 2), "%)", sep = " "),
    subtitle = paste0("CACO-2, ", n_peps, " peptides")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_viridis_d(option = "D", end = 0.9)

print(CACO2_PCA_C6)

# Re-run Caco-2 PCA excluding C6 samples for demonstration
Cvec <- paste0("C", c(1:5,7:12))
HpH_norm_C <- HpH_norm[, Cvec]
HpH_t <- t(HpH_norm_C[complete.cases(HpH_norm_C), ])
n_peps <- ncol(HpH_t)

pca_C <- prcomp(HpH_t, scale. = FALSE)
pca_scores_C <- as.data.frame(pca_C$x)
pca_C$percVarExpl <- ((pca_C$sdev^2) / sum(pca_C$sdev^2)) * 100
names(pca_C$percVarExpl) <- colnames(pca_C$rotation)

pca_scores_C$Number <- as.numeric(substr(Cvec, 2, nchar(Cvec)))
pca_scores_C$ColorCategory <- cut(
  pca_scores_C$Number,
  breaks = c(0, 3, 6, 9, 12),
  labels = c('Control', 'MEKi', 'GFmix', 'MEKi+GFmix'),
  include.lowest = TRUE
)

CACO2_PCA <- ggplot(pca_scores_C, aes(x = PC1, y = PC2, color = ColorCategory)) +
  geom_point(size = 4, shape = 16) +
  labs(
    x = paste("PC1 (", round(pca_C$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC2 (", round(pca_C$percVarExpl["PC2"], 2), "%)", sep = " "),
    subtitle = "Caco-2"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print("CACO2")
print(n_peps)
print(CACO2_PCA)

#---------------------------
# PCA Analysis - DLD-1
#---------------------------
Dvec <- paste0("D", c(1:6,8:12))
HpH_norm_D <- HpH_norm[, grep("D", colnames(HpH_norm))]
HpH_t_D <- t(HpH_norm_D[complete.cases(HpH_norm_D), ])
n_peps <- ncol(HpH_t_D)

pca_D <- prcomp(HpH_t_D, scale. = FALSE)
pca_scores_D <- as.data.frame(pca_D$x)
pca_D$percVarExpl <- ((pca_D$sdev^2) / sum(pca_D$sdev^2)) * 100
names(pca_D$percVarExpl) <- colnames(pca_D$rotation)

pca_scores_D$Number <- as.numeric(substr(Dvec, 2, nchar(Dvec)))
pca_scores_D$ColorCategory <- cut(
  pca_scores_D$Number,
  breaks = c(0, 3, 6, 9, 12),
  labels = c('Control', 'MEKi', 'GFmix', 'MEKi+GFmix'),
  include.lowest = TRUE
)

DLD1_PCA_12 <- ggplot(pca_scores_D, aes(x = PC1, y = PC2, color = ColorCategory)) +
  geom_point(size = 4, shape = 16) +
  geom_text_repel(aes(label = rownames(pca_scores_D)), hjust = 0, vjust = 0, size = 4) +
  labs(
    x = paste("PC1 (", round(pca_D$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC2 (", round(pca_D$percVarExpl["PC2"], 2), "%)", sep = " "),
    subtitle = "DLD-1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_viridis_d(option = "D", end = 0.9)

print(DLD1_PCA_12)

DLD1_PCA_14 <- ggplot(pca_scores_D, aes(x = PC1, y = PC4, color = ColorCategory)) +
  geom_point(size = 4, shape = 17) +
  labs(
    x = paste("PC1 (", round(pca_D$percVarExpl["PC1"], 2), "%)", sep = " "),
    y = paste("PC4 (", round(pca_D$percVarExpl["PC4"], 2), "%)", sep = " "),
    subtitle = "DLD-1"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  guides(color = FALSE, size = FALSE)

print("DLD1")
print(n_peps)
print(DLD1_PCA_14)

#---------------------------
# Combine Selected PCA Plots
#---------------------------
final_plot <- PCA_HCT116 + DLD1_PCA_14 + CACO2_PCA +
  plot_layout(guides = "collect", ncol = 3) &
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, "20240308_PCA.pdf"), plot = final_plot, width = 9, height = 3, units = "in", dpi = 300)

###############################################################################
# End of Script
###############################################################################