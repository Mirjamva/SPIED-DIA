#===============================================================================
# Title:       Cell Line Panel Global data Differential Expression Analysis
# Author:      Mirjam van Bentum
# Date:        2024-12-16
# 
# Description: 
#   This script performs differential expression analysis on normalized and
#   filtered global phosphoproteomic data for multiple cell lines (HCT116, DLD-1, Caco2).
#   The analysis includes:
#     - Limma-based differential expression analysis for each cell line
#     - Hierarchical clustering of significantly regulated phospho-peptides
#     - Visualization through heatmaps and volcano plots
#     - Identification and visualization of regulation per cluster
#     - Manual selection of clusters for synergistic regulation
#     - Generation of heatmaps for selected clusters
#     - Kinase library enrichment analysis to identify enriched kinase signatures
#       based on PTM-SEA and ikip_db databases
#       (requires Cluster_Kinase_Enrichment.R script)
# 
# Inputs:
#   - Normalized and filtered data from global analysis (produced with cellLinePanel_Global_Norm)
#   - Annotation file: HpH_annotation_20230626.csv
#   - Cluster_Kinase_Enrichment.R script (contains functions for kinase enrichment analysis)
#   - Kinase sequence pairs: Cluster_Kinase_sequence_pairs.rds
#   - Flanking sequences annotation: flankingseqsanno
#
# Outputs:
#   - Differential expression results saved as RDS files for each cell line
#   - Heatmaps and cluster zoom plots saved as PDF files for each cell line
#   - Kinase library input tables for enrichment analysis with Cantley et al. kinase signature enrichment tool
#   - Tables compatible for analysis with Cantley et al. kinase signature enrichment
#
# Dependencies:
#   - R packages: tidyverse, data.table, patchwork, limma, circlize, ComplexHeatmap, scales, purrr, stringr
#   - External scripts: Cluster_Kinase_Enrichment.R
#   - Databases: ikip_db, PTM-SEA
#
# Execution Steps:
#   1. Load necessary libraries and functions
#   2. Define file paths
#   3. Define custom functions for sequence modifications and data export
#   4. Load and transform data
#   Per cell line (HCT116, DLD-1, Caco2)
#     5. Perform differential expression analysis for each cell line 
#     6. Filter significant peptides based on F p-value threshold
#     7. Prepare row annotations with log fold changes and adjusted p-values
#     8. Scale data for heatmap visualization
#     9. Define column annotations and color mappings
#    10. Generate and save heatmaps for each cell line
#    11. Perform hierarchical clustering and visualize cluster-specific regulations
#    12. Conduct kinase library enrichment analysis on selected clusters
#    13. Export kinase library data for further enrichment analysis with cantley et al.
#    
#===============================================================================



#### libraries ####
library(tidyverse)
library(data.table)
library(patchwork)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(scales)
library(purrr)

#===============================================================================
# File Paths
#===============================================================================
data_dir <- "Data"
script_dir <- setwd()
output_dir <- "Plots"

#===============================================================================
# Function Definitions
#===============================================================================

# Function to modify sequence for DIANN MQ library
fun_modseq_diannMQlib2pSTY <- function(seq_p) {
  patterns <- c(
    "_" = "",
    "(S|T|Y)\\(Phospho \\(STY\\)\\)" = "p\\1",
    "(S|T|Y)\\(UniMod:21\\)" = "p\\1",
    "\\(SILAC-[RK]-[HL]\\)" = "",
    "\\(Oxidation ?\\(M\\)\\)" = "",
    "\\(Acetyl ?\\(Protein ?N-term\\)\\)" = "",
    "\\(Acetyl ?\\(N-term\\)\\)" = "",
    "\\(UniMod:\\d+\\)" = ""
  )
  str_replace_all(seq_p, patterns) %>% 
    str_replace_all("\\s+", "")
}

# Function to modify sequence for MQ 20 to JPT conversion
fun_modseq_mq20tojpt <- function(seq_p) {
  patterns <- c(
    "(S|T|Y)\\(Phospho ?\\(STY\\)\\)" = "p\\1",
    "\\(Acetyl ?\\(Protein ?N-term\\)\\)" = "",
    "\\(Acetyl ?\\(N-term\\)\\)" = "",
    "\\(Oxidation ?\\(M\\)\\)" = "",
    "\\(Lys8_asmod\\)|\\(Arg10_asmod\\)|\\(Pro6_asmod\\)|\\(Phe6_asmod\\)|\\(Lys8_label\\)|\\(Phe10_asmod\\)" = "",
    "\\(Lys8asmod\\)|\\(Arg10asmod\\)|\\(Pro6asmod\\)|\\(Phe6asmod\\)|\\(Lys8label\\)|\\(Phe10asmod\\)" = "",
    "^_|_$" = ""
  )
  str_replace_all(seq_p, patterns)
}

fun_modseq_kinlib <- function(sequence) {
  sequence %>%
    # Replace S, Y, T with their phosphorylated forms
    str_replace_all("([SYT])(?:\\(UniMod:21\\)|\\(ph\\))", "\\1*") %>%
    # Remove unwanted modifications
    str_remove_all(
      "\\(UniMod:\\d+\\)|" %>% 
        paste0(
          "\\(Acetyl \\(Protein N-term\\)\\)|",
          "\\(Acetyl \\(N-term\\)\\)|",
          "\\(SILAC-[RK]-L\\)|",
          "\\(M\\(Oxidation \\(M\\)\\)\\)|",
          "\\(Oxidation \\(M\\)\\)|",
          "\\)[2-4]*$|",
          "^[^\\(]*\\(|",
          "\\(ox\\)|",
          "\\(ac\\)"
        )
    ) %>%
    # Remove trailing numbers
    str_remove_all("[1-4]*$")
}

# Function to export kinase library data
export_kinlib_data <- function(data, coefficients, cell_line, output_dir) {
  today_date <- Sys.Date()
  coefficient_name <- paste(coefficients, collapse = "_")
  file_name <- sprintf("%s_%s_kinlib_input_%s.txt", today_date, cell_line, coefficient_name)
  full_path <- file.path(output_dir, file_name)
  
  data_selected <- data %>%
    select(kinlibseq, all_of(coefficients)) %>%
    filter(str_count(kinlibseq, "\\*") < 2, !is.na(.data[[coefficients[1]]]))
  
  write.table(data_selected, file = full_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("File saved to: ", full_path)
}

# Function to retrieve top table for a specific coefficient
get_top_table <- function(fit, coef, number, suffix) {
  topTable(fit, coef = coef, adjust = "BH", number = number) %>%
    rename_with(~ paste0(., ".", suffix)) %>%
    rownames_to_column("Modified.Sequence")
}

#===============================================================================
# Source External Scripts
#===============================================================================
source(file.path(script_dir, "Cluster_Kinase_Enrichment.R"))

#===============================================================================
# Load and Transform Data
#===============================================================================
HpH_norm <- readRDS(file.path(data_dir, "CellLinePanel_Global_norm.rds"))

tomatch <- read.delim(file.path(data_dir, "HpH_annotation_20230626.csv"), 
                      header = TRUE, sep = ",", 
                      stringsAsFactors = FALSE, 
                      dec = " ", 
                      check.names = FALSE) %>% 
  mutate(
    `Modified sequence` = gsub("_", "", `Modified sequence`),
    c_psite_id = gsub("_", "", c_psite_id),
    plot_name = paste0(geneid, "_", c_psite_id, "(", `Modified sequence`, ")") %>%
      str_replace_all(c(
        "S\\(Phospho \\(STY\\)\\)" = "(ph)S",
        "Y\\(Phospho \\(STY\\)\\)" = "(ph)Y",
        "T\\(Phospho \\(STY\\)\\)" = "(ph)T",
        "\\(Acetyl \\(Protein N-term\\)\\)" = "(ac)",
        "\\(Acetyl \\(N-term\\)\\)" = "(ac)",
        "M\\(Oxidation \\(M\\)\\)" = "(ox)M"
      )),
    jptseq = fun_modseq_diannMQlib2pSTY(`Modified sequence`),
    plot_jpt_seq = paste0(geneid, "_", c_psite_id, "(", jptseq, ")"),
    plot_genepsite = paste0(geneid, "_", c_psite_id)
  ) %>%
  select(jptseq, plot_genepsite) %>%
  distinct()

#===============================================================================
# Differential Expression Analysis for HCT116 ####
#===============================================================================

# Define column order and subset normalized data
col_order <- paste0("H", 1:12)
norm_limma <- HpH_norm[, col_order]

# Define factors for design matrix as per original Caco-2 specification
Ligandmix <- factor(
  c("BSA", "BSA", "BSA", "BSA", "BSA", "BSA",
    "LigandMix", "LigandMix", "LigandMix",
    "LigandMix", "LigandMix", "LigandMix"),
  levels = c("BSA", "LigandMix")
)

MEKi <- factor(
  c("DMSO", "DMSO", "DMSO", "MEKi", "MEKi", "MEKi",
    "DMSO", "DMSO", "DMSO", "MEKi", "MEKi", "MEKi"),
  levels = c("DMSO", "MEKi")
)

# Create treatment combinations
TS <- paste(Ligandmix, MEKi, sep = ".")
TS <- factor(TS, levels = unique(TS))

# Create design matrix
design <- model.matrix(~0 + TS)
colnames(design) <- levels(TS)

# Fit linear model and apply contrasts
fit <- lmFit(norm_limma, design)
contrast.matrix <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi = LigandMix.MEKi - BSA.MEKi,
  SynSign = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2 = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design
)
fit <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)

# Extract results
results <- decideTests(fit2)
total_ids <- colSums(summary(results))

# Retrieve top tables for all contrasts
LMvsBSAinDMSO <- get_top_table(fit2, 1, total_ids[1], "LMvsBSAinDMSO")
LMvsBSAinMEKi  <- get_top_table(fit2, 2, total_ids[2], "LMvsBSAinMEKi")
SynSign        <- get_top_table(fit2, 3, total_ids[3], "SynSign")
MEKivsDMSOinBSA <- get_top_table(fit2, 4, total_ids[4], "MEKivsDMSOinBSA")
MEKivsDMSOinLM <- get_top_table(fit2, 5, total_ids[5], "MEKivsDMSOinLM")
SynSign2       <- get_top_table(fit2, 6, total_ids[6], "SynSign2")

# Combine differential expression results
df_HCT116 <- norm_limma %>%
  as_tibble() %>%
  add_column(F_pval = fit2$F.p.value, Modified.sequence = rownames(fit2$p.value)) %>%
  mutate(across(everything(), ~ replace(., is.infinite(.), NA))) %>%
  filter(rowSums(is.na(.)) < 4) %>%
  select(Modified.sequence, F_pval, starts_with("H")) %>%
  left_join(SynSign, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinDMSO, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinMEKi, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinBSA, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinLM, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(SynSign2, by = c("Modified.sequence" = "Modified.Sequence"))

# Save differential expression results
saveRDS(df_HCT116, file.path(data_dir, "20241216_CellLinePanel_Global_HCT116_DiffExpr.rds"))

#===============================================================================
# Filter Significant Peptides
#===============================================================================
full_sign <- df_HCT116 %>%
  filter(F_pval < 0.1) %>%
  rename(Modified.Sequence = Modified.sequence)

#===============================================================================
# Prepare Row Annotations
#===============================================================================
# Prepare row annotations directly from full_sign using already merged data
row_annotation <- full_sign %>%
  # Adjust p-values using Benjamini-Hochberg method and apply -log10 transformation
  mutate(
    `Interaction - Pval` = -log10(p.adjust(P.Value.SynSign, method = "BH")),
    `GFmix w/o MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinDMSO, method = "BH")),
    `GFmix w MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinMEKi, method = "BH")),
    `MEKi w/o GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinBSA, method = "BH")),
    `MEKi w GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinLM, method = "BH"))
  ) %>%
  # Rename logFC columns to more descriptive names
  rename(
    `Interaction - logFC` = logFC.SynSign,
    `GFmix w/o MEKi - logFC` = logFC.LMvsBSAinDMSO,
    `GFmix w MEKi - logFC` = logFC.LMvsBSAinMEKi,
    `MEKi w/o GFmix - logFC` = logFC.MEKivsDMSOinBSA,
    `MEKi w GFmix - logFC` = logFC.MEKivsDMSOinLM
  ) %>%
  # Set Modified.Sequence as row names
  column_to_rownames("Modified.Sequence") %>%
  # Select only the relevant annotation columns
  select(
    `Interaction - logFC`, `Interaction - Pval`,
    `GFmix w/o MEKi - logFC`, `GFmix w/o MEKi - Pval`,
    `GFmix w MEKi - logFC`, `GFmix w MEKi - Pval`,
    `MEKi w/o GFmix - logFC`, `MEKi w/o GFmix - Pval`,
    `MEKi w GFmix - logFC`, `MEKi w GFmix - Pval`
  )

row_annotation_HCT116_4f <- row_annotation

# check pvals manually
HCT116_adjpvals <- cbind(row_annotation, full_sign)
modseq_tbl <- gsub("\\(SILAC-[R|K]-L\\)|[0-9]*$", "", HCT116_adjpvals$Modified.Sequence)
jptseq_tbl <- fun_modseq_diannMQlib2pSTY(modseq_tbl)
precursor_charge <- str_extract(HCT116_adjpvals$Modified.Sequence, "[0-9]*$")
tomatch %>% 
  filter(jptseq %in% jptseq_tbl) %>%
  group_by(jptseq, plot_genepsite) %>%
  slice(1) %>% ungroup() -> df
df[match(jptseq_tbl, df$jptseq),] %>%
  pull(plot_genepsite) -> rownames_heatmap
HCT116_adjpvals <- cbind(HCT116_adjpvals, HGNCids = rownames_heatmap) 

HCT116_adjpvals %>% 
  filter(grepl("JUN_S63", HGNCids)) %>% glimpse
HCT116_adjpvals %>% 
  filter(grepl("SPAG9", HGNCids)) %>% glimpse
HCT116_adjpvals %>% 
  filter(grepl("MAP3K3", HGNCids)) %>% glimpse

#===============================================================================
# Scale Data for Heatmap
#===============================================================================
z.mat <- full_sign %>%
  select(starts_with("H")) %>%
  as.matrix() %>% t() %>% 
  scale(center = TRUE, scale = TRUE) %>%
  t()

#===============================================================================
# Column Annotations
#===============================================================================
column_annotation <- data.frame(
  GFmix = rep(c("BSA", "GFmix"), each = 6),
  MEKi = c(rep("DMSO", 3), rep("MEKi", 3), rep("DMSO", 3), rep("MEKi", 3))
)
rownames(column_annotation) <- col_order

col_annotate_mapping <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green3", MEKi = "purple")
)

#===============================================================================
# Color Mappings
#===============================================================================
color_mapping_row <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w/o MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w/o GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(c(min(row_annotation$`Interaction - Pval`, na.rm = TRUE), 3), c("white", "black")),
  `GFmix w/o MEKi - Pval` = colorRamp2(c(min(row_annotation$`GFmix w/o MEKi - Pval`, na.rm = TRUE), 3), c("white", "black")),
  `GFmix w MEKi - Pval` = colorRamp2(c(min(row_annotation$`GFmix w MEKi - Pval`, na.rm = TRUE), 3), c("white", "black")),
  `MEKi w/o GFmix - Pval` = colorRamp2(c(min(row_annotation$`MEKi w/o GFmix - Pval`, na.rm = TRUE), 3), c("white", "black")),
  `MEKi w GFmix - Pval` = colorRamp2(c(min(row_annotation$`MEKi w GFmix - Pval`, na.rm = TRUE), 3), c("white", "black"))
)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))

#===============================================================================
# Create and Save Heatmap
#===============================================================================
set.seed(123)

heatmap_H_LFQ_int_up <- ComplexHeatmap::Heatmap(
  z.mat,
  name = "Abundance",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(df = column_annotation, col = col_annotate_mapping),
  right_annotation = HeatmapAnnotation(
    df = row_annotation, 
    which = "row", 
    col = color_mapping_row,
    annotation_legend_param = list(
      `Interaction - logFC` = list(title = "log2(FC)"),
      `GFmix w MEKi - Pval` = list(title = "-log10(pval)")
    ),
    show_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  ),
  row_names_gp = gpar(fontsize = 1),
  cluster_rows = T,
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  row_split = 14,
  col = col_fun,
  column_order = col_order
)

pdf(file.path(output_dir, "20241216_heat_HCT116_Global_all.pdf"), width = 7, height = 14)
draw(heatmap_H_LFQ_int_up)
clustrows <- row_order(heatmap_H_LFQ_int_up)
dev.off()

#===============================================================================
# Cluster Zoom and Enrichment
#===============================================================================

#first way
df <- clustrows %>% 
  enframe(name = "sublist_index") %>%
  mutate(cluster = as.integer(sublist_index)) %>%
  unnest_longer(value, indices_include = TRUE) %>%
  rename(element_in_cluster = value_id, value = value) %>%
  arrange(value)

z.mat %>%
  data.frame() %>%
  cbind(df) %>%
  mutate(Modified.sequence = full_sign$Modified.Sequence) -> plotcluster_mean

# second way, doesn't work yet
plotcluster_mean <- z.mat %>%
  as_tibble() %>%
  bind_cols(
    cluster = df$cluster,
    Modified.Sequence = full_sign$Modified.Sequence
  )

saveRDS(plotcluster_mean, file.path(data_dir, "20241216_Global_HCT116_DiffExpr_clustered.rds"))

df_long <- plotcluster_mean %>% 
  pivot_longer(cols = H1:H12, names_to = "variable", values_to = "value") %>%
  mutate(
    category = if_else(variable %in% paste0("H", 1:6), "BSA", "GFmix"),
    treatment = if_else(variable %in% paste0("H", c(1:3,7:9)), "DMSO", "MEKi")
  )

#===============================================================================
# Calculate Statistics
#===============================================================================
stats <- df_long %>%
  group_by(cluster, category, treatment) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = 'drop'
  )

#===============================================================================
# Plot Cluster Means
#===============================================================================
meanclusterplots <- ggplot(stats, aes(x = category, y = mean, color = treatment, group = treatment)) +
  geom_point(position = position_dodge(0.1), size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(y = "Mean Value", x = "Category", color = "Treatment") +
  facet_grid(rows = vars(cluster)) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "vertical"
  ) +
  scale_color_manual(values = c("DMSO" = "green3", "MEKi" = "purple"))

ggsave(file.path(output_dir, "20241216_HCT116_clustmean.pdf"), 
       meanclusterplots, height = 12, width = 1.5)

#===============================================================================
# Kinase Library Enrichment
#===============================================================================
kinase_sequence_pairs <- readRDS(file.path(data_dir, "Cluster_Kinase_sequence_pairs.rds"))
AKT_seqs <- kinase_sequence_pairs %>% filter(str_detect(kinase, "AKT"))
JNK_seqs <- kinase_sequence_pairs %>% filter(str_detect(kinase, "JNK"))

cluster_HCT116_flanking <- plotcluster_mean %>%
  mutate(
    kinlibseq = fun_PrecurIDto_sty(Modified.Sequence),
    kinlibseq = map(kinlibseq, split_sequences)
  ) %>%
  unnest(kinlibseq) %>%
  left_join(flankingseqsanno, by = c("kinlibseq" = "phosphopeptide")) %>%
  mutate(
    flanking = if_else(kinlibseq == "AASEVAGVVANAPsPPESSSLCASK", "GVVANAPSPPESSSL", flanking)
  )

plot_clusterenrichement(c(6), "HCT116", cluster_HCT116_flanking)

View(readRDS(file.path(data_dir, paste0("HCT116_6.rds")))%>%
       mutate(ER_notlog = Freq_interesting / Freq_background))

#===============================================================================
# Zoomed Heatmap
#===============================================================================
HCT116_clust6 <- plotcluster_mean %>%
  filter(cluster == 6)

modseq_tbl <- gsub("\\(SILAC-[R|K]-L\\)|[0-9]*$", "", HCT116_clust6$Modified.Sequence)
jptseq_tbl <- fun_modseq_diannMQlib2pSTY(modseq_tbl)
precursor_charge <- str_extract(HCT116_clust6$Modified.Sequence, "[0-9]*$")
tomatch %>% 
  filter(jptseq %in% jptseq_tbl) %>%
  group_by(jptseq, plot_genepsite) %>%
  slice(1) %>% ungroup() -> df
df[match(jptseq_tbl, df$jptseq),] %>%
  pull(plot_genepsite) -> rownames_heatmap

z.mat_clustzoom <- HCT116_clust6 %>%
  dplyr::select(starts_with("H")) %>%
  as.matrix()
rownames(z.mat_clustzoom) <- rownames_heatmap

row_annotation_cluster6 <- row_annotation %>%
  rownames_to_column("Modified.Sequence") %>%
  filter(`Modified.Sequence` %in% HCT116_clust6$Modified.Sequence) %>%
  select(`Interaction - logFC`, `Interaction - Pval`)

color_mapping_row_cluster6 <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(c(min(row_annotation$`Interaction - Pval`, na.rm = TRUE), 3), c("white", "black"))
)

heatmap_H_LFQ_clustzoom <- Heatmap(
  z.mat_clustzoom,
  name = "HCT116 clust 6",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(df = column_annotation, col = col_annotate_mapping, show_legend = FALSE),
  right_annotation = HeatmapAnnotation(
    df = row_annotation_cluster6, 
    which = "row", 
    col = color_mapping_row_cluster6,
    show_legend = FALSE
  ),
  row_names_gp = gpar(fontsize = 5),
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  col = col_fun,
  column_order = col_order, 
  show_heatmap_legend = FALSE
)

pdf(file.path(output_dir, "heat_Global_HCT116_clust6.pdf"), width = 3.5, height = 10)
draw(heatmap_H_LFQ_clustzoom)
dev.off() 

cbind(row_annotation_cluster6, rownames_heatmap) %>% View

#===============================================================================
# Prepare Kinase Library Input
#===============================================================================
KINASElib_input <- norm_limma %>%
  as_tibble(rownames = "Modified.sequence") %>%
  # Expand pval and FC matrices into separate columns
  bind_cols(
    as_tibble(fit2$p.value, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("pval.", .))
  ) %>%
  bind_cols(
    as_tibble(fit2$coefficients, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("FC.", .))
  ) %>%
  ungroup() %>%
  # Clean the Modified.sequence to create kinlibseq
  mutate(
    kinlibseq = fun_modseq_kinlib(Modified.sequence)
  )

#===============================================================================
# Export Kinase Library Data
#===============================================================================
# Proceed with selecting and exporting the data
KINASElib_export <- KINASElib_input %>%
  select(kinlibseq, starts_with("FC."), starts_with("pval.")) %>%
  as.matrix()

# Define your coefficients list as before
coefficients_list <- list(
  c("FC.SynSign", "pval.SynSign"),
  c("FC.SynSign2", "pval.SynSign2"),
  c("FC.LMvsBSAinDMSO", "pval.LMvsBSAinDMSO"),
  c("FC.LMvsBSAinMEKi", "pval.LMvsBSAinMEKi"),
  c("FC.MEKivsDMSOinBSA", "pval.MEKivsDMSOinBSA"),
  c("FC.MEKivsDMSOinLM", "pval.MEKivsDMSOinLM")
)

# Export Kinase Library Data
walk(coefficients_list, ~ export_kinlib_data(KINASElib_input, .x, "HCT116", data_dir))

#### Caco-2 Differential Expression Analysis ####

#===============================================================================
# Differential Expression Analysis for Caco-2 ####
#===============================================================================

# Define column order specific to Caco-2 (excluding C6 if present)
col_order_caco2 <- paste0("C", c(1:5, 7:12))
norm_limma_caco2 <- HpH_norm[, col_order_caco2]

# Define factors for design matrix as per original Caco-2 specification
Ligandmix_caco2 <- factor(
  c("BSA", "BSA", "BSA", "BSA", "BSA",
    "LigandMix", "LigandMix", "LigandMix",
    "LigandMix", "LigandMix", "LigandMix"),
  levels = c("BSA", "LigandMix")
)

MEKi_caco2 <- factor(
  c("DMSO", "DMSO", "DMSO", "MEKi", "MEKi",
    "DMSO", "DMSO", "DMSO", "MEKi", "MEKi", "MEKi"),
  levels = c("DMSO", "MEKi")
)

# Create treatment combinations
TS_caco2 <- paste(Ligandmix_caco2, MEKi_caco2, sep = ".")
TS_caco2 <- factor(TS_caco2, levels = unique(TS_caco2))

# Create design matrix
design_caco2 <- model.matrix(~0 + TS_caco2)
colnames(design_caco2) <- levels(TS_caco2)

# Fit linear model and apply contrasts
fit_caco2 <- lmFit(norm_limma_caco2, design_caco2)
contrast.matrix_caco2 <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi  = LigandMix.MEKi - BSA.MEKi,
  SynSign        = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2       = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design_caco2
)
fit_caco2 <- contrasts.fit(fit_caco2, contrast.matrix_caco2)
fit2_caco2 <- eBayes(fit_caco2)

# Extract results
results_caco2 <- decideTests(fit2_caco2)
total_ids_caco2 <- colSums(summary(results_caco2))

# Retrieve top tables for all contrasts using the predefined function
LMvsBSAinDMSO_caco2 <- get_top_table(fit2_caco2, 1, total_ids_caco2[1], "LMvsBSAinDMSO")
LMvsBSAinMEKi_caco2  <- get_top_table(fit2_caco2, 2, total_ids_caco2[2], "LMvsBSAinMEKi")
SynSign_caco2        <- get_top_table(fit2_caco2, 3, total_ids_caco2[3], "SynSign")
MEKivsDMSOinBSA_caco2 <- get_top_table(fit2_caco2, 4, total_ids_caco2[4], "MEKivsDMSOinBSA")
MEKivsDMSOinLM_caco2 <- get_top_table(fit2_caco2, 5, total_ids_caco2[5], "MEKivsDMSOinLM")
SynSign2_caco2       <- get_top_table(fit2_caco2, 6, total_ids_caco2[6], "SynSign2")

# Combine differential expression results
df_Caco2 <- norm_limma_caco2 %>%
  as_tibble() %>%
  add_column(F_pval = fit2_caco2$F.p.value, Modified.sequence = rownames(fit2_caco2$p.value)) %>%
  mutate(across(everything(), ~ replace(., is.infinite(.), NA))) %>%
  filter(rowSums(is.na(.)) < 4) %>%
  select(Modified.sequence, F_pval, starts_with("C")) %>%
  left_join(SynSign_caco2, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinDMSO_caco2, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinMEKi_caco2, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinBSA_caco2, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinLM_caco2, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(SynSign2_caco2, by = c("Modified.sequence" = "Modified.Sequence"))

# Save differential expression results
saveRDS(
  df_Caco2,
  file.path(data_dir, "20241216_CellLinePanel_Global_Caco2_DiffExpr.rds")
)
#===============================================================================
# Filter Significant Peptides for Caco-2 ####
#===============================================================================

full_sign_caco2 <- df_Caco2 %>%
  filter(F_pval < 0.1) %>%
  rename(Modified.Sequence = Modified.sequence)

#===============================================================================
# Prepare Row Annotations for Caco-2 ####
#===============================================================================

row_annotation_Caco2 <- full_sign_caco2 %>%
  # Adjust p-values using Benjamini-Hochberg method and apply -log10 transformation
  mutate(
    `Interaction - Pval` = -log10(p.adjust(P.Value.SynSign, method = "BH")),
    `GFmix w/o MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinDMSO, method = "BH")),
    `GFmix w MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinMEKi, method = "BH")),
    `MEKi w/o GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinBSA, method = "BH")),
    `MEKi w GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinLM, method = "BH"))
  ) %>%
  # Rename logFC columns to more descriptive names
  rename(
    `Interaction - logFC` = logFC.SynSign,
    `GFmix w/o MEKi - logFC` = logFC.LMvsBSAinDMSO,
    `GFmix w MEKi - logFC` = logFC.LMvsBSAinMEKi,
    `MEKi w/o GFmix - logFC` = logFC.MEKivsDMSOinBSA,
    `MEKi w GFmix - logFC` = logFC.MEKivsDMSOinLM
  ) %>%
  # Set Modified.Sequence as row names
  column_to_rownames("Modified.Sequence") %>%
  # Select only the relevant annotation columns
  select(
    `Interaction - logFC`, `Interaction - Pval`,
    `GFmix w/o MEKi - logFC`, `GFmix w/o MEKi - Pval`,
    `GFmix w MEKi - logFC`, `GFmix w MEKi - Pval`,
    `MEKi w/o GFmix - logFC`, `MEKi w/o GFmix - Pval`,
    `MEKi w GFmix - logFC`, `MEKi w GFmix - Pval`
  )

row_annotation_Caco2_4f <- row_annotation_Caco2

#===============================================================================
# Scale Data for Heatmap for Caco-2 ####
#===============================================================================

z.mat_caco2 <- full_sign_caco2 %>%
  select(starts_with("C")) %>%
  as.matrix() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  t()

#===============================================================================
# Column Annotations for Caco-2 ####
#===============================================================================

column_annotation_caco2 <- data.frame(
  GFmix = c(rep("BSA", 5), rep("GFmix", 6)),
  MEKi = c(rep("DMSO", 3), rep("MEKi", 2), rep("DMSO", 3), rep("MEKi", 3))
)
rownames(column_annotation_caco2) <- col_order_caco2

col_annotate_mapping_caco2 <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green3", MEKi = "purple")
)

#===============================================================================
# Color Mappings for Caco-2 ####
#===============================================================================

color_mapping_row_caco2 <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w/o MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w/o GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`Interaction - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `GFmix w/o MEKi - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`GFmix w/o MEKi - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `GFmix w MEKi - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`GFmix w MEKi - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `MEKi w/o GFmix - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`MEKi w/o GFmix - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `MEKi w GFmix - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`MEKi w GFmix - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  )
)

col_fun_caco2 <- colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))

#===============================================================================
# Create and Save Heatmap for Caco-2 ####
#===============================================================================

set.seed(123)

heatmap_Caco2_LFQ_all <- ComplexHeatmap::Heatmap(
  z.mat_caco2,
  name = "Abundance",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(
    df = column_annotation_caco2,
    col = col_annotate_mapping_caco2
  ),
  right_annotation = HeatmapAnnotation(
    df = row_annotation_Caco2,
    which = "row",
    col = color_mapping_row_caco2,
    annotation_legend_param = list(
      `Interaction - logFC` = list(title = "log2(FC)"),
      `GFmix w/o MEKi - Pval` = list(title = "-log10(pval)")
    ),
    show_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  ),
  row_names_gp = gpar(fontsize = 1),
  cluster_rows = TRUE,
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  row_split = 8, # Adjust based on clustering results
  col = col_fun_caco2,
  column_order = col_order_caco2
)

# Save the heatmap to a PDF
pdf(
  file.path(output_dir, "20241216_heat_Caco2_Global_all.pdf"),
  width = 7,
  height = 14
)
draw(heatmap_Caco2_LFQ_all)
clustrows_caco2 <- row_order(heatmap_Caco2_LFQ_all)
dev.off()

#===============================================================================
# Cluster Zoom and Enrichment for Caco-2 ####
#===============================================================================

# Prepare cluster information
df_clusters_caco2 <- clustrows_caco2 %>% 
  enframe(name = "sublist_index") %>%
  mutate(cluster = as.integer(sublist_index)) %>%
  unnest_longer(value, indices_include = TRUE) %>%
  rename(element_in_cluster = value_id, value = value) %>%
  arrange(value)

# Combine with scaled data
plotcluster_mean_caco2 <- z.mat_caco2 %>%
  as_tibble() %>%
  bind_cols(df_clusters_caco2) %>%
  mutate(Modified.sequence = full_sign_caco2$Modified.Sequence)

# Save clustered differential expression data
saveRDS(
  plotcluster_mean_caco2,
  file.path(data_dir, "20241216_Global_Caco2_DiffExpr_clustered.rds")
)

# Convert to long format for plotting
df_long_caco2 <- plotcluster_mean_caco2 %>% 
  pivot_longer(cols = C1:C12, names_to = "variable", values_to = "intensity") %>%
  mutate(
    category = case_when(
      variable %in% paste0("C", 1:5) ~ "BSA",
      variable %in% paste0("C", 7:12) ~ "GFmix"
    ),
    treatment = case_when(
      variable %in% paste0("C", c(1:3,7:9)) ~ "DMSO",
      variable %in% paste0("C", c(4:5, 10:12)) ~ "MEKi"
    )
  )

#===============================================================================
# Calculate Statistics for Caco-2 ####
#===============================================================================

stats_caco2 <- df_long_caco2 %>%
  group_by(cluster, category, treatment) %>%
  summarize(
    mean = mean(intensity, na.rm = TRUE),
    sd = sd(intensity, na.rm = TRUE),
    .groups = 'drop'
  )

#===============================================================================
# Plot Cluster Means for Caco-2 ####
#===============================================================================

meanclusterplots_caco2 <- ggplot(stats_caco2, aes(x = category, y = mean, color = treatment, group = treatment)) +
  geom_point(position = position_dodge(0.1), size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(y = "Mean Value", x = "Category", color = "Treatment") +
  facet_grid(rows = vars(cluster)) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "vertical"
  ) +
  scale_color_manual(values = c("DMSO" = "green3", "MEKi" = "purple"))

# Save the cluster means plot
ggsave(
  file.path(output_dir, "20241216_Caco2_clustmean.pdf"), 
  meanclusterplots_caco2, 
  height = 12, 
  width = 1.5
)

#===============================================================================
# Kinase Library Enrichment for Caco-2 ####
#===============================================================================

# Prepare flanking sequences for enrichment
cluster_Caco2_flanking <- plotcluster_mean_caco2 %>%
  mutate(
    kinlibseq = fun_PrecurIDto_sty(Modified.sequence),
    kinlibseq = map(kinlibseq, split_sequences)
  ) %>%
  unnest(kinlibseq) %>%
  left_join(flankingseqsanno, by = c("kinlibseq" = "phosphopeptide")) %>%
  mutate(
    flanking = if_else(kinlibseq == "AASEVAGVVANAPsPPESSSLCASK", "GVVANAPSPPESSSL", flanking)
  )

# Perform cluster enrichment (adjust cluster numbers as needed)
plot_clusterenrichement(c(1), "Caco-2", cluster_Caco2_flanking)
View(readRDS(file.path(data_dir, paste0("Caco-2_1.rds")))%>%
       mutate(ER_notlog = Freq_interesting / Freq_background))

#===============================================================================
# Zoomed Heatmap for Specific Cluster in Caco-2 ####
#===============================================================================

# Example: Zoom into cluster 1
Caco2_clust1 <- plotcluster_mean_caco2 %>%
  filter(cluster == 1)

# Prepare row names for heatmap
modseq_tbl_caco2 <- gsub("\\(SILAC-[R|K]-L\\)|[0-9]*$", "", Caco2_clust1$Modified.sequence)
jptseq_tbl_caco2 <- fun_modseq_diannMQlib2pSTY(modseq_tbl_caco2)
precursor_charge <- str_extract(Caco2_clust1$Modified.sequence, "[0-9]*$")
tomatch %>% 
  filter(jptseq %in% jptseq_tbl_caco2) %>%
  group_by(jptseq, plot_genepsite) %>%
  slice(1) %>% ungroup() -> df
df[match(jptseq_tbl_caco2, df$jptseq),] %>%
  pull(plot_genepsite) -> rownames_heatmap_caco2

z.mat_clustzoom_caco2 <- Caco2_clust1 %>%
  dplyr::select(starts_with("C", ignore.case=F)) %>%
  as.matrix()
rownames(z.mat_clustzoom_caco2) <- rownames_heatmap_caco2

# Prepare row annotations for the cluster
row_annotation_cluster_caco2 <- row_annotation_Caco2 %>%
  rownames_to_column("Modified.Sequence") %>%
  filter(`Modified.Sequence` %in% Caco2_clust1$Modified.sequence) %>%
  select(`Interaction - logFC`, `Interaction - Pval`)

color_mapping_row_cluster_caco2 <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(
    c(min(row_annotation_Caco2$`Interaction - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  )
)

# Create the zoomed heatmap
heatmap_Caco2_clustzoom <- ComplexHeatmap::Heatmap(
  z.mat_clustzoom_caco2,
  name = "Caco2 Cluster 1",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(
    df = column_annotation_caco2,
    col = col_annotate_mapping_caco2,
    show_legend = FALSE
  ),
  right_annotation = HeatmapAnnotation(
    df = row_annotation_cluster_caco2,
    which = "row",
    col = color_mapping_row_cluster_caco2,
    show_legend = FALSE
  ),
  row_names_gp = gpar(fontsize = 5),
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  col = col_fun_caco2,
  column_order = col_order_caco2,
  show_heatmap_legend = FALSE
)

# Save the zoomed heatmap to a PDF
pdf(
  file.path(output_dir, "heat_Global_Caco2_clust1.pdf"), 
  width = 5, 
  height = 8
)
draw(heatmap_Caco2_clustzoom)
dev.off()

#===============================================================================
# Prepare Kinase Library Input for Caco-2 ####
#===============================================================================

KINASElib_input_caco2 <- norm_limma_caco2 %>%
  as_tibble(rownames = "Modified.sequence") %>%
  bind_cols(
    as_tibble(fit2_caco2$p.value, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("pval.", .))
  ) %>%
  bind_cols(
    as_tibble(fit2_caco2$coefficients, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("FC.", .))
  ) %>%
  ungroup() %>%
  mutate(
    kinlibseq = fun_modseq_kinlib(Modified.sequence)
  )

#===============================================================================
# Export Kinase Library Data for Caco-2 ####
#===============================================================================

KINASElib_export_caco2 <- KINASElib_input_caco2 %>%
  select(kinlibseq, starts_with("FC."), starts_with("pval.")) 

# Define coefficients list for Caco-2
coefficients_list_caco2 <- list(
  c("FC.SynSign", "pval.SynSign"),
  c("FC.SynSign2", "pval.SynSign2"),
  c("FC.LMvsBSAinDMSO", "pval.LMvsBSAinDMSO"),
  c("FC.LMvsBSAinMEKi", "pval.LMvsBSAinMEKi"),
  c("FC.MEKivsDMSOinBSA", "pval.MEKivsDMSOinBSA"),
  c("FC.MEKivsDMSOinLM", "pval.MEKivsDMSOinLM")
)

# Export Kinase Library Data for each coefficient
walk(coefficients_list, ~ export_kinlib_data(KINASElib_export_caco2, .x, "Caco-2", data_dir))

# Differential Expression Analysis for DLD-1 ####

#===============================================================================
# Differential Expression Analysis for DLD-1 ####
#===============================================================================

# Define column order specific to DLD-1 (excluding D7)
col_order_dld1 <- paste0("D", c(1:6, 8:12))
norm_limma_dld1 <- HpH_norm[, col_order_dld1]

# Define factors for design matrix as per DLD-1 specification
Ligandmix_dld1 <- factor(
  c("BSA", "BSA", "BSA", "BSA", "BSA", "BSA",
    "LigandMix", "LigandMix",
    "LigandMix", "LigandMix", "LigandMix"),
  levels = c("BSA", "LigandMix")
)

MEKi_dld1 <- factor(
  c("DMSO", "DMSO", "DMSO", "MEKi", "MEKi", "MEKi",
    "DMSO", "DMSO", "MEKi", "MEKi", "MEKi"),
  levels = c("DMSO", "MEKi")
)

# Create treatment combinations
TS_dld1 <- paste(Ligandmix_dld1, MEKi_dld1, sep = ".")
TS_dld1 <- factor(TS_dld1, levels = unique(TS_dld1))

# Create design matrix without intercept
design_dld1 <- model.matrix(~0 + TS_dld1)
colnames(design_dld1) <- levels(TS_dld1)

# Fit linear model and apply contrasts
fit_dld1 <- lmFit(norm_limma_dld1, design_dld1)
contrast.matrix_dld1 <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi  = LigandMix.MEKi - BSA.MEKi,
  SynSign        = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2       = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design_dld1
)
fit_dld1 <- contrasts.fit(fit_dld1, contrast.matrix_dld1)
fit2_dld1 <- eBayes(fit_dld1)

# Extract results
results_dld1 <- decideTests(fit2_dld1)
total_ids_dld1 <- colSums(summary(results_dld1))

# Retrieve top tables for all contrasts using the predefined function
LMvsBSAinDMSO_dld1 <- get_top_table(fit2_dld1, 1, total_ids_dld1[1], "LMvsBSAinDMSO")
LMvsBSAinMEKi_dld1  <- get_top_table(fit2_dld1, 2, total_ids_dld1[2], "LMvsBSAinMEKi")
SynSign_dld1        <- get_top_table(fit2_dld1, 3, total_ids_dld1[3], "SynSign")
MEKivsDMSOinBSA_dld1 <- get_top_table(fit2_dld1, 4, total_ids_dld1[4], "MEKivsDMSOinBSA")
MEKivsDMSOinLM_dld1 <- get_top_table(fit2_dld1, 5, total_ids_dld1[5], "MEKivsDMSOinLM")
SynSign2_dld1       <- get_top_table(fit2_dld1, 6, total_ids_dld1[6], "SynSign2")

# Combine differential expression results
df_DLD1 <- norm_limma_dld1 %>%
  as_tibble() %>%
  add_column(F_pval = fit2_dld1$F.p.value, Modified.sequence = rownames(fit2_dld1$p.value)) %>%
  mutate(across(everything(), ~ replace(., is.infinite(.), NA))) %>%
  filter(rowSums(is.na(.)) < 4) %>%
  select(Modified.sequence, F_pval, starts_with("D")) %>%
  left_join(SynSign_dld1, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinDMSO_dld1, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(LMvsBSAinMEKi_dld1, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinBSA_dld1, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(MEKivsDMSOinLM_dld1, by = c("Modified.sequence" = "Modified.Sequence")) %>%
  left_join(SynSign2_dld1, by = c("Modified.sequence" = "Modified.Sequence"))

# Save differential expression results
saveRDS(
  df_DLD1,
  file.path(data_dir, "20241216_CellLinePanel_Global_DLD1_DiffExpr.rds")
)

#===============================================================================
# Filter Significant Peptides for DLD-1 ####
#===============================================================================

full_sign_dld1 <- df_DLD1 %>%
  filter(F_pval < 0.1) %>%
  rename(Modified.Sequence = Modified.sequence)

#===============================================================================
# Prepare Row Annotations for DLD-1 ####
#===============================================================================

row_annotation_DLD1 <- full_sign_dld1 %>%
  # Adjust p-values using Benjamini-Hochberg method and apply -log10 transformation
  mutate(
    `Interaction - Pval` = -log10(p.adjust(P.Value.SynSign, method = "BH")),
    `GFmix w/o MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinDMSO, method = "BH")),
    `GFmix w MEKi - Pval` = -log10(p.adjust(P.Value.LMvsBSAinMEKi, method = "BH")),
    `MEKi w/o GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinBSA, method = "BH")),
    `MEKi w GFmix - Pval` = -log10(p.adjust(P.Value.MEKivsDMSOinLM, method = "BH"))
  ) %>%
  # Rename logFC columns to more descriptive names
  rename(
    `Interaction - logFC` = logFC.SynSign,
    `GFmix w/o MEKi - logFC` = logFC.LMvsBSAinDMSO,
    `GFmix w MEKi - logFC` = logFC.LMvsBSAinMEKi,
    `MEKi w/o GFmix - logFC` = logFC.MEKivsDMSOinBSA,
    `MEKi w GFmix - logFC` = logFC.MEKivsDMSOinLM
  ) %>%
  # Set Modified.Sequence as row names
  column_to_rownames("Modified.Sequence") %>%
  # Select only the relevant annotation columns
  select(
    `Interaction - logFC`, `Interaction - Pval`,
    `GFmix w/o MEKi - logFC`, `GFmix w/o MEKi - Pval`,
    `GFmix w MEKi - logFC`, `GFmix w MEKi - Pval`,
    `MEKi w/o GFmix - logFC`, `MEKi w/o GFmix - Pval`,
    `MEKi w GFmix - logFC`, `MEKi w GFmix - Pval`
  )

row_annotation_DLD1_4f <- row_annotation_DLD1

#===============================================================================
# Scale Data for Heatmap for DLD-1 ####
#===============================================================================

z.mat_dld1 <- full_sign_dld1 %>%
  select(starts_with("D")) %>%
  as.matrix() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  t()

#===============================================================================
# Column Annotations for DLD-1 ####
#===============================================================================

column_annotation_dld1 <- data.frame(
  GFmix = c(rep("BSA", 6), rep("GFmix", 5)),
  MEKi = c(rep("DMSO", 3), rep("MEKi", 3), rep("DMSO", 2), rep("MEKi", 3))
)
rownames(column_annotation_dld1) <- col_order_dld1

col_annotate_mapping_dld1 <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green3", MEKi = "purple")
)

#===============================================================================
# Color Mappings for DLD-1 ####
#===============================================================================

color_mapping_row_dld1 <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w/o MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `GFmix w MEKi - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w/o GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `MEKi w GFmix - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`Interaction - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `GFmix w/o MEKi - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`GFmix w/o MEKi - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `GFmix w MEKi - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`GFmix w MEKi - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `MEKi w/o GFmix - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`MEKi w/o GFmix - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  ),
  `MEKi w GFmix - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`MEKi w GFmix - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  )
)

col_fun_dld1 <- colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))

#===============================================================================
# Create and Save Heatmap for DLD-1 ####
#===============================================================================

set.seed(123)

heatmap_DLD1_LFQ_all <- ComplexHeatmap::Heatmap(
  z.mat_dld1,
  name = "Abundance",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(
    df = column_annotation_dld1,
    col = col_annotate_mapping_dld1
  ),
  right_annotation = HeatmapAnnotation(
    df = row_annotation_DLD1,
    which = "row",
    col = color_mapping_row_dld1,
    annotation_legend_param = list(
      `Interaction - logFC` = list(title = "log2(FC)"),
      `Interaction - Pval` = list(title = "-log10(pval)")
    ),
    show_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  ),
  row_names_gp = gpar(fontsize = 1),
  cluster_rows = TRUE,
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  row_split = 14, # Adjust based on clustering results
  col = col_fun_dld1,
  column_order = col_order_dld1
)

# Save the heatmap to a PDF
pdf(
  file.path(output_dir, "20240215_heat_DLD1_Global_all.pdf"),
  width = 7,
  height = 14
)
draw(heatmap_DLD1_LFQ_all)
clustrows_dld1 <- row_order(heatmap_DLD1_LFQ_all)
dev.off()

#===============================================================================
# Cluster Zoom and Enrichment for DLD-1 ####
#===============================================================================

# Prepare cluster information
df_clusters_dld1 <- clustrows_dld1 %>% 
  enframe(name = "sublist_index") %>%
  mutate(cluster = as.integer(sublist_index)) %>%
  unnest_longer(value, indices_include = TRUE) %>%
  rename(element_in_cluster = value_id, cluster_value = value) %>%
  arrange(cluster_value)

# Combine with scaled data
plotcluster_mean_dld1 <- z.mat_dld1 %>%
  as_tibble() %>%
  bind_cols(df_clusters_dld1) %>%
  mutate(Modified.sequence = full_sign_dld1$Modified.Sequence)

# Save clustered differential expression data
saveRDS(
  plotcluster_mean_dld1,
  file.path(data_dir, "20241216_Global_DLD1_DiffExpr_clustered.rds")
)

# Convert to long format for plotting
df_long_dld1 <- plotcluster_mean_dld1 %>% 
  pivot_longer(cols = D1:D12, names_to = "variable", values_to = "intensity") %>%
  mutate(
    category = case_when(
      variable %in% paste0("D", 1:6) ~ "BSA",
      variable %in% paste0("D", 8:12) ~ "GFmix"
    ),
    treatment = case_when(
      variable %in% paste0("D", c(1:3, 8:9)) ~ "DMSO",
      variable %in% paste0("D", c(4:6, 10:12)) ~ "MEKi"
    )
  )

#===============================================================================
# Calculate Statistics for DLD-1 ####
#===============================================================================

stats_dld1 <- df_long_dld1 %>%
  group_by(cluster, category, treatment) %>%
  summarize(
    mean = mean(intensity, na.rm = TRUE),
    sd = sd(intensity, na.rm = TRUE),
    .groups = 'drop'
  )

#===============================================================================
# Plot Cluster Means for DLD-1 ####
#===============================================================================

meanclusterplots_dld1 <- ggplot(stats_dld1, aes(x = category, y = mean, 
                                                color = treatment, 
                                                group = treatment)) +
  geom_point(position = position_dodge(0.1), size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1)) +
  labs(y = "Mean Value", x = "Category", color = "Treatment") +
  facet_grid(rows = vars(cluster)) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "vertical"
  ) +
  scale_color_manual(values = c("DMSO" = "green3", "MEKi" = "purple"))

# Save the cluster means plot
ggsave(
  file.path(output_dir, "20240215_DLD1_clustmean.pdf"), 
  meanclusterplots_dld1, 
  height = 12, 
  width = 1.5
)

#===============================================================================
# Kinase Library Enrichment for DLD-1 ####
#===============================================================================

# Prepare flanking sequences for enrichment
cluster_DLD1_flanking <- plotcluster_mean_dld1 %>%
  mutate(
    kinlibseq = fun_PrecurIDto_sty(Modified.sequence),
    kinlibseq = map(kinlibseq, split_sequences)
  ) %>%
  unnest(kinlibseq) %>%
  left_join(flankingseqsanno, by = c("kinlibseq" = "phosphopeptide")) %>%
  mutate(
    flanking = if_else(kinlibseq == "AASEVAGVVANAPsPPESSSLCASK", "GVVANAPSPPESSSL", flanking)
  )

# Perform cluster enrichment (adjust cluster numbers as needed)
plot_clusterenrichement(c(10), "DLD1", cluster_DLD1_flanking)
View(readRDS(file.path(data_dir, paste0("DLD1_10.rds")))%>%
       mutate(ER_notlog = Freq_interesting / Freq_background))


#===============================================================================
# Zoomed Heatmap for Specific Cluster in DLD-1 ####
#===============================================================================

#Zoom into cluster 10
DLD1_clust1 <- plotcluster_mean_dld1 %>%
  filter(cluster == 10)

# Prepare row names for heatmap
modseq_tbl_dld1 <- gsub("\\(SILAC-[R|K]-L\\)|[0-9]*$", "", DLD1_clust1$Modified.sequence)
jptseq_tbl_dld1 <- fun_modseq_diannMQlib2pSTY(modseq_tbl_dld1)
precursor_charge <- str_extract(DLD1_clust1$Modified.sequence, "[0-9]*$")
tomatch %>% 
  filter(jptseq %in% jptseq_tbl_dld1) %>%
  group_by(jptseq, plot_genepsite) %>%
  slice(1) %>% ungroup() -> df
df[match(jptseq_tbl_dld1, df$jptseq),] %>%
  pull(plot_genepsite) -> rownames_heatmap_dld1

z.mat_clustzoom_dld1 <- DLD1_clust1 %>%
  dplyr::select(starts_with("D", ignore.case=F)) %>%
  as.matrix()
rownames(z.mat_clustzoom_dld1) <- rownames_heatmap_dld1

# Extract scaled data for the cluster
z.mat_clustzoom_dld1 <- DLD1_clust1 %>%
  select(starts_with("D")) %>%
  as.matrix()
rownames(z.mat_clustzoom_dld1) <- rownames_heatmap_dld1

# Prepare row annotations for the cluster
row_annotation_cluster_dld1 <- row_annotation_DLD1 %>%
  rownames_to_column("Modified.Sequence") %>%
  filter(`Modified.Sequence` %in% DLD1_clust1$Modified.sequence) %>%
  select(`Interaction - logFC`, `Interaction - Pval`)

color_mapping_row_cluster_dld1 <- list(
  `Interaction - logFC` = colorRamp2(c(-1, 0, 1), c("blue3", "white", "red3")),
  `Interaction - Pval` = colorRamp2(
    c(min(row_annotation_DLD1$`Interaction - Pval`, na.rm = TRUE), 3),
    c("white", "black")
  )
)

# Create the zoomed heatmap
heatmap_DLD1_clustzoom <- ComplexHeatmap::Heatmap(
  z.mat_clustzoom_dld1,
  name = "DLD1 Cluster 10",
  show_row_names = TRUE,
  top_annotation = HeatmapAnnotation(
    df = column_annotation_dld1,
    col = col_annotate_mapping_dld1,
    show_legend = FALSE
  ),
  right_annotation = HeatmapAnnotation(
    df = row_annotation_cluster_dld1,
    which = "row",
    col = color_mapping_row_cluster_dld1,
    show_legend = FALSE
  ),
  row_names_gp = gpar(fontsize = 5),
  row_title_rot = 45,
  row_dend_reorder = TRUE,
  col = col_fun_dld1,
  column_order = col_order_dld1,
  show_heatmap_legend = FALSE
)

# Save the zoomed heatmap to a PDF
pdf(
  file.path(output_dir, "heat_Global_DLD1_clust10.pdf"), 
  width = 5, 
  height = 10
)
draw(heatmap_DLD1_clustzoom)
dev.off()

#===============================================================================
# Prepare Kinase Library Input for DLD-1 ####
#===============================================================================

KINASElib_input_dld1 <- norm_limma_dld1 %>%
  as_tibble(rownames = "Modified.sequence") %>%
  bind_cols(
    as_tibble(fit2_dld1$p.value, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("pval.", .))
  ) %>%
  bind_cols(
    as_tibble(fit2_dld1$coefficients, rownames = "Modified.sequence") %>%
      select(-`Modified.sequence`) %>%
      rename_with(~ paste0("FC.", .))
  ) %>%
  ungroup() %>%
  mutate(
    kinlibseq = fun_modseq_kinlib(Modified.sequence)
  )

#===============================================================================
# Export Kinase Library Data for DLD-1 ####
#===============================================================================

KINASElib_export_dld1 <- KINASElib_input_dld1 %>%
  select(kinlibseq, starts_with("FC."), starts_with("pval.")) 

# Define coefficients list for DLD-1
coefficients_list_dld1 <- list(
  c("FC.SynSign", "pval.SynSign"),
  c("FC.SynSign2", "pval.SynSign2"),
  c("FC.LMvsBSAinDMSO", "pval.LMvsBSAinDMSO"),
  c("FC.LMvsBSAinMEKi", "pval.LMvsBSAinMEKi"),
  c("FC.MEKivsDMSOinBSA", "pval.MEKivsDMSOinBSA"),
  c("FC.MEKivsDMSOinLM", "pval.MEKivsDMSOinLM")
)

# Export Kinase Library Data for each coefficient
walk(
  coefficients_list_dld1, 
  ~ export_kinlib_data(
    KINASElib_input_dld1, # Use the tibble/data frame input
    .x, 
    "DLD1", 
    data_dir
  )
)

