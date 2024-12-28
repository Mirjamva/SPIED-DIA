################################################################################
# Targeted Analysis Cell line panel 
# SPIED-DIA HCT116, DLD-1, Caco2 with GFMix/MEKi
#
# Author: Mirjam
# Date:   2024-12-01
#
# Description:
# This script performs a targeted phosphopeptide analysis from the DIA-NN SILAC dataset.
# analysis of cell-line panel (HCT116, DLD-1, Caco2) treated with/without MEKi &
# GFmix. Produces figures 3A-C, provides data underlying figure 3D/E and supplementary
# figure 7C, 8. Logic of limma model is explained in supplementary figure 7A.
# 
# It includes:
#  - Loading and cleaning data from target analysis annotated with HGNC IDs.
#  - Normalizing data using cyclic loess normalization.
#  - SPIED-DIA identification and quantification of target peptides
#  - Statistical testing of conditions using limma (linear models and contrasts).
#  - Generating heatmaps of significant changes.
#  
#  Per cell line the following approach is followed:
#  
#  - Filter and reshape data
#  - Define design matrices and contrasts
#  - Fit linear models with limma
#  - Extract results and map to HGNC IDs
#  - Generate heatmaps
#
# The script is designed to be executed in an environment where the required 
# input data files are available in the "Data" directory, and outputs (plots and 
# results) are stored in "Plots" and "Results" directories, respectively.

################################################################################
# ---------------------------- #
#         Libraries            #
# ---------------------------- #

# Load necessary libraries
library(tidyverse)
library(data.table)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(ggh4x)
library(pheatmap)

# ---------------------------- #
#      Configuration Params    #
# ---------------------------- #

# Define paths (update these paths as needed)
data_dir <- "Data"
plots_dir <- "Plots"
results_dir <- "Results"

# Ensure directories exist
dir.create(data_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

# File paths
# CSV file with annotated JPT peptides (supplementaray table 1)
annotated_peptides_file <- file.path(data_dir, "jptpeps_20240208_HGNC_annotated.csv")
# TSV file (DIA-NN "report.tsv") containing SILAC target data:
silac_report_file <- file.path(data_dir, "SILAC_targetlib", "report.tsv")

# ---------------------------- #
#       Function Definitions   #
# ---------------------------- #

# Function to clean modified sequences
fun_modseq_diannMQlib2pSTY <- function(seq_p) {
  # Remove all underscores
  seq_p <- gsub("_", "", seq_p)
  
  # Replace S/Y/T with their phosphorylated forms (pS, pY, pT)
  # This captures S/Y/T followed by either (Phospho (STY)) or (UniMod:21)
  seq_p <- gsub("([SYT])\\((Phospho \\(STY\\)|UniMod:21)\\)", "p\\1", seq_p)
  
  # Remove unwanted patterns:
  # (SILAC-R-H), (SILAC-K-H), (SILAC-R-L), (SILAC-K-L),
  # (Oxidation (M)), (Acetyl (Protein N-term)), (UniMod:\d+)
  seq_p <- gsub("\\(SILAC-[RK]-[HL]\\)|\\(Oxidation \\(M\\)\\)|\\(Acetyl \\(Protein N-term\\)\\)|\\(UniMod:\\d+\\)", "", seq_p)
  
  return(seq_p)
}

# Function to create ordered levels for runs
create_ordered_levels <- function() {
  vector_228_C <- paste0("228_C", 1:12)
  vector_309_D <- paste0("309_D", 1:12)
  vector_228_H <- paste0("228_H", 1:12)
  
  ordered_levels <- c(vector_228_C, vector_309_D, vector_228_H)
  return(ordered_levels)
}

# Function to normalize data using cyclic loess
# normalize per cell line
normalize_data <- function(JPT_m, loess_span = 0.7, loess_it = 3) {
  jpt_norm_loess_H <- normalizeCyclicLoess(
    JPT_m[, grep("^H", colnames(JPT_m), value = TRUE)],
    weights = NULL, 
    span = loess_span, 
    iterations = loess_it, 
    method = "fast"
  )
  
  jpt_norm_loess_D <- normalizeCyclicLoess(
    JPT_m[, grep("^D", colnames(JPT_m), value = TRUE)],
    weights = NULL, 
    span = loess_span, 
    iterations = loess_it, 
    method = "fast"
  )
  
  jpt_norm_loess_C <- normalizeCyclicLoess(
    JPT_m[, grep("^C", colnames(JPT_m), value = TRUE)],
    weights = NULL, 
    span = loess_span, 
    iterations = loess_it, 
    method = "fast"
  )
  
  normalized <- cbind(jpt_norm_loess_H, jpt_norm_loess_D, jpt_norm_loess_C) %>%
    as.data.frame() %>% 
    rownames_to_column("Precursor.Id") %>%
    pivot_longer(
      cols = starts_with("C") | starts_with("D") | starts_with("H"), 
      names_to = "Run_simple",
      values_to = "norm_intensity"
    )
  
  return(normalized)
}


# ---------------------------- #
#       Data Loading           #
# ---------------------------- #

# Read annotated peptides
jptpeps_annotated <- fread(
  annotated_peptides_file, 
  sep = ",", 
  header = TRUE, 
  stringsAsFactors = FALSE, 
  check.names = FALSE
) %>% 
  rowwise()

# Extract unique pSTY sequences
ids_ordered <- jptpeps_annotated %>%
  pull(pSTYseq) %>% 
  unique()

# Read SILAC report
report_SILACtarget <- fread(
  silac_report_file,
  sep = "\t"
)

# ----------------------------- #
#       Data Preprocessing      #
# ----------------------------- #

# Preprocess DIA-NN report
jpt_fil <- report_SILACtarget %>%
  ungroup() %>%
  mutate(
    SILAClabel = ifelse(
      grepl("-H|-L", Modified.Sequence), 
      str_extract(Modified.Sequence, "-H|-L"),
      "LFQ"
    ),
    jptseq = fun_modseq_diannMQlib2pSTY(Modified.Sequence),
    Run_simple = gsub("Popeye_20230|_MVB_HStdia_SS2|_[A-Z0-9]*_[0-9]_[0-9]*$", "", Run)
  ) %>%
  filter(!grepl("Popeye_20230220", Run), 
         !grepl("Popeye_20230228_MVB_HStdia_SS2_C6_A6_1_3273", Run)) %>% #remove outlier
  mutate(
    sample_simple = Run_simple,
    sample_simple = case_when(
      grepl("10$|11$|12$", sample_simple) ~ gsub("10$|11$|12$", "_GF_MEKi", sample_simple),
      grepl("1$|2$|3$", sample_simple) ~ gsub("1$|2$|3$", "_BSA_DMSO", sample_simple),
      grepl("4$|5$|6$", sample_simple) ~ gsub("4$|5$|6$", "_GF_DMSO", sample_simple),
      grepl("7$|8$|9$", sample_simple) ~ gsub("7$|8$|9$", "_BSA_MEKi", sample_simple),
      TRUE ~ sample_simple
    )
  )

ordered_levels <- create_ordered_levels()

# ----------------------------------------------- #
#### SPIED-DIA identification & quantification ####
# ----------------------------------------------- #

# Extract SILAC-L IDs
# L ids medium confidently identified
L_ids <- jpt_fil %>%
  filter(SILAClabel == "-L", Channel.Q.Value < 0.5) %>% 
  mutate(
    run_precursor_id = paste0(
      Run, 
      gsub("-L\\)", "-H\\)", Precursor.Id)
    )
  ) %>%
  pull(run_precursor_id)

# Normalize and filter data
jpt_2norm <- jpt_fil %>% 
  filter(
    Channel.Q.Value < 0.05,
    PTM.Q.Value < 0.05,
    SILAClabel == "-H",
    Channel.H > 1000
  ) %>% 
  mutate(
    run_precursor_id = paste0(Run, Precursor.Id),
    cell_line = str_extract(sample_simple, "C|D|H"),
    L_idd = run_precursor_id %in% L_ids
  ) %>% # remove right? filter(grepl("VADPEHDHTGFLTEYVA", Precursor.Id)) %>% 
  group_by(Precursor.Id, cell_line) %>%
  filter(n() > 9) %>%
  mutate(rescalingfct = median(Ms1.Area)) %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(rescaled_int = log10(Channel.L / Channel.H * rescalingfct)) %>% 
  # Target id needs to be identified confidently in the equivalent of at least one experimental group
  group_by(cell_line, Precursor.Id) %>%
  mutate(
    C_count = sum(cell_line == "C" & rescaled_int > 2),
    D_count = sum(cell_line == "D" & rescaled_int > 2),
    H_count = sum(cell_line == "H" & rescaled_int > 2)
  ) %>% 
  ungroup() %>%  
  filter(
    (cell_line == "C" & C_count > 1 & sum(L_idd) > 1) |
      (cell_line == "D" & D_count > 1 & sum(L_idd) > 1) |
      (cell_line == "H" & H_count > 2 & sum(L_idd) > 2)
  ) %>%
  ungroup()

# Save normalized data (intermediate file)
saveRDS(jpt_2norm, 
        file = file.path(results_dir, "target_results_prenorm.rds"), 
        compress = "gzip")

# visualise heavy ID's & light ID's 
jpt_fil %>%
  filter(
    Channel.Q.Value < 0.05,
    PTM.Q.Value < 0.05,
    SILAClabel == "-H",
    Channel.H > 1000
  ) %>%
  mutate(CellLine = case_when(grepl("H", sample_simple) ~ "HCT116",
                              grepl("C", sample_simple) ~ "Caco2",
                              grepl("D", sample_simple) ~ "DLD-1"), 
         Id = "Spike-In") %>%
  dplyr::select(Precursor.Id, jptseq, CellLine, Run_simple, Id) -> H_ids

jpt_2norm %>%
  mutate(CellLine = case_when(grepl("H", sample_simple) ~ "HCT116",
                              grepl("C", sample_simple) ~ "Caco2",
                              grepl("D", sample_simple) ~ "DLD-1"), 
         Id = "Endogenous") %>%
  dplyr::select(Precursor.Id, jptseq, CellLine, Run_simple, Id) -> L_idsfinal

comp_HL <- rbind(H_ids, L_idsfinal)

# Summarize the data, calculating the number of IDs per 'CellLine' and 'Id' category (target and endogenous peptides)
summary_HL <- comp_HL %>%
  group_by(Run_simple, CellLine, Id) %>%
  summarise(IDs = n_distinct(jptseq), .groups = 'drop') %>% 
  group_by(CellLine, Id) %>%
  mutate(Avg_IDs = mean(IDs),  
         Std_Dev = sd(IDs))   

summary_HL$CellLine <- factor(summary_HL$CellLine, levels = c( "HCT116","DLD-1", "Caco2"))  

# Plot the data using ggplot2 with a similar style as provided
compHL_analysis <- ggplot(summary_HL, aes(x = CellLine, y = IDs, color = Id)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.5, draw_quantiles = c(0.5)) +  
  geom_point(position = position_jitterdodge(seed = 1), alpha = 0.5) +  
  scale_color_manual(values = c("Spike-In" = "darkred",
                                "Endogenous" = "darkblue")) +  
  labs(y = "N Identified") +  
  theme_minimal(base_size = 14) + 
  ylim(c(0, 200)) +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "bottom",  
        legend.direction = "vertical",  
        axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave(compHL_analysis, 
       filename=file.path(plots_dir, "Target_HL_Ids.pdf"), 
       width = 2, height = 3.5)

# ---------------------------- #
#        Normalization         #
# ---------------------------- #

# Pivot normalized data for normalization
JPT_m <- jpt_2norm %>% 
  mutate(
    log_int = rescaled_int,
    Run_simple = gsub("^[0-9]*_", "", Run_simple)
  ) %>%
  select(Precursor.Id, Run_simple, log_int) %>%
  pivot_wider(names_from = Run_simple, values_from = log_int) %>%
  column_to_rownames("Precursor.Id") %>% 
  as.matrix()

# Normalize using cyclic loess
normalized_data <- normalize_data(JPT_m)

# Save normalized matrix
saveRDS(normalized_data, 
        file = file.path(results_dir, "target_results_normalized.rds"),
        compress = "gzip")

################################################################################
# Statistical Analysis
# Differential Expression in HCT116, DLD-1, and Caco-2 Cell Lines
#
# In this section, we use the normalized SILAC data to:
#  - Fit linear models (via limma) to identify differentially abundant peptides.
#  - Define contrasts for various conditions (e.g., BSA vs. LigandMix, DMSO vs. MEKi).
#  - Extract top-ranked peptides for each contrast.
#  - Map sequences to HGNC IDs for gene-level interpretation.
#  - Generate heatmaps visualizing significant changes.
################################################################################

# Parameters for filtering and visualization
Fpval_cutoff <- 0.1
universal_font_size <- 10

# Define bold phosphosites for annotation
bold_phosphosites <- c("JNK3_Y223;JNK1_Y185", 
                       "ERK1_Y204", "ERK2_Y187")

# ----------------------------- #
#          HCT116 Analysis      #
# ----------------------------- #

# Step 1: Filter and reshape data for HCT116
# Selecting and formatting the data matrix for linear modeling
hct116_data <- normalized_data %>%
  filter(str_detect(Run_simple, "^H")) %>% 
  group_by(Precursor.Id) %>% 
  filter(n() > 5) %>% 
  ungroup() %>% 
  select(Precursor.Id, Run_simple, norm_intensity) %>% 
  pivot_wider(names_from = Precursor.Id, values_from = norm_intensity) %>%
  arrange(factor(Run_simple, levels = paste0("H", 1:12))) %>%
  column_to_rownames("Run_simple") %>%
  as.matrix() %>%
  t()  

# Step 2: Define experimental design for HCT116
ligandmix_levels <- c(rep("BSA", 6), rep("LigandMix", 6))
meki_levels <- c(rep("DMSO", 3), rep("MEKi", 3),
                 rep("DMSO", 3), rep("MEKi", 3))

Ligandmix <- factor(ligandmix_levels, levels = c("BSA", "LigandMix"))
MEKi <- factor(meki_levels, levels = c("DMSO", "MEKi"))

TS <- factor(paste(Ligandmix, MEKi, sep = "."), 
             levels = unique(paste(ligandmix_levels, meki_levels, sep = ".")))

design <- model.matrix(~0 + TS)
colnames(design) <- levels(TS)

# Step 3: Define contrasts for HCT116
contrast_matrix_hct116 <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi = LigandMix.MEKi - BSA.MEKi,
  SynSign = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2 = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design
)

# Step 4: Fit model and apply contrasts for HCT116
fit_hct116 <- lmFit(hct116_data, design)
fit_hct116 <- contrasts.fit(fit_hct116, contrast_matrix_hct116)
fit2_hct116 <- eBayes(fit_hct116)

# Step 5: Extract Top Tables for Each Contrast
contrasts <- colnames(contrast_matrix_hct116)
top_tables_hct116 <- lapply(seq_along(contrasts), function(i) {
  topTable(fit2_hct116, coef = i, adjust.method = "BH", number = Inf) %>%
    rename_with(~ paste0(., ".", contrasts[i])) %>%
    rownames_to_column("Precursor.Id")
})

combined_top_tables_hct116 <- reduce(top_tables_hct116, left_join, by = "Precursor.Id")

# Step 6: Compile Final Dataframe with F P-Values
df_hct116_final <- as_tibble(hct116_data, rownames = "Precursor.Id") %>%
  mutate(
    F_pval = fit2_hct116$F.p.value,
    Modified.sequence = rownames(fit2_hct116$p.value)
  ) %>% 
  mutate(across(where(is.numeric), ~ replace(., is.infinite(.), NA))) %>% # Replace Inf with NA
  filter(rowSums(is.na(.)) < 6) %>%  # Retain rows with fewer than 6 NAs
  select(Modified.sequence, F_pval, starts_with("H")) %>%
  left_join(combined_top_tables_hct116, by = c("Modified.sequence" = "Precursor.Id"))

# Step 7: Map Modified.sequence to HGNC_id
jptseq_tbl <- fun_modseq_diannMQlib2pSTY(
  gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", rownames(fit2_hct116$p.value))
)

sequence_mapping <- jptpeps_annotated %>%
  filter(pSTYseq %in% jptseq_tbl) %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", pSTYseq)
    )
  ) %>%
  select(Processed_seq, HGNC_id)

df_hct116_final <- df_hct116_final %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", Modified.sequence)
    )
  ) %>%
  left_join(sequence_mapping, by = c("Processed_seq" = "Processed_seq")) %>%
  select(-Processed_seq)

if(any(is.na(df_hct116_final$HGNC_id))) {
  warning("Some sequences could not be mapped to HGNC_id.")
}

write.csv(
  df_hct116_final, 
  file = file.path(results_dir, "df_hct116_final.csv"), 
  row.names = FALSE
)

# ----------------------------- #
#     Preparing the Heatmap     #
# ----------------------------- #

# Step 8: Prepare Data for Heatmap
heatmap_data <- df_hct116_final %>%
  filter(F_pval < Fpval_cutoff) %>%
  group_by(HGNC_id) %>%
  slice_min(order_by = F_pval, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    HGNC_id = case_when(
      F_pval < 0.005 ~ paste0(HGNC_id, "***"),
      F_pval < 0.01  ~ paste0(HGNC_id, "**"),
      F_pval < 0.05  ~ paste0(HGNC_id, "*"),
      TRUE            ~ HGNC_id
    )
  ) %>%
  select(HGNC_id, starts_with("H")) %>%
  column_to_rownames("HGNC_id") 

heatmap_matrix <- as.matrix(heatmap_data)

# Step 9: Scale Row-wise
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix), center = TRUE, scale = TRUE))

# Step 10: Create Heatmap annotations (HCT116)
MEKi <- factor(c("DMSO","DMSO", "DMSO", "MEKi", "MEKi", "MEKi",
                 "DMSO","DMSO","DMSO", "MEKi", "MEKi", "MEKi"), 
               levels = c("DMSO","MEKi"))

column_annotation_hct116 <- data.frame(
  GFmix = c("BSA","BSA","BSA","BSA","BSA","BSA",
            "GFmix","GFmix","GFmix","GFmix","GFmix","GFmix"),
  MEKi = MEKi
)
rownames(column_annotation_hct116) <- colnames(heatmap_matrix_scaled)

annotation_colors_hct116 <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green", MEKi = "purple")
)

column_annotation_object_hct116 <- HeatmapAnnotation(
  df = column_annotation_hct116, 
  col = annotation_colors_hct116, 
  annotation_name_gp = gpar(fontsize = universal_font_size)
)

is_bold_hct116 <- sapply(gsub("\\*", "", rownames(heatmap_matrix_scaled)), function(seq) {
  any(grepl(seq, gsub("\\*", "", bold_phosphosites)))
})

row_annotation_hct116 <- rowAnnotation(
  Modified.sequence = anno_text(
    rownames(heatmap_matrix_scaled), 
    rot = 0, 
    gp = gpar(
      fontface = ifelse(is_bold_hct116, "bold", "plain"), 
      fontsize = universal_font_size
    )
  )
)

# Step 11: Generate and Save Heatmap
heatmap_plot_hct116 <- Heatmap(
  heatmap_matrix_scaled, 
  name = "Abundance", 
  show_row_names = FALSE,
  show_column_names = TRUE, 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  column_title = "HCT116", 
  column_title_gp = gpar(fontsize = universal_font_size + 2, fontface = "bold"),
  top_annotation = column_annotation_object_hct116,
  right_annotation = row_annotation_hct116,
  col = colorRamp2(c(-2, 0, 2), c("blue", "ivory", "red")),
  rect_gp = gpar(col = "white", lwd = 1),
  row_split = 4,  # Adjust based on your data
  height = unit(nrow(heatmap_matrix_scaled) * 5, "mm"),
  width = unit(ncol(heatmap_matrix_scaled) * 5, "mm"),
  heatmap_legend_param = list(
    title = "Abundance", 
    position = "bottom", 
    title_position = "topcenter", 
    direction = "horizontal"
  )
)

heatmap_pdf_path_hct116 <- file.path(
  plots_dir, 
  paste0("heatmapTarget_HCT116_Fpval", Fpval_cutoff, ".pdf")
)

# Save the Heatmap to PDF
pdf(heatmap_pdf_path_hct116, width = 10, height = 10)
draw(heatmap_plot_hct116)
dev.off()

# ---------------------------- #
#          DLD-1 Analysis      #
# ---------------------------- #
# Step 1: Filter and Reshape Data for DLD-1
dld1_data <- normalized_data %>% 
  filter(grepl("^D", Run_simple)) %>%      
  group_by(Precursor.Id) %>% 
  filter(n() > 4) %>%                 
  ungroup() %>% 
  select(Precursor.Id, Run_simple, norm_intensity) %>% 
  pivot_wider(names_from = Precursor.Id, values_from = norm_intensity) %>%
  arrange(factor(Run_simple, levels = paste0("D", c(1:6, 8:12)))) %>%
  column_to_rownames("Run_simple") %>%
  as.matrix() %>%
  t() 

# Step 2: Define Experimental Design
Ligandmix_dld1 <- factor(c("BSA","BSA","BSA","BSA","BSA","BSA",
                           "LigandMix","LigandMix",
                           "LigandMix","LigandMix","LigandMix"),
                         levels = c("BSA","LigandMix"))
MEKi_dld1 <- factor(c("DMSO","DMSO", "DMSO", "MEKi", "MEKi", "MEKi",
                      "DMSO","DMSO", "MEKi", "MEKi", "MEKi"), 
                    levels = c("DMSO","MEKi"))
TS_dld1 <- factor(paste(Ligandmix_dld1, MEKi_dld1, sep = "."), 
                  levels = unique(paste(Ligandmix_dld1, MEKi_dld1, sep = ".")))

design_dld1 <- model.matrix(~0 + TS_dld1)
colnames(design_dld1) <- levels(TS_dld1)

# Step 3: Define Contrasts
contrast_matrix_dld1 <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi = LigandMix.MEKi - BSA.MEKi,
  SynSign = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2 = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design_dld1
)

# Step 4: Fit the Model and Apply Contrasts
fit_dld1 <- lmFit(dld1_data, design_dld1)
fit_dld1 <- contrasts.fit(fit_dld1, contrast_matrix_dld1)
fit2_dld1 <- eBayes(fit_dld1)

# Step 5: Extract Top Tables for Each Contrast
contrasts_dld1 <- colnames(contrast_matrix_dld1)
top_tables_dld1 <- lapply(seq_along(contrasts_dld1), function(i) {
  topTable(fit2_dld1, coef = i, adjust.method = "BH", number = Inf) %>%
    rename_with(~ paste0(., ".", contrasts_dld1[i])) %>%
    rownames_to_column("Precursor.Id")
})

# Combine all top tables into a single dataframe
combined_top_tables_dld1 <- reduce(top_tables_dld1, left_join, by = "Precursor.Id")

# Step 6: Compile Final Dataframe with F P-Values
df_dld1_final <- as_tibble(dld1_data, rownames = "Precursor.Id") %>%
  mutate(
    F_pval = fit2_dld1$F.p.value,
    Modified.sequence = rownames(fit2_dld1$p.value)
  ) %>% 
  mutate(across(where(is.numeric), ~ replace(., is.infinite(.), NA))) %>% 
  filter(rowSums(is.na(.)) < 6) %>%  
  select(Modified.sequence, F_pval, starts_with("D")) %>%
  left_join(combined_top_tables_dld1, by = c("Modified.sequence" = "Precursor.Id"))


# Step 7: Map Modified.sequence to HGNC_id

jptseq_tbl_dld1 <- fun_modseq_diannMQlib2pSTY(
  gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", rownames(fit2_dld1$p.value))
)

sequence_mapping_dld1 <- jptpeps_annotated %>%
  filter(pSTYseq %in% jptseq_tbl_dld1) %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", pSTYseq)
    )
  ) %>%
  select(Processed_seq, HGNC_id)

df_dld1_final <- df_dld1_final %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", Modified.sequence)
    )
  ) %>%
  left_join(sequence_mapping_dld1, by = c("Processed_seq" = "Processed_seq")) %>%
  select(-Processed_seq) 

write.csv(
  df_dld1_final, 
  file = file.path(results_dir, "df_dld1_final.csv"), 
  row.names = FALSE
)

# ---------------------------- #
#     Preparing the Heatmap     #
# ---------------------------- #

# Step 8: Prepare Data for Heatmap
heatmap_data_dld1 <- df_dld1_final %>%
  filter(F_pval < Fpval_cutoff) %>%
  group_by(HGNC_id) %>%
  slice_min(order_by = F_pval, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    HGNC_id = case_when(
      F_pval < 0.005 ~ paste0(HGNC_id, "***"),
      F_pval < 0.01  ~ paste0(HGNC_id, "**"),
      F_pval < 0.05  ~ paste0(HGNC_id, "*"),
      TRUE            ~ HGNC_id
    )
  ) %>%
  select(HGNC_id, starts_with("D")) %>%
  column_to_rownames("HGNC_id") 

heatmap_matrix_dld1 <- as.matrix(heatmap_data_dld1)

# Step 9: Scale Row-wise
heatmap_matrix_dld1_scaled <- t(scale(t(heatmap_matrix_dld1), center = TRUE, scale = TRUE))

# Step 10: Create Annotations

# Define MEKi factor based on experimental design
MEKi_dld1_annotation <- factor(c("DMSO","DMSO", "DMSO", "MEKi", "MEKi", "MEKi",
                                 "DMSO","DMSO", "MEKi", "MEKi", "MEKi"), 
                               levels = c("DMSO","MEKi"))

# Define column annotations
column_annotation_dld1 <- data.frame(
  GFmix = c(rep("BSA", 6), rep("GFmix", 5)),  
  MEKi = MEKi_dld1_annotation
)
rownames(column_annotation_dld1) <- colnames(heatmap_matrix_dld1_scaled)

# Define colors for annotations
annotation_colors_dld1 <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green", MEKi = "purple")
)

# Create ComplexHeatmap annotation object
column_annotation_object_dld1 <- HeatmapAnnotation(
  df = column_annotation_dld1, 
  col = annotation_colors_dld1, 
  annotation_name_gp = gpar(fontsize = universal_font_size)
)

# Create a logical vector for bold sequences, disregarding asterisks
is_bold_dld1 <- sapply(gsub("\\*", "", rownames(heatmap_matrix_dld1_scaled)), function(seq) {
  any(grepl(seq, gsub("\\*", "", bold_phosphosites)))
})

# Define row annotations for bold sequences
row_annotation_dld1 <- rowAnnotation(
  Modified.sequence = anno_text(
    rownames(heatmap_matrix_dld1_scaled), 
    rot = 0, 
    gp = gpar(
      fontface = ifelse(is_bold_dld1, "bold", "plain"), 
      fontsize = universal_font_size
    )
  )
)

# Step 11: Generate and Save Heatmap
heatmap_plot_dld1 <- Heatmap(
  heatmap_matrix_dld1_scaled, 
  name = "Abundance", 
  show_row_names = FALSE,
  show_column_names = TRUE, 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  column_title = "DLD-1", 
  column_title_gp = gpar(fontsize = universal_font_size + 2, fontface = "bold"),
  top_annotation = column_annotation_object_dld1,
  right_annotation = row_annotation_dld1,
  col = colorRamp2(c(-2, 0, 2), c("blue", "ivory", "red")),
  rect_gp = gpar(col = "white", lwd = 1),
  row_split = 4,  # Adjust based on your data
  height = unit(nrow(heatmap_matrix_dld1_scaled) * 5, "mm"),
  width = unit(ncol(heatmap_matrix_dld1_scaled) * 5, "mm"),
  heatmap_legend_param = list(
    title = "Abundance", 
    position = "bottom", 
    title_position = "topcenter", 
    direction = "horizontal"
  )
)


heatmap_pdf_path_dld1 <- file.path(
  plots_dir, 
  "20241202_heatmapTarget_DLD1_Fpval_0.1.pdf"
)

pdf(heatmap_pdf_path_dld1, width = 10, height = 10)
draw(heatmap_plot_dld1)
dev.off()

cat("DLD-1 Heatmap saved to:", heatmap_pdf_path_dld1, "\n")

# ---------------------------- #
####     Caco-2 Analysis    ####
# ---------------------------- #

# Step 1: Filter and Reshape Data for Caco-2
caco2_data <- normalized_data %>% 
  filter(grepl("^C", Run_simple)) %>%         
  filter(!grepl("220", Run_simple)) %>%   
  group_by(Precursor.Id) %>% 
  filter(n() > 4) %>%                          
  ungroup() %>% 
  select(Precursor.Id, Run_simple, norm_intensity) %>% 
  pivot_wider(names_from = Precursor.Id, values_from = norm_intensity) %>%
  arrange(factor(Run_simple, levels = paste0("C", c(1:5, 7:12)))) %>%
  column_to_rownames("Run_simple") %>%
  as.matrix() %>%
  t() 

# Step 2: Define Experimental Design
Ligandmix_caco2 <- factor(c("BSA","BSA","BSA","BSA","BSA",
                            "LigandMix","LigandMix","LigandMix",
                            "LigandMix","LigandMix","LigandMix"),
                          levels = c("BSA","LigandMix"))
MEKi_caco2 <- factor(c("DMSO","DMSO", "DMSO", "MEKi", "MEKi", 
                       "DMSO","DMSO","DMSO", "MEKi", "MEKi", "MEKi"), 
                     levels = c("DMSO","MEKi"))
TS_caco2 <- factor(paste(Ligandmix_caco2, MEKi_caco2, sep = "."), 
                   levels = unique(paste(Ligandmix_caco2, MEKi_caco2, sep = ".")))

design_caco2 <- model.matrix(~0 + TS_caco2)
colnames(design_caco2) <- levels(TS_caco2)

# Step 3: Define Contrasts
contrast_matrix_caco2 <- makeContrasts(
  LMvsBSAinDMSO = LigandMix.DMSO - BSA.DMSO,
  LMvsBSAinMEKi = LigandMix.MEKi - BSA.MEKi,
  SynSign = (LigandMix.MEKi - BSA.MEKi) - (LigandMix.DMSO - BSA.DMSO),
  MEKivsDMSOinBSA = BSA.MEKi - BSA.DMSO,
  MEKivsDMSOinLM = LigandMix.MEKi - LigandMix.DMSO,
  SynSign2 = (LigandMix.MEKi - LigandMix.DMSO) - (BSA.MEKi - BSA.DMSO),
  levels = design_caco2
)

# Step 4: Fit the Model and Apply Contrasts
fit_caco2 <- lmFit(caco2_data, design_caco2)
fit_caco2 <- contrasts.fit(fit_caco2, contrast_matrix_caco2)
fit2_caco2 <- eBayes(fit_caco2)

# Step 5: Extract Top Tables for Each Contrast
contrasts_caco2 <- colnames(contrast_matrix_caco2)
top_tables_caco2 <- lapply(seq_along(contrasts_caco2), function(i) {
  topTable(fit2_caco2, coef = i, adjust.method = "BH", number = Inf) %>%
    rename_with(~ paste0(., ".", contrasts_caco2[i])) %>%
    rownames_to_column("Precursor.Id")
})

# Combine all top tables into a single dataframe
combined_top_tables_caco2 <- reduce(top_tables_caco2, left_join, by = "Precursor.Id")

# Step 6: Compile Final Dataframe with F P-Values
df_caco2_final <- as_tibble(caco2_data, rownames = "Precursor.Id") %>%
  mutate(
    F_pval = fit2_caco2$F.p.value,
    Modified.sequence = rownames(fit2_caco2$p.value)
  ) %>% 
  mutate(across(where(is.numeric), ~ replace(., is.infinite(.), NA))) %>% 
  filter(rowSums(is.na(.)) < 6) %>%  
  select(Modified.sequence, F_pval, starts_with("C")) %>%
  left_join(combined_top_tables_caco2, by = c("Modified.sequence" = "Precursor.Id"))

# ---------------------------- #
#     Mapping to HGNC_id         #
# ---------------------------- #

# Step 7: Map Modified.sequence to HGNC_id

# Process sequences using custom function
jptseq_tbl_caco2 <- fun_modseq_diannMQlib2pSTY(
  gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", rownames(fit2_caco2$p.value))
)

# Create a mapping dataframe
sequence_mapping_caco2 <- jptpeps_annotated %>%
  filter(pSTYseq %in% jptseq_tbl_caco2) %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", pSTYseq)
    )
  ) %>%
  select(Processed_seq, HGNC_id)

# Merge HGNC_id into df_caco2_final
df_caco2_final <- df_caco2_final %>%
  mutate(
    Processed_seq = fun_modseq_diannMQlib2pSTY(
      gsub("\\(SILAC-(K|R)-(H|L)\\)[234]$", "", Modified.sequence)
    )
  ) %>%
  left_join(sequence_mapping_caco2, by = c("Processed_seq" = "Processed_seq")) %>%
  select(-Processed_seq)  

write.csv(
  df_caco2_final, 
  file = file.path(results_dir, "df_caco2_final.csv"), 
  row.names = FALSE
)

# ---------------------------- #
#     Preparing the Heatmap     #
# ---------------------------- #

# Step 8: Prepare Data for Heatmap
heatmap_data_caco2 <- df_caco2_final %>%
  filter(F_pval < Fpval_cutoff) %>%
  group_by(HGNC_id) %>%
  slice_min(order_by = F_pval, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    HGNC_id = case_when(
      F_pval < 0.005 ~ paste0(HGNC_id, "***"),
      F_pval < 0.01  ~ paste0(HGNC_id, "**"),
      F_pval < 0.05  ~ paste0(HGNC_id, "*"),
      TRUE            ~ HGNC_id
    )
  ) %>%
  select(HGNC_id, starts_with("C")) %>%
  column_to_rownames("HGNC_id") 

# Convert to matrix
heatmap_matrix_caco2 <- as.matrix(heatmap_data_caco2)

# Step 9: Scale Row-wise
heatmap_matrix_caco2_scaled <- t(scale(t(heatmap_matrix_caco2), center = TRUE, scale = TRUE))

# Step 10: Create Annotations

# Define MEKi factor based on experimental design
MEKi_caco2_annotation <- factor(c("DMSO","DMSO", "DMSO", "MEKi", "MEKi", 
                                  "DMSO","DMSO","DMSO", "MEKi", "MEKi", "MEKi"), 
                                levels = c("DMSO","MEKi"))

# Define column annotations
column_annotation_caco2 <- data.frame(
  GFmix = c(rep("BSA", 6), rep("GFmix", 5)),  
  MEKi = MEKi_caco2_annotation
)
rownames(column_annotation_caco2) <- colnames(heatmap_matrix_caco2_scaled)

# Define colors for annotations
annotation_colors_caco2 <- list(
  GFmix = c(BSA = "blue", GFmix = "red"),
  MEKi = c(DMSO = "green", MEKi = "purple")
)

# Create ComplexHeatmap annotation object
column_annotation_object_caco2 <- HeatmapAnnotation(
  df = column_annotation_caco2, 
  col = annotation_colors_caco2, 
  annotation_name_gp = gpar(fontsize = universal_font_size)
)

# Create a logical vector for bold sequences, disregarding asterisks
is_bold_caco2 <- sapply(gsub("\\*", "", rownames(heatmap_matrix_caco2_scaled)), function(seq) {
  any(grepl(seq, gsub("\\*", "", bold_phosphosites)))
})

# Define row annotations for bold sequences
row_annotation_caco2 <- rowAnnotation(
  Modified.sequence = anno_text(
    rownames(heatmap_matrix_caco2_scaled), 
    rot = 0, 
    gp = gpar(
      fontface = ifelse(is_bold_caco2, "bold", "plain"), 
      fontsize = universal_font_size
    )
  )
)

# Step 11: Generate and Save Heatmap
heatmap_plot_caco2 <- Heatmap(
  heatmap_matrix_caco2_scaled, 
  name = "Abundance", 
  show_row_names = FALSE,
  show_column_names = TRUE, 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  column_title = "Caco-2", 
  column_title_gp = gpar(fontsize = universal_font_size + 2, fontface = "bold"),
  top_annotation = column_annotation_object_caco2,
  right_annotation = row_annotation_caco2,
  col = colorRamp2(c(-2, 0, 2), c("blue", "ivory", "red")),
  rect_gp = gpar(col = "white", lwd = 1),
  row_split = 4, 
  height = unit(nrow(heatmap_matrix_caco2_scaled) * 5, "mm"),
  width = unit(ncol(heatmap_matrix_caco2_scaled) * 5, "mm"),
  heatmap_legend_param = list(
    title = "Abundance", 
    position = "bottom", 
    title_position = "topcenter", 
    direction = "horizontal"
  )
)

# Define the output PDF file path
heatmap_pdf_path_caco2 <- file.path(
  plots_dir, 
  "20241202_heatmapTarget_Caco2_Fpval_0.1.pdf"
)

# Save the Heatmap to PDF
pdf(heatmap_pdf_path_caco2, width = 10, height = 10)
draw(heatmap_plot_caco2)
dev.off()

cat("Caco-2 Heatmap saved to:", heatmap_pdf_path_caco2, "\n")

# ---------------------------- #
#          End of Script        #
# ---------------------------- #
