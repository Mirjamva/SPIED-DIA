################################################################################
# PHOSPHO SPIED-DIA MANUSCRIPT
# Script to analyze and plot DIA-NN results of a full SILAC dilution series
# Prduces Figure 1D, E, Supplementary Figures 1 A, B, C
# in R 4.3.0. Designed for GitHub reproducibility:
#   - Place "report.tsv" in a "data/" folder. (Can be downloaded from the associated PRIDE repository)
#   - Outputs (plots, CSV files) will be written to an "output/" folder.
#   - Adjust the "data_dir" and "output_dir" variables below if needed.
################################################################################

#########################
# 1) PACKAGES AND SETUP #
#########################

library(tidyverse)
library(stringr)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(fst)
library(patchwork)
library(data.table)

# Relative paths (adjust if needed). The script expects:
#  - data_dir: folder containing your input files (like "report.tsv").
#  - output_dir: folder where figures/CSVs will be saved.
data_dir   <- file.path("..", "data")
output_dir <- file.path("..", "output")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###############################################
# 2) FUNCTION TO CONVERT MOD SEQ TO pS/pT/pY  #
###############################################

fun_modseq_diannMQlib2pSTY <- function(seq_p) {
  seq_p <- gsub("_", "", seq_p)
  jptseq <- gsub("S\\(Phospho \\(STY\\)\\)", "pS", seq_p)
  jptseq <- gsub("Y\\(Phospho \\(STY\\)\\)", "pY", jptseq)
  jptseq <- gsub("T\\(Phospho \\(STY\\)\\)", "pT", jptseq)
  jptseq <- gsub("S\\(UniMod:21\\)", "pS", jptseq)
  jptseq <- gsub("Y\\(UniMod:21\\)", "pY", jptseq)
  jptseq <- gsub("T\\(UniMod:21\\)", "pT", jptseq)
  jptseq <- gsub("\\(SILAC-R-H\\)|\\(SILAC-K-H\\)|\\(SILAC-R-L\\)|\\(SILAC-K-L\\)", "", jptseq)
  jptseq <- gsub("\\(Oxidation \\(M\\)\\)|\\(Acetyl \\(Protein N-term\\)\\)", "", jptseq)
  jptseq <- gsub("\\(UniMod:\\d+\\)", "", jptseq)
  jptseq
}

###########################################################
# 3) READ AND ORGANIZE INPUT DATA (report.tsv + df_exprat)#
###########################################################

report_path <- file.path(data_dir, "report.tsv")
reporttsv   <- fread(report_path, sep = "\t")  # Fast reading with data.table

# Example expression ratio table
df_exprat <- data_frame(
  DIL = c("1_1","1_3","1_7","1_15","1_31","1_63","1_127","1_255",
          "1_511","1_1023","1_2047"),
  ExpRat = c(0, log2(3), log2(7), log2(15), log2(31), log2(63),
             log2(127), log2(255), log2(511), log2(1023), log2(2047)),
  lowerlim = c(0.6,1.6,2.6,3.6,4.6,5.6,6.6,7.6,8.6,9.6,10.6),
  upperlim = c(1.4,2.4,3.4,4.4,5.4,6.4,7.4,8.4,9.4,10.4,11.4)
)

# Convert to data.table for speed
setDT(reporttsv)

# Clean up certain columns (e.g. remove substrings from 'Run', unify 'Modified.Sequence', etc.)
reporttsv[, run_simple        := gsub("Popeye_20221125_MVB_HStdia_SIL_", "", Run)]
reporttsv[, run_simple        := gsub("_[A-Z0-9]*_[0-9]_[0-9]*$", "", run_simple)]
reporttsv[, run_simple        := gsub("100ng_*", "", run_simple)]
reporttsv[, DIL               := gsub("_[123]$", "", run_simple)]
reporttsv[, SILAClabel        := str_extract(Modified.Sequence, "-H|-L")]
reporttsv[, SILAClabel        := gsub("-", "", SILAClabel)]
reporttsv[, Modified.Sequence := gsub("\\(UniMod:21\\)", "(Phospho (STY))", Modified.Sequence)]

# Merge expression ratio info (df_exprat) based on "DIL"
report_organised <- merge(reporttsv, df_exprat, by = "DIL", all.x = TRUE)

# Filter for phospho-peptides with Q < 0.05
report_organised %>%
  filter(Channel.Q.Value < 0.05, grepl("Phospho", Modified.Sequence)) ->
  report_CQVphosfil

# Factor for DIL levels
DIL_levels <- c("1_1","1_3","1_7","1_15","1_31","1_63","1_127","1_255")
dil_labels <- c("1:1","1:3","1:7","1:15","1:31","1:63","1:127","1:255")

report_CQVphosfil$DIL <- factor(report_CQVphosfil$DIL, levels = DIL_levels)
report_organised$DIL  <- factor(report_organised$DIL,  levels = DIL_levels)

################################################################################
# 4) SUPPLEMENTARY FIGURE 1A - CORRELATIONS
################################################################################

# Summarize Ms1.Area by precursor & DIL
tmp_normalapproach <- report_CQVphosfil %>%
  group_by(Precursor.Id, DIL, SILAClabel) %>%
  filter(n() > 1) %>%
  summarise(Ms1.Area = median(Ms1.Area)) %>%
  mutate(DIL = gsub("_", ":", DIL))

tmp_normalapproach$DIL <- factor(tmp_normalapproach$DIL, levels = dil_labels)
tmp_normalapproach %>% filter(!is.na(DIL)) -> tmp_normalapproach

tmp_normalapproach %>%
  mutate(Precursor.Id = gsub("-L|-H", "", Precursor.Id)) ->
  tmp_normalapproach_2join

# spike in approach
report_CQVphosfil%>% 
  filter(SILAClabel == "L") %>% 
  mutate(DIL = gsub("_", ":", DIL)) %>%
  group_by(Precursor.Id, DIL) %>%
  filter(n() > 1) %>%
  summarise(Channel.L = median(Channel.L, na.rm = T), 
            Channel.H = median(Channel.H, na.rm = T)) %>%
  pivot_longer(
    cols = c(Channel.L, Channel.H),
    names_to = "SILAClabel.MS1",
    values_to = "MS1.Intensity"
  ) -> plot_spikein

# Join with the normal approach
# and keep only SILAClabel == "H" with Ms1.Area not NA
df_corrplot <- plot_spikein %>%
  mutate(
    SILAClabel   = gsub("Channel.", "", SILAClabel.MS1),
    Precursor.Id = gsub("-L|-H", "", Precursor.Id)
  ) %>%
  select(-SILAClabel.MS1) %>%
  left_join(tmp_normalapproach_2join, by = c("Precursor.Id","DIL","SILAClabel")) %>%
  filter(SILAClabel == "H", !is.na(Ms1.Area)) %>%
  mutate(DIL = gsub("_", ":", DIL))

df_corrplot$DIL <- factor(df_corrplot$DIL, levels = c("1:1","1:3","1:7","1:15","1:31","1:63","1:127"))

calc_r_squared <- function(data) {
  valid_data <- data[
    !is.na(data$MS1.Intensity) & !is.nan(data$MS1.Intensity) & !is.infinite(data$MS1.Intensity) &
      !is.na(data$Ms1.Area)      & !is.nan(data$Ms1.Area)      & !is.infinite(data$Ms1.Area)      &
      data$MS1.Intensity > 0     & data$Ms1.Area > 0,
  ]
  if (
    nrow(valid_data) > 1 &&
    length(unique(valid_data$MS1.Intensity)) > 1 &&
    length(unique(valid_data$Ms1.Area)) > 1
  ) {
    model <- lm(log10(MS1.Intensity) ~ log10(Ms1.Area), data = valid_data)
    summary(model)$r.squared
  } else {
    NA
  }
}

r_squared_df <- df_corrplot %>%
  group_by(DIL) %>%
  summarize(
    R_squared = calc_r_squared(cur_data()),
    n         = n(),
    .groups   = 'drop'
  ) %>%
  mutate(
    label   = paste0("R^2 = ", round(R_squared, 3)),
    n_label = paste0("n = ", n),
    x_pos   = 1,
    y_pos   = 5.5,
    y_pos_n = 5
  )

publication_ready_plot <- ggplot(df_corrplot, aes(x = log10(Ms1.Area), y = log10(MS1.Intensity))) +
  geom_point(alpha = 0.05) +
  facet_grid(cols = vars(DIL)) +
  theme_minimal() +
  labs(x = "H intensity (Normal)", y = "H intensity (Spike-in)") +
  geom_text(data = r_squared_df,
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 0, color = "blue", size = 3) +
  geom_text(data = r_squared_df,
            aes(x = x_pos, y = y_pos_n, label = n_label),
            hjust = 0, vjust = 0, color = "blue", size = 3)

ggsave(
  file.path(output_dir, "2025_Corr_plot_PhosOnly.jpeg"),
  publication_ready_plot, width = 11, height = 2, units = "in"
)

write.csv(
  df_corrplot,
  file.path(output_dir, "SourceData_SuppFigure1A.csv")
)

################################################################################
# 5) FIGURE 1D - PHOSPHOSITE-LEVEL QUANTIFICATION (NORMAL VS. SPIED)
################################################################################

PTMSITECONF <- 0.5

plot_SILACratios_normalfilter <- report_organised %>%
  filter(Channel.Q.Value < 0.05,
         PTM.Q.Value < 0.05,
         grepl("Phospho", Modified.Sequence),
         DIL %in% DIL_levels,
         PTM.Site.Confidence > PTMSITECONF) %>%
  mutate(phosseq = fun_modseq_diannMQlib2pSTY(Modified.Sequence)) %>%
  group_by(phosseq, DIL, SILAClabel, ExpRat, lowerlim, upperlim) %>%
  filter(n_distinct(run_simple) > 1) %>%
  ungroup() %>%
  mutate(
    Modified.Sequence = gsub("-K|-R|-L|-H", "", Modified.Sequence),
    Precursor.Id      = gsub("-K|-R|-L|-H", "", Precursor.Id)
  ) %>%
  pivot_wider(
    names_from   = SILAClabel,
    values_from  = Ms1.Area,
    names_prefix = "Ms1.",
    id_cols      = c(Modified.Sequence, Precursor.Id, run_simple, DIL, ExpRat, lowerlim, upperlim, phosseq)
  ) %>%
  mutate(SILACratio = log2(Ms1.L / Ms1.H), method = "Normal DIA-SILAC") %>%
  filter(is.finite(SILACratio), !is.na(SILACratio), SILACratio != Inf, SILACratio != -Inf) %>%
  select(Modified.Sequence, Precursor.Id, DIL, ExpRat, lowerlim, upperlim, SILACratio, method, phosseq) %>%
  group_by(phosseq, DIL, ExpRat, lowerlim, upperlim, method) %>%
  summarize(SILACratio = median(SILACratio), .groups = "drop")

H_ids <- report_organised %>%
  filter(SILAClabel == "H", Channel.Q.Value < 0.5) %>%
  mutate(run_precursor_id = paste0(Run, gsub("-H\\)", "-L\\)", Precursor.Id))) %>%
  pull(run_precursor_id)

plot_SILACratios_spikein <- report_organised %>%
  mutate(phosseq = fun_modseq_diannMQlib2pSTY(Modified.Sequence)) %>%
  filter(Channel.Q.Value < 0.05,
         PTM.Q.Value < 0.05,
         SILAClabel == "L",
         Channel.L > 1000,
         DIL %in% DIL_levels,
         PTM.Site.Confidence > PTMSITECONF) %>%
  mutate(
    run_precursor_id = paste0(Run, Precursor.Id),
    H_idd            = run_precursor_id %in% H_ids
  ) %>%
  group_by(phosseq, DIL, SILAClabel, ExpRat, lowerlim, upperlim) %>%
  filter(n_distinct(run_simple) > 1) %>%
  group_by(phosseq, DIL) %>%
  filter(sum(H_idd) > 1) %>%
  ungroup() %>%
  mutate(
    SILACratio = log2(Channel.L / Channel.H),
    method     = "SPIED-DIA"
  ) %>%
  filter(is.finite(SILACratio), !is.na(SILACratio), SILACratio != Inf, SILACratio != -Inf) %>%
  mutate(
    int_4bin = as.character(
      cut(Channel.L, breaks = quantile(Channel.L, probs = 0:4/4, na.rm = TRUE, type = 5),
          include.lowest = TRUE)
    ),
    DIL = factor(DIL, levels = DIL_levels)
  ) -> plotSILACratios_spikein_forLOQ
plotSILACratios_spikein_forLOQ %>%
  rowwise() %>%
  select(Modified.Sequence, Precursor.Id, DIL, ExpRat, lowerlim, upperlim, SILACratio, phosseq) %>%
  mutate(method = "DIA-SPIED") %>%
  group_by(DIL, ExpRat, lowerlim, upperlim, method, phosseq) %>%
  summarize(SILACratio = median(SILACratio), .groups = "drop")

plot_SILAC_comb <- rbind(plot_SILACratios_normalfilter, plot_SILACratios_spikein)
plot_SILAC_comb$DIL <- factor(plot_SILAC_comb$DIL, levels = DIL_levels)

bar_data_normal <- plot_SILAC_comb %>%
  group_by(DIL, method) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  complete(DIL, method, fill = list(count = 0))

custom_colors <- c("orange","blue")
custom_labels <- c("SPIED-DIA","Normal SILAC-DIA")

violin_plot_normal <- ggplot(plot_SILAC_comb, aes(x = DIL, y = SILACratio, fill = method)) +
  geom_boxplot(
    size = 0.2, color = "black", outlier.alpha = 0,
    position = position_dodge(preserve = "single")
  ) +
  geom_errorbarh(
    aes(x = DIL, y = ExpRat, xmin = lowerlim, xmax = upperlim),
    inherit.aes = FALSE, alpha = 0.5, linetype = 11, colour = "black", size = 0.25
  ) +
  labs(y = "log2(SILACratio)", x = "none") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  ylim(-2, 10) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  guides(fill = FALSE)

bar_plot_normal <- ggplot(plot_SILAC_comb, aes(x = DIL, fill = method)) +
  geom_bar(position = position_dodge(width = 0.8), alpha = 0.7, color = "black", size = 0.2) +
  geom_text_repel(
    stat = 'count', aes(label = ..count.., group = method),
    position = position_dodge(width = 0.8),
    vjust = -0.25, size = 3, angle = 90
  ) +
  labs(x = "Dilution Factor", y = "Identified p-sites") +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.title    = element_blank(),
    legend.position = "bottom",
    legend.direction= "vertical"
  ) +
  scale_fill_manual(values = custom_colors, labels = custom_labels) +
  scale_x_discrete(limits = DIL_levels, labels = dil_labels)

final_plot <- violin_plot_normal + bar_plot_normal +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, paste0("2025Comparison_fig3_PTMSC", PTMSITECONF, "psite.pdf")),
  final_plot, width = 3, height = 5.5, units = "in", dpi = 300
)

write.csv(
  plot_SILAC_comb,
  file.path(output_dir, "SourceData_Figure1D.csv")
)

############################
# 6) FIGURE 1E - NOISE, LOQ
############################

simple_output <- plotSILACratios_spikein_forLOQ %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    H_intensity = log10(mean(Channel.H, na.rm = TRUE)),
    L_intensity = log10(mean(Channel.L, na.rm = TRUE))
  ) %>%
  select(DIL, SILAClabel, Modified.Sequence, SILACratio, H_intensity, L_intensity, Run) %>%
  group_by(DIL, SILAClabel) %>%
  mutate(
    int_10bin   = as.character(cut_number(L_intensity, 5, labels = FALSE)),
    H_int_10bin = as.character(cut_number(H_intensity, 10, labels = FALSE))
  ) %>%
  group_by(SILAClabel) %>%
  mutate(
    H_int_10bin_total = as.character(cut_interval(H_intensity, 20, labels = FALSE)),
    L_int_10bin_total = as.character(cut_interval(L_intensity, 20, labels = FALSE))
  ) %>%
  ungroup()

simple_output$DIL <- factor(simple_output$DIL, levels = DIL_levels)
df_exprat$DIL     <- factor(df_exprat$DIL, levels = DIL_levels)
simple_output$int_10bin <- factor(simple_output$int_10bin, levels = c("1","2","3","4","5"))

plot_tmp <- simple_output %>%
  group_by(DIL, int_10bin) %>%
  summarize(
    median_SILACratio  = median(SILACratio, na.rm = TRUE),
    sd_SILACratiobin   = sd(SILACratio),
    .groups            = "drop"
  ) %>%
  left_join(df_exprat[, c("DIL", "ExpRat")], by = "DIL")

plot_tmp$DIL <- factor(plot_tmp$DIL, levels = DIL_levels)

high_contrast_colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

dotline_10bin <- ggplot(plot_tmp, aes(x = DIL, y = median_SILACratio, group = as.factor(int_10bin),
                                      color = as.factor(int_10bin))) +
  scale_color_manual(values = high_contrast_colors) +
  geom_point() +
  geom_line() +
  geom_line(aes(y = ExpRat), colour = "black", linetype = "dotted") +
  labs(y = "Median SILACratio (log2)") +
  theme_minimal() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())

hist_10bin <- ggplot(filter(simple_output)) +
  geom_histogram(aes(x = L_intensity, fill = as.factor(int_10bin)), bins = 100) +
  coord_flip() +
  labs(x = "L Intensity (log10)", y = "# Precursors", fill = "Ref. Intensity Bin") +
  theme_minimal() +
  scale_fill_manual(values = high_contrast_colors)

plot_CVs <- simple_output %>%
  ungroup() %>%
  group_by(Modified.Sequence, DIL) %>%
  summarize(
    mean_intbin   = mean(as.numeric(int_10bin)),
    sd_SILACratio = sd(SILACratio),
    .groups       = "drop"
  ) %>%
  group_by(DIL) %>%
  mutate(int_10bin = as.character(cut_interval(mean_intbin, 5, labels = FALSE)))

plot_CVs$DIL       <- factor(plot_CVs$DIL,       levels = DIL_levels)
plot_CVs$int_10bin <- factor(plot_CVs$int_10bin, levels = c("1","2","3","4","5"))

sd_plot <- ggplot(plot_CVs) +
  geom_boxplot(aes(x = DIL, y = sd_SILACratio, fill = as.factor(int_10bin)),
               outlier.alpha = 0, size = 0.2) +
  ylim(0, 2.5) +
  theme_minimal() +
  scale_fill_manual(values = high_contrast_colors) +
  theme(legend.position = "none") +
  scale_x_discrete(limits = DIL_levels, labels = dil_labels) +
  labs(x = "Dilution", y = "sd log2(SILACratio)")

final_plot_noise <- dotline_10bin + hist_10bin + sd_plot +
  plot_layout(guides = "collect", ncol = 2, widths = c(3,1), heights = c(3,2)) &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "Figure1E_NoisePhospho.jpeg"),
  final_plot_noise, width = 6.5, height = 5, units = "in", dpi = 300
)

write.csv(
  plot_tmp,
  file.path(output_dir, "SourceData_Figure1E_Linegraph.csv")
)
write.csv(
  simple_output,
  file.path(output_dir, "SourceData_Figure1E_Boxplot.csv")
)
write.csv(
  plot_CVs,
  file.path(output_dir, "SourceData_Figure1E_Hist.csv")
)

################################################################################
# 7) SUPPLEMENTARY FIGURE 1B - FILTER BY PTM.SITE.CONFIDENCE
################################################################################

report_DT <- setDT(report_organised)

plotSILACratios_spikein_noPTMSC <- report_DT[
  Channel.Q.Value < 0.05 & PTM.Q.Value < 0.05 &
    SILAClabel == "L" & Channel.L > 1000 & DIL %in% DIL_levels
][
  , run_precursor_id := paste0(Run, Precursor.Id)
][
  , H_idd := run_precursor_id %in% H_ids
][
  , .SD[.N > 1 & sum(H_idd) > 1], by = .(Precursor.Id, DIL)
][
  , SILACratio := log2(Channel.L / Channel.H)
][
  is.finite(SILACratio) & !is.na(SILACratio) & SILACratio != Inf & SILACratio != -Inf
][
  , int_4bin := as.character(
    cut(Channel.L, breaks = quantile(Channel.L, probs = 0:4/4, na.rm = TRUE, type = 5),
        include.lowest = TRUE)
  )
][
  , DIL := factor(DIL, levels = DIL_levels)
]

setDF(plotSILACratios_spikein_noPTMSC)

total_results <- data.frame(
  PTM.Q.Value_filter = numeric(),
  R2_alldil_psite    = numeric(),
  n_alldil_psite     = numeric(),
  R2_4dil_psite      = numeric(),
  n_4dil_psite       = numeric(),
  R2_alldil_isobaric = numeric(),
  n_alldil_isobaric  = numeric(),
  R2_4dil_isobaric   = numeric(),
  n_4dil_isobaric    = numeric(),
  stringsAsFactors   = FALSE
)

report_CQVphosfil %>%
  filter(SILAClabel == "L") %>%
  mutate(jptseq = fun_modseq_diannMQlib2pSTY(Modified.Sequence)) %>%
  mutate(
    nophosseq   = gsub("\\(Phospho \\(STY\\)\\)", "", Modified.Sequence),
    n_phosgroups= str_count(Modified.Sequence, "\\(Phospho \\(STY\\)\\)"),
    phos_sites  = str_extract_all(Modified.Sequence, "\\(Phospho \\(STY\\)\\)")
  ) %>%
  group_by(Run, Stripped.Sequence, n_phosgroups) %>%
  filter(n_distinct(jptseq) > 1) %>%
  ungroup() %>%
  mutate(key = paste0(Precursor.Id, Run)) ->
  phospho_localisation

PTMQVAL_filters <- c(0.999, 0.99, 0.95, 0.9, 0.75, 0.7, 0.6, 0.5, 0.49, 0.45, 0.25, 0.1, 0)

for (PTMQVAL in PTMQVAL_filters) {
  simple_output <- plotSILACratios_spikein_noPTMSC %>%
    filter(PTM.Site.Confidence > PTMQVAL, SILAClabel == "L") %>%
    mutate(SILACratio = log2(Channel.L / Channel.H)) %>%
    filter(is.finite(SILACratio), !is.na(SILACratio))
  
  simple_output_pepseq <- simple_output %>%
    group_by(Precursor.Id, Run, DIL, ExpRat) %>%
    summarise(SILACratio = median(SILACratio), .groups = "drop")
  
  model            <- lm(ExpRat ~ SILACratio, data = simple_output_pepseq)
  R2_alldil_psite  <- summary(model)$r.squared
  n_alldil_psite   <- nrow(simple_output_pepseq)
  
  simple_output_4DIL <- simple_output_pepseq %>%
    filter(DIL %in% c("1_1","1_3","1_7","1_15"))
  
  model        <- lm(ExpRat ~ SILACratio, data = simple_output_4DIL)
  R2_4dil_psite<- summary(model)$r.squared
  n_4dil_psite <- nrow(simple_output_4DIL)
  
  simple_output_iso <- simple_output %>%
    mutate(key = paste0(Precursor.Id, Run)) %>%
    filter(key %in% phospho_localisation$key)
  
  model              <- lm(ExpRat ~ SILACratio, data = simple_output_iso)
  R2_alldil_isobaric <- summary(model)$r.squared
  n_alldil_isobaric  <- nrow(simple_output_iso)
  
  simple_output_4DIL <- simple_output_iso %>%
    filter(DIL %in% c("1_1","1_3","1_7","1_15"))
  
  model               <- lm(ExpRat ~ SILACratio, data = simple_output_4DIL)
  R2_4dil_isobaric    <- summary(model)$r.squared
  n_4dil_isobaric     <- nrow(simple_output_4DIL)
  
  total_results <- rbind(total_results, data.frame(
    PTM.Q.Value_filter   = PTMQVAL,
    R2_alldil_psite      = R2_alldil_psite,
    n_alldil_psite       = n_alldil_psite,
    R2_4dil_psite        = R2_4dil_psite,
    n_4dil_psite         = n_4dil_psite,
    R2_alldil_isobaric   = R2_alldil_isobaric,
    n_alldil_isobaric    = n_alldil_isobaric,
    R2_4dil_isobaric     = R2_4dil_isobaric,
    n_4dil_isobaric      = n_4dil_isobaric
  ))
}

total_results_long <- total_results %>%
  pivot_longer(
    cols = starts_with("R2_") | starts_with("n_"),
    names_to = c("metric", "dilution", "type"),
    names_pattern = "(R2|n)_(alldil|4dil)_(psite|isobaric)",
    values_to = "value"
  ) %>%
  mutate(
    metric   = ifelse(metric == "R2", "R2", "n"),
    dilution = ifelse(dilution == "alldil", "alldil", "4dil"),
    type     = ifelse(type == "psite", "All", "isobaric")
  )

write.csv(
  total_results_long,
  file.path(output_dir, "SourceData_SuppFigure1B.csv")
)

qval_num_simple <- c(0.999,0.9,0.75,0.7,0.6,0.5,0.49,0.25,0.1,0)

total_results_r2 <- filter(total_results_long, metric == "R2")
R2_total_phosseq <- ggplot(total_results_r2, aes(x = PTM.Q.Value_filter, y = value,
                                                 color = dilution, linetype = type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("alldil" = "red","4dil" = "blue")) +
  scale_linetype_manual(values = c("isobaric"="dotdash","All"="solid")) +
  scale_x_continuous(breaks = qval_num_simple) +
  labs(y = "R2 Value") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")

total_results_r2 <- filter(total_results_long, metric == "n")
n_total_phosseq  <- ggplot(total_results_r2, aes(x = PTM.Q.Value_filter, y = value,
                                                 color = dilution, linetype = type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("alldil"="red","4dil"="blue")) +
  scale_linetype_manual(values = c("isobaric"="dotdash","All"="solid")) +
  scale_x_continuous(breaks = qval_num_simple) +
  labs(y = "# quantified", x = "PTM.Site.Confidence > x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

final_plot_Supp1B <- R2_total_phosseq + n_total_phosseq +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "SuppFigure1B_FilterPTMSiteConfidence.pdf"),
  final_plot_Supp1B, width = 7, height = 5, units = "in", dpi = 300
)

################################################################################
# 8) SUPPLEMENTARY FIGURE 1C - ISOBARIC PEPTIDES (PHOSPHO LOCALISATION)
################################################################################

dt <- report_CQVphosfil[SILAClabel == "L" & Channel.Q.Value < 0.05]
dt[, jptseq       := fun_modseq_diannMQlib2pSTY(Modified.Sequence)]
dt[, nophosseq    := gsub("\\(Phospho \\(STY\\)\\)", "", Modified.Sequence)]
dt[, n_phosgroups := str_count(Modified.Sequence, "\\(Phospho \\(STY\\)\\)")]
dt[, phos_sites   := str_extract_all(Modified.Sequence, "\\(Phospho \\(STY\\)\\)")]

group_sum   <- dt[, .(n_jptseq = uniqueN(jptseq)), by = .(Run, Stripped.Sequence, n_phosgroups)]
groups_keep <- group_sum[n_jptseq > 1]
phospho_localisation <- dt[groups_keep, on = .(Run, Stripped.Sequence, n_phosgroups)]
phospho_localisation[, key := paste0(Precursor.Id, Run)]
phospho_localisation <- as.data.frame(phospho_localisation)

phospho_localisation_switch <- phospho_localisation %>%
  filter(DIL %in% c("1_1","1_15")) %>%
  group_by(nophosseq, n_phosgroups, Run, Precursor.Charge) %>%
  mutate(
    paired_jptseq = map_chr(jptseq, ~ {
      same_charge_jptseqs <- jptseq[Precursor.Charge == first(Precursor.Charge)]
      other_jptseqs       <- setdiff(same_charge_jptseqs, .x)
      if (length(other_jptseqs) > 0) sample(other_jptseqs, 1) else .x
    })
  ) %>%
  ungroup()

mapping_table <- phospho_localisation_switch %>%
  group_by(nophosseq, n_phosgroups, Run, Precursor.Charge, paired_jptseq) %>%
  summarise(
    Channel.H.original = first(Channel.H),
    IM.original        = first(IM),
    RT.original        = first(RT),
    .groups            = "drop"
  )

phospho_localisation_switch <- phospho_localisation_switch %>%
  left_join(mapping_table,
            by=c("nophosseq","n_phosgroups","Run","Precursor.Charge","jptseq"="paired_jptseq")) %>%
  mutate(
    Channel.H.FALSE = ifelse(jptseq == paired_jptseq, Channel.H, Channel.H.original),
    IM.FALSE        = ifelse(jptseq == paired_jptseq, IM, IM.original),
    RT.FALSE        = ifelse(jptseq == paired_jptseq, RT, RT.original)
  ) %>%
  select(Modified.Sequence, nophosseq, n_phosgroups, Run, run_simple, DIL,
         Precursor.Charge, jptseq, paired_jptseq, Channel.L, Channel.H, Channel.H.FALSE,
         IM, RT, IM.FALSE, RT.FALSE) %>%
  rowwise() %>%
  filter(Channel.H != Channel.H.FALSE, IM != IM.FALSE, RT != RT.FALSE)

phospho_localisation_switch_long <- phospho_localisation_switch %>%
  mutate(
    SILAC_Ratio_True  = Channel.L / Channel.H,
    SILAC_Ratio_False = Channel.L / Channel.H.FALSE
  ) %>%
  pivot_longer(
    cols = c("SILAC_Ratio_True","SILAC_Ratio_False"),
    names_to   = "Pair_Type",
    values_to  = "SILAC_Ratio"
  ) %>%
  mutate(Pair_Type = ifelse(Pair_Type == "SILAC_Ratio_True","TRUE","FALSE")) %>%
  group_by(DIL) %>%
  mutate(n_modseq = paste0("n = ", n_distinct(Modified.Sequence))) %>%
  ungroup()

line_heights        <- data.frame(DIL = c("1_1","1_15"), y_intercept = c(0, 4))
common_xbreaks_dens <- c(0, 0.3, 0.6, 0.9, 1.2)

scatter_plot <- ggplot(phospho_localisation_switch_long,
                       aes(x=log10(Channel.L), y=log2(SILAC_Ratio), colour=Pair_Type)) +
  geom_label(
    mapping = aes(x=-Inf, y=Inf, label=n_modseq),
    color = "black", hjust=-0.05, vjust=1.05, size=2,
    label.padding=unit(0.15,"lines"), label.size=unit(0.15,"lines")
  ) +
  geom_hline(data=line_heights, aes(yintercept=y_intercept),
             linetype="dashed", color="blue") +
  geom_point(alpha=0.05) +
  scale_colour_manual(values=c("TRUE"="green","FALSE"="red2")) +
  facet_grid(rows=vars(DIL)) +
  theme_minimal() +
  guides(colour=FALSE) +
  labs(y="L/H (log2)", x="L Intensity (log10)", subtitle="PTM.Site.Confidence > 0") +
  xlim(2,6) +
  theme(strip.text.y=element_blank())

density_plot <- ggplot(phospho_localisation_switch_long, aes(y=log2(SILAC_Ratio), fill=Pair_Type)) +
  geom_hline(data=line_heights, aes(yintercept=y_intercept),
             linetype="dashed", color="blue") +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("TRUE"="green","FALSE"="red2")) +
  facet_grid(rows=vars(DIL)) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
  labs(fill="Positional isomer H/L pairs") +
  scale_x_continuous(breaks=common_xbreaks_dens)

combined_plot <- scatter_plot + density_plot +
  plot_layout(guides="collect", widths=c(5,2)) &
  theme(legend.position="bottom")

# Additional filter for PTM.Site.Confidence > 0.5 for demonstration
phospho_localisation_switch_long_PTM0 <- phospho_localisation_switch_long

phospho_localisation_switch <- phospho_localisation %>%
  filter(DIL %in% c("1_1","1_15"), PTM.Site.Confidence > 0.5, PTM.Q.Value < 0.05) %>%
  mutate(jptseq=fun_modseq_diannMQlib2pSTY(Modified.Sequence)) %>%
  group_by(nophosseq, n_phosgroups, Run, Precursor.Charge) %>%
  mutate(
    paired_jptseq=map_chr(jptseq, ~ {
      same_charge_jptseqs <- jptseq[Precursor.Charge == first(Precursor.Charge)]
      other_jptseqs       <- setdiff(same_charge_jptseqs, .x)
      if (length(other_jptseqs)>0) sample(other_jptseqs,1) else .x
    })
  ) %>%
  ungroup()

mapping_table <- phospho_localisation_switch %>%
  group_by(nophosseq, n_phosgroups, Run, Precursor.Charge, paired_jptseq) %>%
  summarise(
    Channel.H.original=first(Channel.H),
    IM.original       =first(IM),
    RT.original       =first(RT),
    .groups="drop"
  )

phospho_localisation_switch <- phospho_localisation_switch %>%
  left_join(mapping_table,
            by=c("nophosseq","n_phosgroups","Run","Precursor.Charge","jptseq"="paired_jptseq")) %>%
  mutate(
    Channel.H.FALSE=ifelse(jptseq==paired_jptseq,Channel.H,Channel.H.original),
    IM.FALSE=ifelse(jptseq==paired_jptseq,IM,IM.original),
    RT.FALSE=ifelse(jptseq==paired_jptseq,RT,RT.original)
  ) %>%
  select(Modified.Sequence,nophosseq,n_phosgroups,Run,run_simple,DIL,
         Precursor.Charge,jptseq,paired_jptseq,Channel.L,Channel.H,Channel.H.FALSE,
         IM,RT,IM.FALSE,RT.FALSE) %>%
  rowwise() %>%
  filter(Channel.H!=Channel.H.FALSE, IM!=IM.FALSE, RT!=RT.FALSE)

phospho_localisation_switch_long <- phospho_localisation_switch %>%
  mutate(
    SILAC_Ratio_True = Channel.L / Channel.H,
    SILAC_Ratio_False= Channel.L / Channel.H.FALSE
  ) %>%
  pivot_longer(
    cols=c("SILAC_Ratio_True","SILAC_Ratio_False"),
    names_to="Pair_Type",
    values_to="SILAC_Ratio"
  ) %>%
  mutate(Pair_Type=ifelse(Pair_Type=="SILAC_Ratio_True","TRUE","FALSE")) %>%
  group_by(DIL) %>%
  mutate(n_modseq=paste0("n = ", n_distinct(Modified.Sequence))) %>%
  ungroup()

scatter_plot_fil <- ggplot(phospho_localisation_switch_long,
                           aes(x=log10(Channel.L), y=log2(SILAC_Ratio), colour=Pair_Type)) +
  geom_label(
    mapping=aes(x=-Inf, y=Inf, label=n_modseq),
    color="black", hjust=-0.05, vjust=1.05, size=2,
    label.padding=unit(0.15,"lines"), label.size=unit(0.15,"lines")
  ) +
  geom_hline(data=line_heights, aes(yintercept=y_intercept),
             linetype="dashed", color="blue") +
  geom_point(alpha=0.05) +
  scale_colour_manual(values=c("TRUE"="green","FALSE"="red2")) +
  facet_grid(rows=vars(DIL)) +
  theme_minimal() +
  guides(colour=FALSE) +
  labs(y="L/H (log2)", x="L Intensity (log10)", subtitle="PTM.Site.Confidence > 0.5") +
  xlim(2,6) +
  theme(strip.text.y=element_blank())

density_plot_fil <- ggplot(phospho_localisation_switch_long, aes(y=log2(SILAC_Ratio), fill=Pair_Type)) +
  geom_hline(data=line_heights, aes(yintercept=y_intercept), linetype="dashed", color="blue") +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("TRUE"="green","FALSE"="red2")) +
  facet_grid(rows=vars(DIL)) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
  labs(fill="Positional isomer H/L pairs") +
  scale_x_continuous(breaks=common_xbreaks_dens)

combined_plot_1C <- scatter_plot + density_plot + scatter_plot_fil + density_plot_fil +
  plot_layout(guides="collect", widths=c(5,2)) &
  theme(legend.position="bottom")

ggsave(
  file.path(output_dir, "SuppFigure1C_isobaricPairs_RatioScatter_comb.pdf"),
  combined_plot_1C, width=5, height=5, units="in", dpi=300
)

phospho_localisation_switch_long_PTM0 %>%
  mutate(PTMfil="0") -> PTMswitch1
phospho_localisation_switch_long %>%
  mutate(PTMfil="0.5") -> PTMswitch2

write.csv(
  rbind(PTMswitch1, PTMswitch2),
  file.path(output_dir, "SourceData_SuppFigure1C.csv")
)

################################################################################
# End of Script
################################################################################
