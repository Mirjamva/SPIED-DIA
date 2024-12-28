################################################################################
# Script Name:     BM3_Ecoli_Background_Analysis.R
# Author:          Mirjam van Bentum
# Date:            2024-12-12
#
# Description:
# This R script performs data analysis for a dilution series (BM3) experiment set in an E. coli 
# background. It processes phospho-peptide data, applies sequence modifications, 
# and generates multiple visualizations (boxplots, jitter plots, and histograms) 
# that summarize the identification and quantification of peptides under 
# different conditions. Produces Figure 1I, and supplementary figure 3/4
#
# Dependencies:
# - R version >= 4.0
# - Packages: data.table, tidyverse, patchwork, viridis
#
# Input Data:
# - "20241104JPTlib_noNA.csv": A CSV file containing raw JPT library data.
# - "report.tsv": A TSV file containing DIA-NN report data.
#
# Output:
# - PDF plots saved to the "Plots/" directory.
#
# To run:
# 1. Place input files "20241104JPTlib_noNA.csv", "sky_output.csv", and 
#    "report.tsv" in the working directory or specify their paths below.
# 2. Ensure the "Plots/" directory exists or create it before running.
# 3. Run the script in R.
################################################################################


# load libraries
library(data.table)  # Efficient data handling and processing
library(tidyverse)   # Collection of R packages for data science
library(patchwork)   # Combining ggplot2 plots
library(viridis)     # Color palettes for visualization

# Set paths for input and output -----------------------------------------------
# Adjust these paths as needed.
path_plot_output <- "Plots"
diann_input_file <- "20241104JPTlib_noNA.csv"
jpt_report_file  <- "report.tsv"

# Define functions -------------------------------------------------------------
# Function to modify peptide sequences from DIA-NN/MQ library to pSTY format
fun_modseq_diannMQlib2pSTY <- function(seq_p) {
  seq_p <- gsub("_", "", seq_p)
  jptseq <- seq_p
  modifications <- list(
    "S\\(Phospho \\(STY\\)\\)" = "pS",
    "Y\\(Phospho \\(STY\\)\\)" = "pY",
    "T\\(Phospho \\(STY\\)\\)" = "pT",
    "S\\(UniMod:21\\)" = "pS",
    "Y\\(UniMod:21\\)" = "pY",
    "T\\(UniMod:21\\)" = "pT",
    "\\(SILAC-R-H\\)" = "",
    "\\(SILAC-K-H\\)" = "",
    "\\(SILAC-R-L\\)" = "",
    "\\(SILAC-K-L\\)" = "",
    "\\(Oxidation \\(M\\)\\)" = "",
    "\\(Acetyl \\(Protein N-term\\)\\)" = "",
    "\\(UniMod:\\d+\\)" = ""
  )
  for (pattern in names(modifications)) {
    jptseq <- gsub(pattern, modifications[[pattern]], jptseq)
  }
  return(jptseq)
}

# Function to modify sequences to a simplified precursor format
fun_modseq_MQlib2simpleprec <- function(seq_p) {
  # Define patterns and replacements in a named vector
  subs <- c(
    "_" = "",
    "S\\(Phospho \\(STY\\)\\)" = "pS", "Y\\(Phospho \\(STY\\)\\)" = "pY", "T\\(Phospho \\(STY\\)\\)" = "pT",
    "S\\(UniMod:21\\)" = "pS", "Y\\(UniMod:21\\)" = "pY", "T\\(UniMod:21\\)" = "pT",
    "\\(SILAC-R-H\\)" = "", "\\(SILAC-K-H\\)" = "", "\\(SILAC-R-L\\)" = "", "\\(SILAC-K-L\\)" = "",
    "\\(Oxidation \\(M\\)\\)" = "", "\\(Acetyl \\(Protein N-term\\)\\)" = "",
    "\\(UniMod:1\\)" = "(ac)", "\\(UniMod:35\\)" = "(ox)"
  )
  
  # Apply all substitutions in a loop
  for (pat in names(subs)) {
    seq_p <- gsub(pat, subs[pat], seq_p)
  }
  
  seq_p
}

# Function to modify Skyline sequences to a simple precursor format
fun_modseq_sky2simpleprec <- function(seq_p){
  seq_p <- gsub("_", "", seq_p)
  jptseq <- gsub("S\\[\\+80\\]","pS", seq_p)
  jptseq <- gsub("Y\\[\\+80\\]","pY", jptseq)
  jptseq <- gsub("T\\[\\+80\\]","pT", jptseq)
  jptseq <- gsub("\\[\\+57\\]","", jptseq)
  jptseq <- gsub("^([A-Z])\\[\\+42\\]", "(ac)\\1", jptseq)
  jptseq <- gsub("\\[\\+16\\]", "(ox)", jptseq)
  return(jptseq)
}

# Data input and processing ----------------------------------------------------
diann_inputlib <- read.csv(diann_input_file)

# Extract JPT sequences from the input library
diann_inputlib %>%
  mutate(
    jptseq = gsub(pattern = "\\(SILAC-(R|K)-H\\)", replacement = "", ModifiedPeptide),
    jptseq = fun_modseq_diannMQlib2pSTY(jptseq)
  ) %>%
  pull(jptseq) %>%
  unique() -> jptseqs_inlib

sky_output <- fread(sky_output_file)

reporttsv_JPTlib <-
  fread(jpt_report_file, sep = "\t") %>%
  mutate(
    SILAClabel = ifelse(grepl("-H|-L", Modified.Sequence), 
                        str_extract(Modified.Sequence, "-H|-L"), "LFQ"),
    jptseq = fun_modseq_diannMQlib2pSTY(Modified.Sequence),
    Run_simple = gsub("Popeye_20240912_SHN_100ngEcoli_50fmolJTP_|_DIA|_[A-Z0-9]*_[0-9]_[0-9]*$", "", Run),
    sample = gsub("K_[123]", "K", Run_simple)
  )

# Parameters for filtering -----------------------------------------------------
PSC <- 0.5 # phosphosite localisation confidence
L_idd_min <- 1 #number precursors needed to be identified with meidum confidence in light

# Identify target precursors with medium confidence --------------------------------
reporttsv_JPTlib %>%
  ungroup() %>% rowwise() %>%
  filter(Channel.Q.Value < 0.5,
         SILAClabel == "-L", 
         jptseq %in% jptseqs_inlib) %>% 
  mutate(run_precursor_id = paste0(Run, gsub("-L\\)", "-H\\)",
                                             Precursor.Id))) %>%
  pull(run_precursor_id) -> L_ids_1D

# SPIED-DIA identification/quantification -----------------------------------------------------
reporttsv_JPTlib %>%
  filter(jptseq %in% jptseqs_inlib) %>% 
  filter(Channel.Q.Value < 0.05,
         PTM.Q.Value < 0.05,
         PTM.Site.Confidence > PSC,
         SILAClabel == "-H") %>% 
  mutate(run_precursor_id = paste0(Run, Precursor.Id),
         L_idd = run_precursor_id %in% L_ids_1D) %>% 
  group_by(Precursor.Id, sample) %>% 
  filter(n() > 1,
         sum(L_idd) > L_idd_min) %>% 
  group_by(Precursor.Id) %>% 
  mutate(rescalingfct = median(Ms1.Area)) %>% 
  ungroup() %>% rowwise() %>% 
  mutate(rescaled_int = log10(Channel.L/Channel.H*rescalingfct)) -> SiS2_Quant

SiS2_Quant$sample <- factor(SiS2_Quant$sample, 
                            levels = c("400ngHEK", "200ngHEK", "100ngHEK", 
                                       "50ngHEK", "25ngHEK", "0ngHEK"))

SiS2_Quant %>%
  filter(sample == "400ngHEK") %>% 
  group_by(Precursor.Id) %>% 
  filter(mean(rescaled_int) > 3) %>%
  pull(Precursor.Id) -> plot_ids

# Create plots -----------------------------------------------------------------
# Supplementary Figure 3/4
reporttsv_JPTlib %>%
  filter(Precursor.Id %in% plot_ids)%>% 
  mutate(Precursor.Id = fun_modseq_MQlib2simpleprec(Precursor.Id)) %>%
  mutate(Ratio_LH = log2(Channel.L / Channel.H)) %>%
  group_by(Precursor.Id) %>%
  mutate(log2_Ratio_LH_adjusted = Ratio_LH - mean(Ratio_LH[sample == "400ngHEK"], na.rm = TRUE)) %>% 
  ungroup() -> test_plot

order_precursorID <- test_plot %>%
  filter(sample == "400ngHEK") %>%
  group_by(Precursor.Id) %>%
  summarise(mean_Channel.H = mean(Channel.L, na.rm = TRUE)) %>%  # Ensure NA values are handled
  arrange(mean_Channel.H) %>%  # Sort by mean Channel.H
  pull(Precursor.Id)

test_plot$Precursor.Id <- factor(test_plot$Precursor.Id, levels = order_precursorID)
test_plot$sample <- factor(test_plot$sample, levels = c("400ngHEK", "200ngHEK", "100ngHEK", 
                                                        "50ngHEK", "25ngHEK", "0ngHEK"))

# Plot 1: Prepare the log10(Channel.L) plot vs precursor id
log10_plot <- test_plot %>%
  filter(sample == "400ngHEK") %>%
  ggplot(aes(y = Precursor.Id, x = log10(Channel.L))) +  # Swap x and y
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(y = "Precursor+charge",  # Adjust y-axis label for rotated axis
       x = "log10(L intensity)") +  # Adjust label for rotated axis
  theme(axis.text.y = element_text(size = 8),  # Y-axis text for readability
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),  # X-axis for log10 intensity
        plot.margin = margin(5, 5, 5, 5))  # Adjust margins

# Plot 2: Adjusted Channel.L/Channel.H ratio jitter plot
jitter_allprecursors <- ggplot(test_plot, aes(y = Precursor.Id, x = log2_Ratio_LH_adjusted, color = sample)) +  # Swap x and y
  geom_jitter(alpha = 0.8, width = 0, height = 0.2) +  # Adding jitter, swapped width/height for horizontal orientation
  theme_minimal() +
  labs(y = NULL,  # Remove y-axis label
       x = "Adjusted log2(Channel.L / Channel.H)") +  # Adjust label for rotated axis
  theme(axis.text.y = element_blank(),  # Remove y-axis labels to prevent overlap
        axis.ticks.y = element_blank(),  # Remove ticks on y-axis
        axis.text.x = element_text(angle = 70, hjust = 1),  # Rotate x-axis labels
        legend.title = element_text(size = 10),
        legend.position = "bottom") +  # Legend below the plot
  scale_color_viridis(discrete = TRUE)  # Use the viridis palette

# Combine plots
combined_plot <- log10_plot + jitter_allprecursors + 
  plot_layout(widths = c(1, 5))  

# Step 4: Save the final plot
string_plot <- paste0(path_plot_output, "spied_dia_quant.pdf")

ggsave(string_plot, 
       plot = combined_plot, width = 7, height = 8, units = "in", dpi = 300)

# Plot Figure 1I
# Summarize data for quantification --------------------------------------------
SILACratio400ng <- SiS2_Quant %>%
  filter(sample == "400ngHEK") %>% 
  group_by(Precursor.Id) %>% 
  summarise(rescaled_int = mean(rescaled_int)) %>%
  filter(rescaled_int > 3) 

# Data preparation for other samples
SiS2_Quant %>%
  filter(Precursor.Id %in% SILACratio400ng$Precursor.Id) %>%
  group_by(Precursor.Id, sample) %>%
  summarise(sd = sd(rescaled_int),
            mean = mean(rescaled_int),
            cv = sd/mean,
            rescaled_int = median(rescaled_int)) %>%
  ungroup() %>%
  right_join(SILACratio400ng, by = "Precursor.Id", suffix = c("DIL", "400ng")) %>% 
  mutate(SILACratio_norm = rescaled_intDIL-rescaled_int400ng,
         int_4bin = as.character(cut_number(rescaled_int400ng, 4, 
                                            labels = FALSE)),
         sample = gsub("HEK", "", sample),
         sample = factor(sample, levels = c("400ng", "200ng", "100ng", 
                                            "50ng", "25ng", "0ng"))) -> plot_norm_ratios

# Median value calculation per sample
plot_norm_ratios_median <- plot_norm_ratios %>%
  group_by(sample) %>%
  summarise(median_value = median(SILACratio_norm, na.rm = TRUE))

color_blind_friendly_palette <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

# Data for segments in the plot
segments_data <- data.frame(sample =  c("400ng", "200ng", "100ng", 
                                        "50ng", "25ng", "0ng"),
                            y_position = log2(c(1, 0.5, 0.25, 0.125, 0.0625, 1))) %>%
  mutate(sample = factor(sample, levels = c("400ng", "200ng", "100ng", 
                                            "50ng", "25ng", "0ng")))

# quantification boxplot
plot101 <- ggplot(data = plot_norm_ratios) +
  geom_boxplot(aes(x = sample,
                   y = log2(10^SILACratio_norm),
                   fill = int_4bin),
               linewidth = 0.2,
               color = "black",
               outlier.shape = NA) +
  geom_segment(data = segments_data,
               aes(x = as.numeric(as.factor(sample)) - 0.3, 
                   xend = as.numeric(as.factor(sample)) + 0.3,
                   y = y_position, yend = y_position),
               size = 1, linetype = "dotted", color = "darkred") +
  scale_x_discrete(limits = unique(segments_data$sample)) +
  labs(y = "Rel. Intensity (log2)",
       fill = "Intensity bin")+
  scale_y_continuous(limits = c(-5, 0.5)) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) 
plot101

#intensity distribution in 400ng condition
plot102 <- ggplot(data = filter(plot_norm_ratios, 
                                sample == "400ng"),
                  aes(x = rescaled_int400ng, fill = int_4bin)) +
  geom_histogram() +
  theme_minimal() +
  guides(fill = F)+
  theme(axis.title.y = element_blank())+
  labs(x = "L Intensity 400 ng (log10)")
plot102

# bar plot identified precursors. 
plot103 <- ggplot(filter(plot_norm_ratios, 
                         is.finite(SILACratio_norm)),
                  aes(x = sample)) +
  geom_bar(position = "dodge") + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "# Precursors")
plot103

Q <- quantile(plot_norm_ratios$cv, probs=c(.25, .75), na.rm = TRUE)
iqr <- IQR(plot_norm_ratios$cv, na.rm = TRUE)

# Define bounds for what constitutes an outlier
lower_bound <- Q[1] - 1.5 * iqr
upper_bound <- Q[2] + 1.5 * iqr

# Filter out the outliers
filtered_data <- plot_norm_ratios %>%
  filter(cv >= lower_bound & cv <= upper_bound) %>%
  mutate(CV = 100 * cv)

# Plot without outliers
plot_filtered <- ggplot(filtered_data, 
                        aes(x = rescaled_intDIL, y = CV)) +
  geom_point(aes(colour = sample), alpha = 0.4, size = 0.3) +
  geom_smooth(method = "loess", 
              se = F, span = 1, colour = "black") +
  theme_minimal() +
  labs(x = "L Intensity (log10)",
       y = "CV",
       colour = "Dilution")
plot_filtered

final_plot <- plot101 + plot_filtered+ plot103  +plot102 + 
  plot_layout(guides = "collect",
              ncol = 2, heights = c(1,1)) &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing = margin())

string_plot <- paste0(path_plot_output, "boxplot_quantifies_thesis_PSC_",PSC,".pdf")

ggsave(string_plot, 
       plot = final_plot, width = 5, height = 5, units = "in", dpi = 300)
print(final_plot)

final_plot <- plot101 + plot103 + plot102 +
  plot_layout(guides = "collect",
              ncol = 3, heights = c(1),
              widths = c(3,2,2)) &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing = margin())

string_plot <- paste0(path_plot_output, "boxplot_quantifies_main_PSC_",
                      PSC, "LidMin_", L_idd_min ,".pdf")
ggsave(string_plot, 
       plot = final_plot, width = 8, height = 3, units = "in", dpi = 300)
print(final_plot)

# Number of Heavy IDs along dilution series ------------------------------------
reporttsv_JPTlib %>%
  filter(jptseq %in% jptseqs_inlib) %>% 
  filter(Channel.Q.Value < 0.05,
         PTM.Q.Value < 0.05,
         PTM.Site.Confidence > PSC,
         SILAClabel == "-H") %>% 
  group_by(Run, sample) %>% 
  summarise("Phospho-sequences" = n_distinct(jptseq),
            Precursors = n_distinct(Precursor.Id)) %>%
  mutate(sample = gsub("HEK", "", sample)) -> summ_Hids

# Reshape data to long format for plotting
data_long <- reshape2::melt(summ_Hids, id.vars = c('Run', 'sample'), variable.name = 'group', value.name = 'value')

# Calculate summary statistics: mean and standard deviation
summary_stats <- data_long %>%
  group_by(sample, group) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = 'drop')

data_long$sample <- factor(data_long$sample, 
                           levels = c('400ng', '200ng', '100ng', '50ng', '25ng', '0ng'))
summary_stats$sample <- factor(summary_stats$sample, 
                               levels = c('400ng', '200ng', '100ng', '50ng', '25ng', '0ng'))
# Create the plot

H_ids_plot <- ggplot(data_long, aes(x = sample, color = group)) +
  geom_errorbar(data = summary_stats, 
                aes(x = sample, ymin = mean - sd, ymax = mean + sd), width = 0.2,
                color = "black") +  # Error bars for SD
  geom_jitter(aes(y = value), 
              width = 0.2, size = 2, alpha = 0.6) +  # Jitter plot for individual points
  labs(y = '# identified \n(heavy)') +
  scale_y_continuous(limits = c(0, NA)) +  # Set y-axis to start from 0
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank())

string_plot <- paste0(path_plot_output, "Hids_PSC_",PSC, "LidMin_", L_idd_min ,"_forThesis.pdf")
ggsave(string_plot, 
       plot = H_ids_plot, width = 2.3, height = 2.3, units = "in", dpi = 300)
print(final_plot)
