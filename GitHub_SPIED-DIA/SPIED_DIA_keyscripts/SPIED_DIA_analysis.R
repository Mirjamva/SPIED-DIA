################################################################################
# Title: SPIED-DIA analysis
#
# Author: Mirjam van Bentum
# Date:   2024-12-01
#
# Description:
# This scripts performs the most important steps of SPIED-DIA analysis. 
# it is by no means ready to run without adjustments for paths, experimental design etc.
# 
# it requires a table with ordered heavy peptides and the DIA-NN output of a SPIED-DIA analysis as outlined in the GitHub readme
# 
# It includes:
#  - Loading and cleaning data from DIA-NN analysis
#  - SPIED-DIA identification and quantification of target peptides


################################################################################
# ---------------------------- #
#         Libraries            #
# ---------------------------- #

# Load necessary libraries
library(tidyverse)
library(data.table)

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
annotated_peptides_file <- file.path(data_dir, "ordered_peptides.csv")
# TSV file (DIA-NN "report.tsv") containing SILAC target data:
silac_report_file <- file.path(data_dir,  "report.tsv")


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
    # experiment specific settings, left in for example purposes. 
    Run_simple = gsub("Popeye_20230|_MVB_HStdia_SS2|_[A-Z0-9]*_[0-9]_[0-9]*$", "", Run)
  ) %>%
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
  # depending on application Channel.Q.Value might need to be reduced
  mutate(
    run_precursor_id = paste0(
      Run, 
      gsub("-L\\)", "-H\\)", Precursor.Id)
    )
  ) %>%
  pull(run_precursor_id)

# Normalize and filter data
jpt_2norm <- jpt_fil %>% 
  # identify heavy confident identfications
  filter(
    Channel.Q.Value < 0.05,
    PTM.Q.Value < 0.05,
    SILAClabel == "-H",
    Channel.H > 1000
  ) %>% 
  #makr light medium confident identifications
  mutate(
    run_precursor_id = paste0(Run, Precursor.Id),
    cell_line = str_extract(sample_simple, "C|D|H"),
    L_idd = run_precursor_id %in% L_ids
  ) %>% 
  # heavy precursor needs to be consistently identified 
  group_by(Precursor.Id, cell_line) %>%
  filter(n() > 9) %>%
  # calculate rescaled intensities. 
  mutate(rescalingfct = median(Ms1.Area)) %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(rescaled_int = log10(Channel.L / Channel.H * rescalingfct)) %>% 
  # Target id needs to be identified with medium confidence in the equivalent of at least one experimental group
  # example code, change as needed:
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

# continue with normalization, differential expression analysis etc. see scripts manuscript for example workflow