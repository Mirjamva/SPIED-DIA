# Title: Convert MaxQuant Output to DIA-NN Library
# Author: Mirjam van Bentum
# Date: 01.11.2024
# Description:
#   This script reads MaxQuant output files and converts them into a DIA-NN compatible library.
#   requires DDA analysis of heavy peptide standard only, processed in MaxQuant.
#   It processes peptide sequences, applies modifications, handles duplicate matches, and outputs a CSV file.
#   Based on MQ version 2.4.0.0 and DIA-NN version 1.8.2_beta11
#   R 4.3.0 
#   
#   There are two ways to set up the MQ search.
#   1. K/R SILAC labeled: standard K/R SILAC labeling, MQ will search for only for fully labeled peptides meaning all K & R residues found in a peptide are heavy labeled. 
#   2. manually added modifications: some heavy peptides contain missed cleavages or SIL labeled amino acids other then K/R. 
#      by manually adding the SIL labels as modificaitons you allow MaxQuant to search for peptides containing more then one K|R to search for different combinations. often you only want a heavy labeled C-terminus
#   
#   (instructions copied from Tjorven Hinzke tutorial, 2016 MaxQuant Summer School presentations (MaxQuant version 1.5.5.1))
#    in MaxQuant:
#     - Configuration  Modification  Add
#     - Name: give arbitrary (meaningful) name for the modification you want to add 
#     eg: Lys8_asmod, Arg10_asmod, Pro6_asmod, Phe6_asmod
#     - Description: shortly describe the modification
#     - Composition: specify compositional change by clicking “Change” and then using the
#     periodic table or the dropdown menu (e.g. H(-6)C(-6)O(-1) for dehydroalanine from tyrosine) – select elements after each other and set loss/yield (field next to composition: shows the monoisotopic mass)
#     - Position: can modification be anywhere or is it fixed to (not) occur at certain positions
#     - Type: is modification e.g. a standard modification or a label (important for handling of the modification in MaxQuant – labels are used for quantification)
#       (in this case, specify as standard modification!!!!)
#     - New terminus: specify whether the modification adds a new C- or N-terminus
#     - Specificities: specify amino acids which are modified, neutral losses (can e.g. happen for phosphorylation)
#     - modification will be saved in MaxQuant and always be available
#     !!! close and re-open MaxQuant for the changes to take effect
#     also see: https://cox-labs.github.io/coxdocs/andromeda_modifications_table.html


# Load required libraries
library(tidyverse)
library(data.table)

# -------------------- User-Defined Parameters --------------------

# Set the paths to your input files
ordered_peptides_file <- "library_annotated.csv" # file with peptides in heavy library
msms_file <- "combined/txt/msms.txt"
evidence_file <- "combined/txt/evidence.txt"
output_file <- "converted_library.csv"

# -------------------- Functions --------------------

# Function to convert MaxQuant modified sequences to pSTY format
# pSTY format: simple peptide sequence ignoring all modifications except phospho
fun_modseq_mq20topSTY <- function(seq_p) {
  pSTYseq <- seq_p
  
  # Replace phosphorylations
  pSTYseq <- gsub("S\\(Phospho[[:space:]]?\\(STY\\)\\)", "pS", pSTYseq)
  pSTYseq <- gsub("Y\\(Phospho[[:space:]]?\\(STY\\)\\)", "pY", pSTYseq)
  pSTYseq <- gsub("T\\(Phospho[[:space:]]?\\(STY\\)\\)", "pT", pSTYseq)
  
  # Remove acetylations and oxidations
  pSTYseq <- gsub("\\(Acetyl[[:space:]]?\\(Protein[[:space:]]?N-term\\)\\)", "", pSTYseq)
  pSTYseq <- gsub("\\(Oxidation[[:space:]]?\\(M\\)\\)", "", pSTYseq)
  
  # Remove labels and modifications
  pSTYseq <- gsub("^_|_$", "", pSTYseq)
  #pSTYseq <- gsub("\\(Lys8_?asmod\\)|\\(Arg10_?asmod\\)|\\(Pro6_?asmod\\)|\\(Phe6_?asmod\\)|\\(Phe10_?asmod\\)|\\(Lys8_?label\\)", "", pSTYseq)
  
  return(pSTYseq)
}

# -------------------- Data Processing --------------------

# Read annotated JPT peptides
ordered_peptides <- read_csv(ordered_peptides_file) %>% 
  rowwise()
ids_ordered <- ordered_peptides$pSTYseq #generate pSTYseq column if required!

# Read MaxQuant msms.txt file
MQ_msms <- read_tsv(msms_file, 
                                col_types = cols(.default = "c")) # Read all columns as character

# Read MaxQuant evidence.txt file
MQ_evidence <- read_tsv(evidence_file, 
                              col_types = cols(.default = "c")) %>%
  mutate(Score = as.numeric(Score))

# Extract ion mobility information from evidence.txt
ev_IM_info <- MQ_evidence %>%
  select(id, `1/K0`)

# Process msms data to create DIA-NN library
temp <- MQ_msms %>%
  # Apply modification function to 'Modified sequence' column
  mutate(pSTY_modseq = fun_modseq_mq20topSTY(`Modified sequence`)) %>%
  
  # Keep only sequences in 'ids_ordered'
  filter(pSTY_modseq %in% ids_ordered) %>%
  # filter as necessecary in case of manually set modifications.
  
  # For each modified sequence, keep the one with the highest Score
  # !!! critical point: keep in all peptide-isoforms
  # peptides carrying Met_ox or n-terminal acetylation might have the highest score in DDA but perform bad quantitatively.
  group_by(`Modified sequence`) %>% 
  filter(Score == max(Score)) %>% 
  ungroup() %>%
  
  # Calculate precursor mass shift
  mutate(
    precursor_mzshift = str_extract(`Modified sequence`, pattern = "R_|K_"),
    precursor_mzshift = gsub("K_", "8.01419", precursor_mzshift),
    precursor_mzshift = as.numeric(gsub("R_", "10.0082", precursor_mzshift)), 
    `m/z` = as.numeric(`m/z`),
    mz_light = `m/z` - (precursor_mzshift / Charge)
  ) %>% 
  
  # Separate 'Matches', 'Intensities', and 'Masses' columns into rows
  separate_rows(Matches, Intensities, Masses, sep = ";") %>% 
  rowwise() %>%
  
  # Create new columns for DIA-NN library
  mutate(
    #ModifiedPeptide = gsub("\\((Lys8|Arg10)_asmod\\)", "", `Modified sequence`), #if necessary
    PrecursorCharge = as.numeric(Charge),
    PrecursorMz = as.numeric(mz_light),
    Tr_recalibrated = as.numeric(`Retention time`),
    ProductMz = as.numeric(Masses),
    LibraryIntensity = as.numeric(Intensities),
    UniprotID = Sequence,
    ProteinName = Sequence,
    Genes = Sequence,
    Proteotypic = 0,
    FragmentCharge = ifelse(
      str_detect(Matches, "\\(\\d+\\+\\)"),
      as.numeric(str_match(Matches, "\\((\\d+)\\+\\)")[,2]),
      1
    ),
    FragmentType = str_extract(Matches, "a|b|c|x|y|z"),
    FragmentType = gsub("a", "b", FragmentType),
    FragmentType = case_when(Matches == "pY" ~ "pY", TRUE ~ FragmentType),
    FragmentSeriesNumber = str_extract(gsub("-H2O|-NH3", "", Matches), "[0-9]+")
  )

# -------------------- Handle Duplicate Matches --------------------

# Identify duplicate fragment matches and resolve them
newfragnames <- temp %>%
  group_by(id, Matches) %>%
  filter(n() > 1) %>% 
  mutate(Difference = as.numeric(Masses) - lag(as.numeric(Masses))) %>% 
  ungroup() %>% 
  rowwise() %>%
  filter(Difference == min(Difference, na.rm = TRUE)) %>%
  mutate(
    identifier = paste0(id, Matches, Masses, Intensities),
    Matches = paste0(Matches, "_ploss")
  ) %>%
  select(-Difference)

# Remove duplicates and combine with updated fragment names
msms2diann <- temp %>%
  mutate(identifier = paste0(id, Matches, Masses, Intensities)) %>% 
  filter(!(identifier %in% newfragnames$identifier)) %>% 
  bind_rows(newfragnames) %>%
  ungroup() %>%
  select(-identifier) %>%
  
  # Add fragment loss type
  mutate(
    FragmentLossType = str_extract(Matches, "H2O_ploss|NH3_ploss|H2O|NH3|\\*"),
    FragmentLossType = gsub("\\*", "ploss", FragmentLossType),
    
    # Adjust ProductMz for light ions
    ProductMz_light = case_when(
      FragmentType == "y" ~ ProductMz - precursor_mzshift,
      FragmentType == "b" ~ ProductMz,
      FragmentType == "pY" ~ ProductMz
    ),
    ProductMz = ProductMz_light
  ) %>% 
  
  # Add ion mobility information
  left_join(ev_IM_info, by = c("Evidence ID" = "id")) %>%
  rename(IonMobility = `1/K0`) %>%
  
  # Select columns for DIA-NN library
  select(
    ModifiedPeptide, PrecursorCharge, 
    PrecursorMz, Tr_recalibrated, ProductMz,
    LibraryIntensity, UniprotID, ProteinName, 
    Genes, Proteotypic, FragmentCharge, 
    FragmentType, FragmentSeriesNumber, FragmentLossType, 
    IonMobility
  )

# filter out pY ions
msms2diann %>%
  filter(!is.na(FragmentType)) ->  msms2diann

# -------------------- Output --------------------

# Write the DIA-NN library to a CSV file
write_csv(msms2diann, output_file)