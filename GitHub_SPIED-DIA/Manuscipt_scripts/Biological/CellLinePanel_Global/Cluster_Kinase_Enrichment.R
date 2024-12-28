#SS2 PTM-SEA and other enrichment analyses, after diff expression
library(tidyverse)
library(cmapR)
library(ggrepel)
library(patchwork)
library(scales)

# function copied from broadinstitute/ssGSEA.2.0 (to minimize download of dependencies)
Read.GeneSets.db2 <- function (gs.db, thres.min = 2, thres.max = 2000) {
  ## read gmt files
  temp <- readLines(gs.db)
  temp <- strsplit(temp, '\t')
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  
  ## filter gene sets according to size
  rm.idx <- which(temp.size.G < thres.min | temp.size.G > thres.max)
  if(length(rm.idx) > 0){
    temp <- temp[-rm.idx]
    temp.size.G <- temp.size.G[-rm.idx]
  }
  
  max.Ng <- length(temp)         ## number of signature sets
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  max.size.G <- max(temp.size.G) ## maximal size
  
  gs <- lapply(temp, function(x)x[3:length(x)])
  gs.names <- sapply(temp, function(x)x[1])
  gs.desc <- sapply(temp, function(x)x[2])
  
  ## check whether gene sets are unique
  gs.unique <- lapply(gs, unique)
  gs.unique.size.G <- sapply(gs.unique, length)
  gs.not.unique.idx <- which(gs.unique.size.G < temp.size.G)
  if( length(gs.not.unique.idx) > 0 ){
    warning("\n\nDuplicated gene set members detected. Removing redundant members from:\n\n", paste(gs.names[gs.not.unique.idx], collapse='\n'))
    gs <- gs.unique
    temp.size.G <- gs.unique.size.G 
  }
  size.G <- temp.size.G
  names(gs) <- names(gs.names) <- names(gs.desc) <- names(size.G) <- gs.names
  
  
  return(list(N.gs = max.Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc,
              size.G = size.G, max.N.gs = max.Ng))
}

  
# set path to relevant database (downloaded from ptm-sea and ikip-db)
GSDB <- Read.GeneSets.db2("/ssGSEA2.0-master/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt",
                                    thres.min = 10, thres.max = 2000)
# Assuming gs is your original list of kinases with sequences
kinases_list <- GSDB$gs

####Flatten the kinase list into a data frame
kinase_sequence_pairs <- do.call(rbind, lapply(names(kinases_list), function(kinase) {
  data.frame(kinase = kinase, sequence = kinases_list[[kinase]], stringsAsFactors = FALSE)
})) %>%
  mutate(sequence = gsub("-p;[du]$", "", sequence),
         kinase = sub("KINASE-iKiP_(.*)", "\\1(i)", kinase),
         kinase = sub("KINASE-PSP_(.*)", "\\1(P)", kinase),
         kinase = gsub("RPS6KA[45]\\/|MAPK[0-9]*\\/|Akt1\\/|mTOR\\/|CSNK2A2\\/|CSNK2A1\\/",
                       "", kinase),
         kinase = gsub("\\/MAPK[0-9]*|Akt1\\/|mTOR\\/|\\/PRKACA|\\/CHEK1|\\/CSNK2A2|\\/CSNK2A1|\\/PRKAA1",
                       "", kinase),
         kinase = gsub("PRKACA", "PKACA", kinase),
         kinase = gsub("RPS6KA1", "p90RSK",kinase),
         kinase = gsub("p90RSK\\/p90RSK", "p90RSK", kinase),
         kinase = gsub("MAPK11", "P38B", kinase),
         kinase = gsub("MAPK12", "P38G", kinase),
         kinase = gsub("MAPK13", "P38D", kinase),
         kinase = gsub("MAPK14", "P38A", kinase),
         kinase = gsub("RPS6KA4.MSK2", "MSK2", kinase),
         kinase = gsub("MAPK9.JNK2", "JNK2", kinase),
         kinase = gsub("MAPK8.JNK1", "JNK1", kinase),
         kinase = gsub("CSNK2A1.CK2A1", "CK2A1", kinase),
         kinase = gsub("PKCB\\/PRKCB", "PKCB", kinase)) %>%

  filter(!grepl("PERT", kinase),
         !grepl("DISEASE", kinase),
         !grepl("PATH", kinase))
saveRDS(kinase_sequence_pairs, "Data/Cluster_Kinase_sequence_pairs.rds")
kinase_sequence_pairs <- readRDS("Data/Cluster_Kinase_sequence_pairs.rds")

#file with flanking sequences per phospho-site in library (not included! )
flankingseqsanno <- read_csv("HpH_flankingseq_20240116.csv")

# Function to find kinases for given sequences
find_kinases_for_sequences <- function(sequence_vector, kinase_sequence_pairs) {
  sapply(sequence_vector, function(seq) {
    matched_kinases <- kinase_sequence_pairs$kinase[kinase_sequence_pairs$sequence %in% seq]
    paste(unique(matched_kinases), collapse=", ")
  })
}

fun_PrecurIDto_sty <- function(seq_p){
  kinlibseq <- gsub("S\\(Phospho \\(STY\\)\\)","s", seq_p)
  kinlibseq <- gsub("Y\\(Phospho \\(STY\\)\\)","y", kinlibseq)
  kinlibseq <- gsub("T\\(Phospho \\(STY\\)\\)","t", kinlibseq)
  kinlibseq <- gsub("S\\(UniMod:21\\)","s", kinlibseq)
  kinlibseq <- gsub("Y\\(UniMod:21\\)","y", kinlibseq)
  kinlibseq <- gsub("T\\(UniMod:21\\)","t", kinlibseq)
  kinlibseq <- gsub("\\(UniMod:\\d+\\)", "", kinlibseq)
  kinlibseq <- gsub("\\(ph\\)S", "s", kinlibseq)
  kinlibseq <- gsub("\\(ph\\)T", "t", kinlibseq)
  kinlibseq <- gsub("\\(ph\\)Y", "y", kinlibseq)
  kinlibseq <- gsub("\\(Acetyl \\(Protein N-term\\)\\)","", kinlibseq)
  kinlibseq <- gsub("\\(Acetyl \\(N-term\\)\\)","", kinlibseq)
  kinlibseq <- gsub("\\(SILAC-[R|K]-L\\)","", kinlibseq)
  kinlibseq <- gsub("\\(M\\(Oxidation \\(M\\)\\)", "", kinlibseq)
  kinlibseq <- gsub("\\(Oxidation \\(M\\)\\)", "", kinlibseq)
  kinlibseq <- gsub("\\)[2-4]*$|^[^\\(]*\\(|\\(ox\\)|\\(ac\\)", "", kinlibseq)
  kinlibseq <- gsub("[2-4]*$", "", kinlibseq)
  kinlibseq <- gsub("Acetyl \\(Protein N-term\\)\\)","", kinlibseq)
  kinlibseq <- gsub("Acetyl \\(N-term\\)\\)","", kinlibseq)
  kinlibseq <- gsub("\\(SILAC-[R|K]-L","", kinlibseq)
  kinlibseq <- gsub("SILAC-[R|K]-L\\)","", kinlibseq)
  kinlibseq <- gsub("M\\(Oxidation \\(M\\)\\)", "", kinlibseq)
  return(kinlibseq)
}
split_sequences <- function(sequence) {
  # Find the positions of the phosphorylated residues
  positions <- gregexpr("[sty]", sequence, ignore.case = FALSE)[[1]]
  
  # Create a vector to store the individual sequences
  individual_sequences <- vector()
  
  # Loop through each position and create the modified sequence
  for (pos in positions) {
    single_phosphorylated <- chartr("sty", "STY", sequence)
    single_phosphorylated <- substr_replace(single_phosphorylated, tolower(substr(single_phosphorylated, pos, pos)), pos, pos)
    individual_sequences <- c(individual_sequences, single_phosphorylated)
  }
  
  return(individual_sequences)
}
substr_replace <- function (string, replacement, start, stop) {
  substr(string, start, stop) <- replacement
  return(string)
}

cluster_tocheck <- 6


plot_clusterenrichement <- function(cluster_tocheck, CellLine, cluster_HCT116_flanking){
  cluster_HCT116_flanking %>%
    filter(cluster %in% cluster_tocheck) %>% 
    pull(flanking) -> flanking_incluster
  
  cluster_HCT116_flanking %>%
    filter(!(cluster %in% cluster_tocheck)) %>% 
    pull(flanking) -> flanking_background
  
  # find kinases for which targets were found in cluster of interest vs background
  interesting_kinases <- find_kinases_for_sequences(flanking_incluster, kinase_sequence_pairs)
  background_kinases <- find_kinases_for_sequences(flanking_background, kinase_sequence_pairs)
  
  # Convert the matched kinases from comma-separated strings to lists
  interesting_kinases_list <- str_split(interesting_kinases, ", ")
  background_kinases_list <- str_split(background_kinases, ", ")
  
  # number of targets found per kinase in background vs cluster of interest
  kinase_counts_interesting <- table(unlist(interesting_kinases_list))
  kinase_counts_background <- table(unlist(background_kinases_list))
  
  # Convert to data frames
  df_interesting <- as.data.frame(kinase_counts_interesting)
  df_background <- as.data.frame(kinase_counts_background)
  
  # Correct column names for joining
  names(df_interesting) <- c("kinase", "Seq_interesting")
  names(df_background) <- c("kinase", "Seq_background")
  
  # number sequences found in background vs target cluser. 
  total_interesting_sequences <- length(flanking_incluster)
  total_background_sequences <- length(flanking_background)
  
  # Join and calculate enrichment
  df_enrichment <- left_join(df_interesting, df_background, by = "kinase") %>%
    mutate(
      #ensure NA's are zero's  after joining
      N_interesting = ifelse(is.na(Seq_interesting), 0, Seq_interesting),
      N_background = ifelse(is.na(Seq_background), 0, Seq_background), 
      #frequency of kinase-targets (per target) in cluster of interest vs background
      Freq_interesting = N_interesting/total_interesting_sequences,
      Freq_background = N_background/total_background_sequences,
      #calculate enrichment ratio
      enrichment_ratio = log10(Freq_interesting / Freq_background)
    ) %>% 
    rowwise() %>%
    mutate(enrichment_significance = fisher.test(matrix(c(N_interesting, total_interesting_sequences - N_interesting,
                                                          N_background, total_background_sequences - N_background),
                                                        ncol = 2))$p.value)
  
  df_enrichment <- df_enrichment %>%
    mutate(neg_log_p_value = -log10(enrichment_significance)) %>%
    arrange(desc(neg_log_p_value)) %>%
    filter(str_length(kinase) > 1,
           N_background>0,
           N_interesting>2)
  
  path_id <- paste0(CellLine, "_", cluster_tocheck)
  
  saveRDS(df_enrichment, file.path(data_dir, paste0(path_id,".rds")))
  
  # Select top N kinases for visualization
  top_n_kinases <- 10
  df_top_enriched <- head(df_enrichment, top_n_kinases)
  
  width_help <- max(str_length(df_top_enriched$kinase))
  
  cluster_tocheck <- paste0(cluster_tocheck, collapse = "_")
  title <- paste0(CellLine, " cluster ", cluster_tocheck)
  
  
  ggplot(df_top_enriched, aes(x = enrichment_ratio, y = reorder(kinase, enrichment_ratio),
                              size = N_interesting, color = neg_log_p_value)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
    scale_size_continuous(range = c(1, 5), 
                          breaks = pretty_breaks(n = 4), # Adjust this to control the number of breaks in the legend
                          limits = c(0, NA)) + # This sets the min
    # guides(color = guide_legend(title.position = "top"),
    labs(subtitle = title,
         x = "Enrichment Ratio (log10)",
         size = "N precursor") +
    theme_minimal() +
    guides(size = guide_legend(title.position = "top")) +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.subtitle = element_text(size = 15, face = "bold"),
          legend.position = "right",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.title.align = 0.5,
          legend.spacing = unit(0.05, "cm"))
  
  plot_path<-file.path(output_dir, paste0("Enrichmentlog10", path_id,".pdf")) 
  
  ggsave(filename = plot_path,
         height = 3.5, 
         width = (width_help/16) + 5)
}

