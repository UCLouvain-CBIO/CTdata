## Code to prepare `methylation_in_embryo` dataset goes here

library("tidyverse")
library("GenomicRanges")
library("data.table")
library("SummarizedExperiment")

load("../extdata/CpG_positions_all_prom_hg19.rda")


## Data from Single Cell DNA Methylome Sequencing of Human Preimplantation
## Embryos  (Zhu et al. 2018)
## 500 scWGBS files were downloaded from GEO (accession: GSE81233).
## These .bed files have been subsetted to only contain CpG methylation.


################################################################################
# Manage sample metadata to create a coldata
################################################################################

## We are using two files : information from their supplemental file and
## metadata from the downloading (with file names).

metadata <- read_tsv("../../../CTdata_extdata/Zhu_metadata")

##  Add file names to metadata
coldata <-
  tibble(sample =
           list.files("../../../CTdata_extdata/GSE81233_subsetted_data/",
                      pattern = "*.gz")) %>%
  mutate(Accession = gsub("_.*$", x = sample, replace = '')) %>%
  left_join(metadata) %>%
  mutate(Sample_Name = gsub("GSM\\w*_", x = sample, replace = '')) %>%
  mutate(Sample_Name = gsub(".CpGmet.bed.gz", x = Sample_Name, replace = '')) %>%
  dplyr::select(sample, Sample_Name, everything())

# In metadata GSM2773562-GSM2773565 are referred as "Mor09-Q" (source name)
# While the name of the files contain "Mor09-W" rather than "Mor09-Q"
# Confusing => Removing samples
coldata <- coldata[!coldata$sample %in%
                     grep('Mor09-W', x = coldata$sample, value = TRUE), ]

# Vilus files were not in metadata, adding the info manually
coldata$bulk_or_single_cell[coldata$Sample_Name %in%
                              c("9W-PBAT-villus-1", "9W-PBAT-villus-2")] <- "bulk cells"
coldata$cell_type[coldata$Sample_Name %in%
                    c("9W-PBAT-villus-1", "9W-PBAT-villus-2")] <- "Villus"
coldata$source_name[coldata$Sample_Name %in%
                      c("9W-PBAT-villus-1", "9W-PBAT-villus-2")] <- "Villus"


##  Add metadata from Supplemental Table and modify to join to coldata
suppl <- readxl::read_xlsx(
  "../../../CTdata_extdata/41588_2017_7_MOESM3_ESM.xlsx", sheet = 1)

suppl$Sample_Name[suppl$Sample_Name == "Villus_9W_rep1"] <- "9W-PBAT-villus-1"
suppl$Sample_Name[suppl$Sample_Name == "Villus_9W_rep2"] <- "9W-PBAT-villus-2"
suppl$Sample_Name[suppl$Sample_Name == "Villus_11W_rep1"] <- "11W-PBAT-villus-1"
suppl$Sample_Name[suppl$Sample_Name == "Villus_11W_rep2"]<- "11W-PBAT-villus-2"
suppl$Sample_Name[suppl$Sample_Name == "Heart_7W_rep1"] <- "7W-PBAT-heart-1"
suppl$Sample_Name[suppl$Sample_Name == "Heart_7W_rep2"] <- "7W-PBAT-heart-2"
suppl$Sample_Name[suppl$Sample_Name == "Heart_11W_rep1"] <- "11W-PBAT-heart-1"
suppl$Sample_Name[suppl$Sample_Name == "Heart_11W_rep2"] <- "11W-PBAT-heart-2"

coldata <- coldata %>%
  left_join(suppl)

# Problem with a sample and the others are duplicated : remove
coldata <- coldata %>%
  filter(!sample %in% c("GSM2773698_Te05-W-b99.CpGmet.bed.gz",
                        "GSM2986409_scBS-O-bst-4.CpGmet.bed.gz",
                        "GSM2986393_Mor08-Q-s167.CpGmet.bed.gz",
                        "GSM2986394_Mor08-Q-s168.CpGmet.bed.gz"))



## Complete coldata when NA, standardize or adding info

coldata$Stage[coldata$Sample_Name %in% c("scBS-Sp12", "scBS-Sp14")] <- "Sperm"
coldata$Stage[coldata$Sample_Name %in% c("scBS-Morula-3bulk",
                                         "scBS-Morula-5bulk")] <- "Morula"

coldata$`Single Cell or Bulk cell` <- coldata$bulk_or_single_cell
coldata$stage_level <- coldata$Stage

coldata <- coldata %>%
  dplyr::rename(low_quality = "Low Quality (True or False)")

# Add timing to stage
coldata$stage_level[coldata$`Developmental Time` %in%
                      c("early zygotic stage (10-11h after ICSI)")] <- "PN-early"
coldata$stage_level[coldata$`Developmental Time` %in%
                      c("mid- zygotic stage (22-23h after ICSI)")] <- "PN-mid"
coldata$stage_level[coldata$`Developmental Time` %in%
                      c("late zygotic stage (25h or even later after ICSI)")] <- "PN-late"
coldata$stage_level[coldata$cell_type %in% c("Heart")] <- "Heart"
coldata$stage_level[coldata$cell_type %in% c("Villus")] <- "Villus"

# Simplify
coldata$Stage[coldata$Stage %in% c("Blastocyst (ICM)", "Blastocyst (TE)",
                                   "Blastocyst (Mixed)")] <- "Blastocyst"
coldata$Stage[coldata$Stage %in% "PN"] <- "Zygote"

# Add info when missing stage
coldata$Stage[is.na(coldata$Stage) & coldata$cell_type == "ICM"] <- "Blastocyst"

# Fill NA
coldata$`Developmental Time`[coldata$cell_type == "Sperm"] <- "Sperm"
coldata$`Developmental Time`[coldata$Stage == "Blastocyst"] <- "Blastocyst"
coldata$`Developmental Time`[coldata$Stage == "Morula"] <- "Morula"
coldata$`Genotype of the embryo`[coldata$source_name == "Icm10-Q"] <- "XX"
coldata$stage_level[coldata$source_name == "Icm10-Q"] <- "Blastocyst (ICM)"

# Order using factors
coldata$Stage <- factor(coldata$Stage, levels = c("GV Oocyte", "MII Oocyte",
                                                  "Sperm", "Zygote", "2-cell",
                                                  "4-cell", "8-cell", "Morula",
                                                  "Blastocyst",
                                                  "Post-implantation"))

coldata$stage_level <- factor(coldata$stage_level, levels =
                                 c("GV Oocyte", "MII Oocyte", "Sperm",
                                   "PN-early", "PN-mid" , "PN-late", "2-cell",
                                   "4-cell", "8-cell", "Morula",
                                   "Blastocyst (Mixed)", "Blastocyst (ICM)",
                                   "Blastocyst (TE)", "Heart", "Villus"))



# 'Sample_name' is not consistent
coldata$individual <- gsub("-\\w*$", x= coldata$Sample_Name, '')
coldata$individual[coldata$Sample_Name == "scBS-Morula-3bulk"] <- "scBS-Morula-3"
coldata$individual[coldata$Sample_Name == "scBS-Morula-5bulk"] <- "scBS-Morula-5"

# Fill NA
coldata$`Genotype of the embryo`[coldata$individual == "scBS-E-PN11"]<- "XY"
coldata$`Genotype of the embryo`[coldata$individual == "scBS-E-PN12"] <- "XX"
coldata$`Genotype of the embryo`[coldata$individual == "scBS-E-PN8"] <- "XX"
coldata$`Genotype of the embryo`[coldata$individual == "scBS-E-PN9"]<- "XY"
coldata$`Genotype of the embryo`[coldata$individual == "scBS-2C-2"] <- NA
coldata$`Genotype of the embryo`[coldata$individual == "scBS-hSP"] <- NA
coldata$`Genotype of the embryo`[coldata$individual == "scBS-Morula-3"] <- "XY"


################################################################################
# Extract methylation from CpG located in all promoters from each file
################################################################################

samples <- unique(coldata$Sample_Name)

# Create all position to keep all CpG even if NA
methylation_values <- tibble(chr =
                               as.vector(CpG_positions_all_prom_hg19@seqnames),
                             CpG_pos =
                               CpG_positions_all_prom_hg19@ranges@start)

setDTthreads(percent = 90) # to maximise cpu usage

for (sample_id in samples) {

  print(which(coldata$Sample_Name == sample_id))

  met_file <- paste0("../../../CTdata_extdata/GSE81233_subsetted_data/",
                     filter(coldata, Sample_Name == sample_id)$sample)

  met <- as_tibble(fread(met_file)) %>%
    dplyr::rename(chromosome_name = '#Chr', strand = Chain, start = Pos) %>%
    mutate(stop = start) %>%
    dplyr::select(chromosome_name, strand, start, stop, Met, UnMet)

  met_GR <- makeGRangesFromDataFrame(met,
                                     keep.extra.columns = TRUE,
                                     seqnames.field = "chromosome_name",
                                     start.field = "start",
                                     end.field = "stop")

  rm(met)

  # Keep CpG from promoters
  tmp <- mergeByOverlaps(CpG_positions_all_prom_hg19, met_GR)
  rm(met_GR)

  # Calculate mean met by CpG position (merge + and - values)
  tmp_met <- as_tibble(tmp) %>%
    dplyr::select(CpG_positions_all_prom_hg19.seqnames,
                  CpG_pos, Met, UnMet) %>%
    dplyr::rename(chr = CpG_positions_all_prom_hg19.seqnames) %>%
    group_by(chr, CpG_pos) %>%
    dplyr::summarize(met = sum(Met), unmet = sum(UnMet)) %>%
    mutate(methylation = round(met / (met + unmet) * 100, 1)) %>%
    dplyr::select(chr, CpG_pos, methylation) %>%
    ungroup()
  rm(tmp)

  names(tmp_met) <- c("chr", "CpG_pos", sample_id)

  methylation_values <- full_join(methylation_values, tmp_met)
  rm(tmp_met)
  gc()
}

rm(sample_id)


################################################################################
# Save as a RangedSE
################################################################################

methylation_values_GR <- makeGRangesFromDataFrame(methylation_values %>%
                                                    dplyr::select(chr,
                                                                  CpG_pos) %>%
                                                    dplyr::rename(start =
                                                                    CpG_pos) %>%
                                                    mutate(stop = start + 1),
                                                  keep.extra.columns = FALSE,
                                                  seqnames.field = "chr",
                                                  start.field = "start",
                                                  end.field = "stop")

coldata <- dplyr::select(as.data.frame(coldata), -command)
rownames(coldata) <- coldata$Sample_Name


methylation_in_embryo <- SummarizedExperiment(assays = methylation_values[, coldata$Sample_Name],
                                              rowRanges = methylation_values_GR,
                                              colData = coldata)

save(methylation_in_embryo, file = "../../eh_data/methylation_in_embryo.rda",
     compress = "xz",
     compression_level = 9)

