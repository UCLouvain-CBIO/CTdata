## Code to prepare `methylation_in_FGC` dataset goes here

library("tidyverse")
library("GenomicRanges")
library("data.table")
library("SummarizedExperiment")

load("../extdata/CpG_positions_all_prom_hg19.rda")


## Data from Dissecting the epigenomic dynamics of human fetal germ cell
## development at single-cell resolution (Li et al. 2021)
## 461 scWGBS files were downloaded from GEO (accession: GSE107714).


################################################################################
# Manage sample metadata
################################################################################

files <- list.files("../../../CTdata_extdata/GSE107714/", pattern = "txt.gz")
coldata <- tibble(files = files,
                  sample = gsub("GSM\\d*_", '', files)) %>%
  mutate(type = gsub('_\\d*W.*', '', sample)) %>%
  mutate(type = gsub('\\w_', '', type)) %>%
  mutate(sample = gsub(".methylation_calling.txt.gz", '', sample)) %>%
  mutate(sex = gsub('_FGC.*', '', sample)) %>%
  mutate(sex = gsub('_Soma.*', '', sex)) %>%
  mutate(time_week = gsub('W.*', '', sample)) %>%
  mutate(time_week = gsub('.*_', '', time_week)) %>%
  mutate(time_week = as.numeric(time_week)) %>%
  mutate(individual = gsub('.*W_', '', sample)) %>%
  mutate(individual = gsub('_.*', '', individual)) %>%
  mutate(cell = gsub('.*sc', 'sc', sample)) %>%
  mutate(group = paste0(sex, '_', type, '_', time_week, 'W')) %>%
  mutate(group = factor(group,
                        levels = c("F_FGC_8W", "F_FGC_9W", "F_FGC_10W",
                                   "F_FGC_17W", "F_FGC_21W", "M_FGC_6W",
                                   "M_FGC_7W", "M_FGC_17W", "M_FGC_21W",
                                   "M_FGC_24W", "F_Soma_8W", "F_Soma_9W",
                                   "F_Soma_17W", "F_Soma_21W", "M_Soma_6W",
                                   "M_Soma_7W", "M_Soma_17W", "M_Soma_21W"))) %>%
  mutate(sex_time = paste0(sex, '_', time_week, "_", type)) %>%
  mutate(sex_time = factor(sex_time,
                           levels = c("F_8_FGC", "F_9_FGC", "F_10_FGC",
                                      "F_17_FGC", "F_21_FGC", "M_6_FGC",
                                      "M_7_FGC", "M_17_FGC", "M_21_FGC",
                                      "M_24_FGC", "M_7_Soma", "F_8_Soma"))) %>%
  mutate(sex_type = paste0(sex, "_", type)) %>%
  mutate(sex_type = factor(sex_type,
                           levels = c("F_FGC", "M_FGC", "F_Soma", "M_Soma")))


metadata <- readxl::read_xlsx("../../../CTdata_extdata/GSE107714/41422_2020_401_MOESM8_ESM.xlsx")
names(metadata)[1] <- "sample"

coldata <- coldata %>%
  left_join(metadata) %>%
  dplyr::select(- Gender)

coldata <- coldata %>%
  filter(`Pass Quality Control (True or False)`) %>%
  mutate(cellType = case_when(`Cell Type` == "FGC2" &
                                sex == "M" ~ "M_mitotic_FGC",
                              `Cell Type` == "FGC3" &
                                sex == "M"~ "M_mitotic_arrested_FGC",
                              `Cell Type` == "FGC1" &
                                sex == "M" ~ "M_FGC",
                              `Cell Type` == "FGC1" &
                                sex == "F" ~ "F_FGC",
                              `Cell Type` == "FGC2" &
                                sex == "F" ~ "F_mitotic_FGC",
                              `Cell Type` == "FGC3" &
                                sex == "F" ~ "F_meiotic_FGC",
                              `Cell Type` == "FGC4" &
                                sex == "F"~ "F_oogenesis_FGC",
                              `Cell Type` == "Soma" &
                                sex == "F" ~ "F_soma",
                              `Cell Type` == "Soma" &
                                sex == "M" ~ "M_soma")) %>%
  mutate(cellType = factor(cellType,
                           levels = c("F_FGC", "F_meiotic_FGC",
                                      "F_oogenesis_FGC", "M_mitotic_FGC",
                                      "M_mitotic_arrested_FGC", "F_soma",
                                      "M_soma")))



################################################################################
# Extract methylation from CpG located in all promoters from each file
################################################################################

# We don't need to keep all soma samples, juste 1 male and 1 female is enough

keep <-
  coldata[c(grep("FGC", coldata$group),
            which(coldata$group %in% c("F_Soma_8W", "M_Soma_7W"))), ]

samples <- unique(keep$sample)

# Create all position to keep all CpG even if NA
methylation_values <- tibble(chr =
                               as.vector(CpG_positions_all_prom_hg19@seqnames),
                             CpG_pos =
                               CpG_positions_all_prom_hg19@ranges@start)

setDTthreads(percent = 90) # to maximise cpu usage

for (sample_id in samples) {

  print(which(keep$sample == sample_id))

  met_file <- paste0("../../../CTdata_extdata/GSE107714/",
                     filter(coldata, sample == sample_id)$files)

  met <- as_tibble(fread(met_file)) %>%
    filter(Type == "CpG") %>%
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

coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample

coldata <- coldata[coldata$sample %in% names(methylation_values), ]

methylation_in_FGC <- SummarizedExperiment(assays = methylation_values[, coldata$sample],
                                           rowRanges = methylation_values_GR,
                                           colData = coldata)

save(methylation_in_FGC, file = "../../eh_data/methylation_in_FGC.rda",
     compress = "xz",
     compression_level = 9)

