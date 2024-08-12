## Code to prepare `methylation_in_tissues` dataset goes here

library(GenomicRanges)
library(SummarizedExperiment)
library(Biostrings)
library(tidyverse)
library(liftOver)
library(AnnotationHub)
load("../../../CTdata_extdata/all_genes_prelim.rda")

# Download human genome fasta file
# download.file("http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
#               destfile = "../extdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
hg <- rtracklayer::import(
  "../../../CTdata_extdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  format = "fasta")

################################################################################
# Create a GRange with all CpG positions in hg38
################################################################################

## Initiate CT_methylation tibble
## Use whole genome => identify all CpGs positions
CpG_positions <- tibble()
chromosomes <- c(1:22, "X", "Y")
for (chr_number in chromosomes) {
  chromosome_ind <- grep(hg@ranges@NAMES,
                         pattern = paste0('dna:chromosome chromosome:GRCh38:',
                                          chr_number, ':'))
  # finds the row corresponding to selected chromosome
  chr_seq <- hg[chromosome_ind][[1]] # gets the sequence of selected chromosome
  chr_CG <- matchPattern("CG", chr_seq)
  # gets the coords of CpG dinucleotides on selected chromosome
  tmp <- as_tibble(as.data.frame(chr_CG@ranges)) %>%
    mutate(chr = paste0("chr", chr_number)) %>%
    mutate(CpG_pos = start) %>%
    dplyr::rename(Cytosine_fw = start, Cytosine_rev = end) %>%
    dplyr::select(chr, CpG_pos, Cytosine_fw, Cytosine_rev)
  CpG_positions <- rbind(CpG_positions, tmp)
}

CpG_positions_gr <- makeGRangesFromDataFrame(CpG_positions,
                                             keep.extra.columns = TRUE,
                                             seqnames.field = "chr",
                                             start.field = "Cytosine_fw",
                                             end.field = "Cytosine_rev")





################################################################################
## Create a table with methylation values of all CpG in normal tissues
################################################################################
met_in_tissues <- as_tibble(CpG_positions_gr) %>%
  dplyr::rename(chr = seqnames) %>%
  dplyr::select(chr, CpG_pos)
Encode_tissues <- c("adipose", "colon", "esophagus", "heart", "intestine",
                    "lung", "muscle", "pancreas", "placenta", "skin", "stomach",
                    "testis", "thyroid")

## WGBS bed files downloaded from ENCODE website
## adipose ENCFF318AMC.bed
## colon ENCFF157POM.bed
## esophagus ENCFF625GVK.bed
## heart ENCFF536RSX.bed
## intestine ENCFF241AQC.bed
## lung ENCFF039JFT.bed
## muscle ENCFF121ZES.bed
## pancreas ENCFF763RUE.bed
## placenta ENCFF437OKM.bed
## skin ENCFF219GCQ.bed
## stomach ENCFF497YOO.bed
## testis ENCFF715DMX.bed
## thyroid ENCFF223LJW.bed

for (cell in Encode_tissues) {

  bismark_file <- paste0("../../../CTdata_extdata/", cell, ".bed")
  bismark <- read_tsv(bismark_file, col_types = cols(.default = "?", X1 = 'c'),
                      col_names = FALSE)

  bismark <- bismark %>%
    dplyr::select(X1, X2, X6, X10, X11) %>%
    dplyr::rename(chr = X1, Start = X2, strand = X6,
                  reads = X10, score = X11) %>%
    filter(reads > 1) %>%
    mutate(start = Start + 1) %>%  # change coords to 1-based system
    mutate(CpG_pos = case_when(strand == '+' ~ start,
                               strand == '-' ~ start - 1)) %>%
    dplyr::select(chr, CpG_pos, strand, score)

  # Merge methylation value of fw-score and rev-score
  bismark_merge <- as_tibble(CpG_positions_gr) %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::select(chr, CpG_pos) %>%
    left_join(bismark %>%
                filter(strand == '+'),
              by = c("chr", "CpG_pos")) %>%
    dplyr::rename(fw_score = score) %>%
    dplyr::select(- strand) %>%
    left_join(bismark %>%
                filter(strand == '-'),
              by = c("chr", "CpG_pos")) %>%
    dplyr::rename(rev_score = score) %>%
    dplyr::select(- strand) %>%
    group_by(chr, CpG_pos) %>%
    dplyr::summarise(score = mean(c(fw_score, rev_score), na.rm = TRUE))

  names(bismark_merge) <- c("chr", "CpG_pos", cell)

  # Add methylation values for the tissue in CpG_methylation tibble
  met_in_tissues <- met_in_tissues %>%
    left_join(bismark_merge, by = c("chr", "CpG_pos"))
  rm(bismark_merge)
}

## Use sperm WGBS from SRR15427118
## raw data was processed with bismark (version 0.20.0)
bismark_file <- "../../../CTdata_extdata/SRR15427118_1_val_1_bismark_bt2_pe.bismark.cov.gz"
bismark <- read_tsv(bismark_file, col_types = cols(.default = "?", X1 = 'c'),
                    col_names = FALSE)

bismark <- bismark %>%
  dplyr::rename(chr = X1, start = X2, end = X3, score = X4,
                met = X5, unmet = X6) %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  mutate(nReads = met + unmet) %>%
  filter(nReads > 1) %>%
  dplyr::select(-nReads)

CG_pos_strand <- as_tibble(CpG_positions_gr) %>%
  dplyr::rename(chr = seqnames) %>%
  left_join(bismark, by = c("chr" = "chr", "start" = "start")) %>%
  dplyr::mutate(end = start + 1) %>%
  dplyr::select(chr, start, end, width, strand, CpG_pos, score, met, unmet)
CG_neg_strand <- as_tibble(CpG_positions_gr) %>%
  dplyr::rename(chr = seqnames) %>%
  left_join(bismark, by = c("chr" = "chr", "end" = "start")) %>%
  dplyr::select(chr, start, end, width, strand, CpG_pos, score, met, unmet)

bismark_merge <- rbind(CG_pos_strand, CG_neg_strand) %>%
  group_by(chr, CpG_pos) %>%
  dplyr::summarise(sperm = mean(c(score), na.rm = TRUE))

met_in_tissues <- met_in_tissues %>%
  left_join(bismark_merge)


#  Save as a RangedSummarizeExperiment
met_normal_tissues_gr <-
  makeGRangesFromDataFrame(met_in_tissues,
                           keep.extra.columns = FALSE,
                           seqnames.field = "chr",
                           start.field = "CpG_pos",
                           end.field = "CpG_pos")

CpG_methylation_in_tissues <- SummarizedExperiment(
  assays = as.matrix(met_in_tissues[, -(1:2)]),
  rowRanges = met_normal_tissues_gr)


## Keep only regions around each gene TSS
promoter_regions <- all_genes_prelim %>%
  mutate(start = transcription_start_site - 5000) %>%
  mutate(stop = transcription_start_site + 5000) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
  mutate(strand = case_when(strand == 1 ~ "+", strand == -1 ~ "-"))

prom_gr <- makeGRangesFromDataFrame(promoter_regions,
                                    keep.extra.columns = TRUE,
                                    seqnames.field = "chromosome_name",
                                    start.field = "start",
                                    end.field = "stop")


methylation_in_tissues <- subsetByOverlaps(CpG_methylation_in_tissues, prom_gr)

save(methylation_in_tissues, file = "../../eh_data/methylation_in_tissues.rda",
     compress = "xz",
     compression_level = 9)
