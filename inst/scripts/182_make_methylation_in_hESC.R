## Code to prepare `methylation_in_hESC` dataset goes here

library(GenomicRanges)
library(SummarizedExperiment)
library(Biostrings)
library(tidyverse)
library(liftOver)
library(AnnotationHub)
load("../../eh_data/all_genes.rda")
load("../../../CTdata_extdata/hESC_RNAseq_coldata.rda")


# Download human genome fasta file
# download.file("http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
#               destfile = "/data/cbio/GENOMES/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
hg <- rtracklayer::import(
  "../../../CTdata_extdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  format = "fasta")

################################################################################
# Create a GRange with all CpG positions located in all genes promoters
################################################################################

## Initiate all_methylation tibble
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


## Keep only regions around each TSS
all_promoter_regions <- all_genes %>%
  mutate(start = transcription_start_site - 5000) %>%
  mutate(stop = transcription_start_site + 5000) %>%
  mutate(chromosome_name = paste0("chr", chr)) %>%
  mutate(strand = case_when(strand == 1 ~ "+", strand == -1 ~ "-"))

all_list_prom_gr <- makeGRangesFromDataFrame(all_promoter_regions,
                                             keep.extra.columns = TRUE,
                                             seqnames.field = "chromosome_name",
                                             start.field = "start",
                                             end.field = "stop")


prom_CpG_positions_gr <- subsetByOverlaps(CpG_positions_gr, all_list_prom_gr)


################################################################################
## Create a table with methylation values of CpG located in promoters
## in hESC
################################################################################
met_in_hESC <- as_tibble(prom_CpG_positions_gr) %>%
  dplyr::rename(chr = seqnames) %>%
  dplyr::select(chr, CpG_pos)
Encode_cells <- c("H1", "HUES64", "H9")

## WGBS data from embryonic stem cells were downloaded from Encode.
## H1  (ENCFF601NBW.bed.gz)
## H9  (ENCFF975CVU.bed.gz)
## HUES64 (ENCFF770UYJ.bed.gz)


for (cell in Encode_cells) {

  bismark_file <- paste0("../../../CTdata_extdata/", cell, ".bed.gz")
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
  bismark_merge <- as_tibble(prom_CpG_positions_gr) %>%
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
  rm(bismark)

  # Add methylation values for the tissue in CpG_methylation tibble
  met_in_hESC <- met_in_hESC %>%
    left_join(bismark_merge, by = c("chr", "CpG_pos"))
  rm(bismark_merge)
}


################################################################################
## Add coldata and save as a Ranged SE
################################################################################

coldata <- as.data.frame(coldata) %>%
  mutate(cell_type = "hESC") %>%
  mutate(genotype = factor(genotype, levels = c("XX", "XY"))) %>%
  mutate(dataset = "Encode") %>%
  mutate(data = "Methylation") %>%
  dplyr::select(-files) %>%
  column_to_rownames("sample")


methylation_hESC_gr <- makeGRangesFromDataFrame(met_in_hESC,
                                                   keep.extra.columns = FALSE,
                                                   seqnames.field = "chr",
                                                   start.field = "CpG_pos",
                                                   end.field = "CpG_pos")

methylation_in_hESC <- SummarizedExperiment(
  assays = as.matrix(met_in_hESC[, -(1:2)]),
  rowRanges = methylation_hESC_gr,
  colData = coldata[colnames(met_in_hESC[, -(1:2)]), ])

save(methylation_in_hESC,
     file = "../../eh_data/methylation_in_hESC.rda",
     compress = "xz",
     compression_level = 9)

