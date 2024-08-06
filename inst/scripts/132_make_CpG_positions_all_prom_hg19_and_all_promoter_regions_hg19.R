## Code to prepare `CpG_positions_all_prom_hg19` and `all_promoter_regions_hg19`
## intermediate datasets goes here

library(tidyverse)
library(GenomicRanges)
library(liftOver)
library(AnnotationHub)
load("../../eh_data/all_genes.rda")


# Download human genome fasta file

# download.file("http://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz",
#               destfile = "/storage/research/dduv/cbio-lg/cluster/GENOMES/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz")

hg <- rtracklayer::import("../../../CTdata_extdata/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz",
                          format = "fasta")

################################################################################
# Create a Grange with all CpG positions in hg19
################################################################################

## Initiate CT_methylation tibble
## Use whole genome => identify all CpGs positions
CpG_positions <- tibble()
chromosomes <- c(1:22, "X", "Y")
for (chr_number in chromosomes){
  chromosome_ind <- grep(hg@ranges@NAMES,
                         pattern = paste0('dna:chromosome chromosome:GRCh37:',
                                          chr_number,':'))
  # finds the row corresponding to selected chromosome
  chr_seq <- hg[chromosome_ind][[1]] # gets the sequence of selected chromosome
  chr_CG <- Biostrings::matchPattern("CG", chr_seq)
  # gets the coords of CpG dinucleotides on selected chromosome
  tmp <- as_tibble(as.data.frame(chr_CG@ranges))%>%
    mutate(chr = paste0("chr", chr_number)) %>%
    mutate(CpG_pos = start) %>%
    dplyr::rename(Cytosine_fw = start, Cytosine_rev = end) %>%
    dplyr::select(chr, CpG_pos, Cytosine_fw, Cytosine_rev)
  CpG_positions <- rbind(CpG_positions, tmp)
}

CpG_positions_GR <- makeGRangesFromDataFrame(CpG_positions,
                                             keep.extra.columns = TRUE,
                                             seqnames.field = "chr",
                                             start.field = "Cytosine_fw",
                                             end.field = "Cytosine_rev")

################################################################################
# Create a table with all genes TSS in hg19 to lift it to hg19
################################################################################

# download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz",
#               destfile = "../../../CTdata_extdata/hg38ToHg19.over.chain.gz")
# Has to be unziped for next line to work

chain <- rtracklayer::import("../../../CTdata_extdata/hg38ToHg19.over.chain")

all_genes_liftover <- all_genes
all_genes_liftover$chr <- paste0("chr", all_genes_liftover$chr)
all_genes_liftover$strand[all_genes_liftover$strand == 1] <- '+'
all_genes_liftover$strand[all_genes_liftover$strand == -1] <- '-'

all_genes_TSS_liftover_GR <- makeGRangesFromDataFrame(all_genes_liftover,
                                                      keep.extra.columns = TRUE,
                                                      seqnames.field = "chr",
                                                      start.field = "transcription_start_site",
                                                      end.field = "transcription_start_site")

# convert TSS to hg19
all_genes_TSS_liftover <- unlist(liftOver(all_genes_TSS_liftover_GR, chain))

# genes were lost in the liftover procedure: incl 2 CT "FAM230C" and "FLJ36000"
# all_genes$external_gene_name[!all_genes$external_gene_name %in%
#                                all_genes_TSS_liftover$external_gene_name]

liftover_TSS <- tibble(external_gene_name = all_genes_TSS_liftover$external_gene_name,
                       TSS_liftover = all_genes_TSS_liftover@ranges@start)

################################################################################
# Define promoter regions for all genes in hg19
################################################################################

all_genes_hg19 <- all_genes %>%
  left_join(liftover_TSS) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, family, chr, strand,
                transcription_start_site, TSS_liftover, everything())

nt_up <- 1000
nt_down <- 1000

all_promoter_regions_hg19 <- all_genes_hg19 %>%
  mutate(start = case_when(strand == 1 ~ TSS_liftover - nt_up,
                           strand == -1 ~ TSS_liftover - nt_down)) %>%
  mutate(stop = case_when(strand == 1 ~ TSS_liftover + nt_down,
                          strand == -1 ~ TSS_liftover + nt_up)) %>%
  mutate(chromosome_name = paste0("chr", chr)) %>%
  mutate(strand = case_when(strand == 1 ~ "+",
                            strand == -1 ~ "-"))
all_promoter_regions_hg19_GR <- makeGRangesFromDataFrame(all_promoter_regions_hg19 %>%
                                                           filter(!is.na(start)),
                                                         keep.extra.columns = TRUE,
                                                         seqnames.field = "chromosome_name",
                                                         start.field = "start",
                                                         end.field = "stop")
all_promoter_regions_hg19 <- all_promoter_regions_hg19_GR

save(all_promoter_regions_hg19,
     file = "../extdata/all_promoter_regions_hg19.rda")


################################################################################
# Create a GRange with CpGs in promoters
################################################################################

CpG_positions_all_prom_hg19 <- subsetByOverlaps(CpG_positions_GR,
                                                all_promoter_regions_hg19_GR)
save(CpG_positions_all_prom_hg19,
     file = "../extdata/CpG_positions_all_prom_hg19.rda")
