## Code to prepare `embryo_sce_Zhu` dataset goes here

library(tidyverse)
library(biomaRt)
library(SingleCellExperiment)


## Data from Single Cell DNA Methylome Sequencing of Human Preimplantation
## Embryos (Zhu et al. 2018)
## 50 FPKM files were downloaded from GEO (accession: GSE81233).


# FPKM expression values import
files <- list.files("../extdata/GSE81233_RAW/", full.names = TRUE )
count_data <- read_tsv(files[1]) %>%
  dplyr::select(gene_id, FPKM)
sample <- gsub(".genes.fpkm_tracking.gz",
               x = gsub('.*_Icm', x = files[1], 'Icm'), '')

# Some genes have several expression values corresponding to several isoforms
# Remove duplicated and keep the max FPKM per gene
count_data <- count_data %>%
  group_by(gene_id) %>%
  summarize(max_FPKM = max(FPKM))
names(count_data) <- c("gene", sample)

for (file in files[-1]){
  tmp <- read_tsv(file) %>%
    dplyr::select(gene_id, FPKM)
  sample <- gsub(".genes.fpkm_tracking.gz",
                 x = gsub('.*_Icm', x = file, 'Icm'), '')
  tmp <- tmp %>%
    group_by(gene_id) %>%
    summarize(max_FPKM = max(FPKM))
  names(tmp) <- c("gene", sample)
  count_data <- dplyr::left_join(count_data, tmp)
}


coldata <- tibble(sample = colnames(count_data[,-1])) %>%
  mutate(embryo = gsub(x = sample, pattern = "-\\w*$", ""))

## Determine gender of blastocysts

# identify Y-linked genes highly expressed
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

attributes_vector <-  c("ensembl_gene_id", "external_gene_name",
                        "ensembl_transcript_id", "external_transcript_name",
                        "chromosome_name", "strand", "transcript_start",
                        "transcript_end", "transcription_start_site",
                        "transcript_length", "transcript_tsl",
                        "transcript_gencode_basic", "transcript_appris",
                        "transcript_mane_select", "transcript_biotype",
                        "transcript_is_canonical" )

transcripts_infos <-  as_tibble(biomaRt::getBM(attributes = attributes_vector,
                                               mart = ensembl))
Y_linked <- transcripts_infos %>%
  filter(external_gene_name %in% count_data$gene) %>%
  filter(chromosome_name == "Y") %>%
  pull(external_gene_name) %>%
  unique()

# Y-linked genes  expressed in the dataset
mat <- as.matrix(count_data[,-1])
rownames(mat) <- count_data$gene
Y_linked_expressed <- Y_linked[Y_linked %in% rownames(mat)[rowSums(mat) > 10]]

# After analysis, RPS4Y1 seems to be the best gene to determine gender
# If the mean expression of RPS4Y1 is higher than 50 FPKM, it's a male sample

embryo_gender <- count_data %>%
  filter(gene %in% c("RPS4Y1")) %>%
  pivot_longer(names_to = "sample", values_to = "FPKM", -gene) %>%
  left_join(coldata) %>%
  group_by(embryo) %>%
  summarize(mean_RPS4Y1 = mean(FPKM)) %>%
  mutate(sex = ifelse(mean_RPS4Y1 > 50, "M", "F")) %>%
  dplyr::select(embryo, sex)

coldata <- coldata %>%
  left_join(embryo_gender) %>%
  mutate(genotype = ifelse(sex == "M", "XY", "XX")) %>%
  mutate(stage = "blastocyst")

embryo_sce_Zhu <- SingleCellExperiment(assays = mat[, coldata$sample],
                                       colData = coldata)

save(embryo_sce_Zhu, file = "../../eh_data/embryo_sce_Zhu.rda",
     compress = "xz",
     compression_level = 9)

