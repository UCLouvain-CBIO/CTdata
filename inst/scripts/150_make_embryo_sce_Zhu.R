## Code to prepare `embryo_sce_Zhu` dataset goes here

library(tidyverse)


## Data from Single Cell DNA Methylome Sequencing of Human Preimplantation
## Embryos (Zhu et al. 2018)
## 50 FPKM files wer downloaded from GEO (accession: GSE81233).


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




## Determine gender of blastocysts


# identify Y-linked genes highly expressed
load("~/cluster/CBIO-templates/ensembl_2_geneName/transcripts_infos.rda")
Y_linked <- transcripts_infos %>%
  filter(external_gene_name %in% count_data$gene) %>%
  filter(chromosome_name == "Y") %>%
  pull(external_gene_name) %>% unique()
# Y-linked genes  expressed in the dataset
mat <- as.matrix(count_data[,-1])
rownames(mat) <- count_data$gene
Y_linked_expressed <- Y_linked[Y_linked %in% rownames(mat)[rowSums(mat) > 10]]
count_data %>%
  filter(gene %in% Y_linked_expressed) %>%
  pivot_longer(names_to = "sample", values_to = "FPKM", -gene) %>%
  ggplot(aes(x = sample, y = FPKM)) +
  geom_col() +
  facet_wrap(~ gene)
# RPS4Y1 seems the best gene to determine gender

coldata <- tibble(sample = colnames(count_data[,-1])) %>%
  mutate(embryo = gsub(x = sample, pattern = "-\\w*$", ""))
coldata


count_data %>%
  filter(gene %in% c("RPS4Y1")) %>%
  pivot_longer(names_to = "sample", values_to = "FPKM", -gene) %>%
  left_join(coldata) %>%
  ggplot(aes(x = sample, y = FPKM)) +
  geom_col() +
  facet_wrap(~ embryo, scales = "free_x", ncol = 8) +
  theme(axis.text.x = element_text(size = 6, angle = 90))

embryo_gender <- count_data %>%
  filter(gene %in% c("RPS4Y1")) %>%
  pivot_longer(names_to = "sample", values_to = "FPKM", -gene) %>%
  left_join(coldata) %>%
  group_by(embryo) %>%
  summarize(mean_RPS4Y1 = mean(FPKM)) %>%
  mutate(gender = ifelse(mean_RPS4Y1 > 50, "M", "F")) %>%
  dplyr::select(embryo, gender)

coldata <- coldata %>%
  left_join(embryo_gender)




# save sce


embryo_sce <- SingleCellExperiment(assays = mat[, coldata$sample],
                                   colData = coldata)
embryo_sce

save(embryo_sce, file = "../data/embryo_sce_Zhu_datset.rda")
