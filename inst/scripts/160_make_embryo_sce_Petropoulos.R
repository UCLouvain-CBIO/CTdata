## Code to prepare `embryo_sce_Zhu` dataset goes here

library(tidyverse)
library(biomaRt)
library(SingleCellExperiment)


## Data from Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in
## Human Preimplantation Embryos (Petropulous et al, 2014)


## Processed data in RPKM was downloaded here
## https://www.ebi.ac.uk/biostudies/files/E-MTAB-3929/rpkm.txt

## Sample metadata was downloaded here
## https://www.ebi.ac.uk/biostudies/files/E-MTAB-3929/E-MTAB-3929.sdrf.txt

counts <- read_tsv("../extdata/Petropoulos_rpkm.txt")

mat <- as.matrix(counts[,-1])
rownames(mat) <- counts$...1

metadata <- read_tsv("../extdata/Petropoulos_metadata.txt")
coldata <- as.data.frame(metadata)
rownames(coldata) <- coldata$`Source Name`

coldata$day <- coldata$`Characteristics[developmental stage]`
coldata$day <- gsub(pattern = "embryonic day ", x = coldata$day,
                    replacement = "E")
coldata$stage <- ifelse(coldata$day == 'E3', "Morula", "Blastocyst")
coldata$stage <- factor(coldata$stage, levels = c("Morula", "Blastocyst"))
coldata$individual <- coldata$`Characteristics[individual]`
coldata$lineage <- coldata$`Characteristics[inferred lineage]`
coldata$singleCellWellQuality <- coldata$`Characteristics[single cell well quality]`
colnames(coldata)[1] <- "sample"
coldata <- coldata[, c(1,46:50)]

embryo_sce_Petropoulos <- SingleCellExperiment(assays = mat[, rownames(coldata)],
                                               colData = coldata)


## Determine embryo sex
## Inference of Embryonic sex : Identify Y linked to use for sex determination
## Genes should be highly expressed in some cells, with a clear cutoff in
## approximatively 50% of sample. Determine sex separately for each stage (day)
## as expression levels of Y-linked genes can be quite different at each stage


# Identify Y genes highly expressed in the dataset
# Y-linked genes highly expressed in the dataset and highly variable
# (should be highly expressed in males but not in females)

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
  filter(external_gene_name %in% rownames(embryo_sce_Petropoulos)) %>%
  filter(chromosome_name == "Y") %>%
  pull(external_gene_name) %>%
  unique()

Y_linked <- Y_linked[Y_linked %in% rownames(embryo_sce_Petropoulos)]

# Remove Y linked genes almost not detected in the dataset
x <- tibble(gene = Y_linked,
            detected_in_percent_cells =
              rowSums(assay(embryo_sce_Petropoulos)[Y_linked, ] > 0) /
              ncol(embryo_sce_Petropoulos),
            sumY = rowSums(assay(embryo_sce_Petropoulos)[Y_linked, ]),
            var = rowVars(assay(embryo_sce_Petropoulos)[Y_linked, ])) %>%
  arrange(desc(sumY)) %>%
  filter(detected_in_percent_cells > 0.1)



# Identify genes expressed with a clear cutoff in approximatively 50% of sample
# Genes expressed in much more than 50% of cells could be pseudoautosomal genes
# also detected in females

# Rapid check using RPS4Y1 gene (that seems to be a very good marker)
# Determine sex based on RPS4Y1 gene only and check the levels of
# other Y-linked genes to see if they are also good markers of male embryos

sex_approximation <- as_tibble(assay(embryo_sce_Petropoulos)["RPS4Y1", ],
                               rownames = "sample") %>%
  mutate(sex_approx = case_when(value > 0 ~ "M",
                                value == 0 ~ "F")) %>%
  arrange(value) %>%
  mutate(rank = seq_along(colnames(embryo_sce_Petropoulos)))

bad_markers <- c("SNORA70", "TMSB4Y", "RPS4Y2", "RBMY2FP")

Y_linked_selected <- x %>%
  filter(!gene %in% bad_markers)

# Used 11 genes for sex determination


## E3 embryos

day <- "E3"
sce <- embryo_sce_Petropoulos[, embryo_sce_Petropoulos$day == day]
sce$sumY <- colSums(assay(sce)[Y_linked_selected$gene, ])
ranks <- as_tibble(colData(sce)) %>%
  arrange(sumY) %>%
  mutate(rank = 1:ncol(sce)) %>%
  mutate(sex = case_when(sumY > 100 ~ "M",
                            sumY < 25 ~ "F",
                            sumY < 100 & sumY > 25 ~ NA))

# 2 embryos have conflicting results (E3.2 and E3.3)

# Flag embryos with conflicting sex results
# (i.e an individual where some cells would be flagged as F and other as M)

ambigous <- ranks %>%
  group_by(individual, sex) %>%
  dplyr::count() %>%
  na.omit() %>%
  ungroup() %>%
  filter(duplicated(individual)) %>%
  pull(individual)
ranks <- ranks %>%
  mutate(ambigous = ifelse(individual %in% ambigous, TRUE, FALSE))

# For non ambigous embryo, replace the NA in sex by the value of other cells
# from the same individual
F_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "F") %>%
  pull(individual)

M_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "M") %>%
  pull(individual)

ranks$sex[ranks$individual %in% F_individual] <- "F"
ranks$sex[ranks$individual %in% M_individual] <- "M"

individual_sex <- ranks

## E4 embryos

day <- "E4"
sce <- embryo_sce_Petropoulos[, embryo_sce_Petropoulos$day == day]
sce$sumY <- colSums(assay(sce)[Y_linked, ])
ranks <- as_tibble(colData(sce)) %>%
  arrange(sumY) %>%
  mutate(rank = 1:ncol(sce)) %>%
  mutate(sex = case_when(sumY > 300 ~ "M",
                            sumY < 100 ~ "F",
                            sumY < 300 & sumY > 100 ~ NA) )

# Flag embryos with conflicting sex results
# (i.e an individual where some cells would be flagged as F and other as M)
ambigous <- ranks %>%
  group_by(individual, sex) %>%
  dplyr::count() %>%
  na.omit() %>% ungroup() %>%
  filter(duplicated(individual)) %>%
  pull(individual)
# No ambiguous

ranks <- ranks %>%
  mutate(ambigous = ifelse(individual %in% ambigous, TRUE, FALSE))

# For non ambigous embryo, replace the NA in sex by the value of other cells
# from the same individual
F_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "F") %>%
  pull(individual)

M_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "M") %>%
  pull(individual)

ranks$sex[ranks$individual %in% F_individual] <- "F"
ranks$sex[ranks$individual %in% M_individual] <- "M"

individual_sex <- rbind(individual_sex, ranks)


## E5 embryos

day <- "E5"
sce <- embryo_sce_Petropoulos[, embryo_sce_Petropoulos$day == day]
sce$sumY <- colSums(assay(sce)[Y_linked, ])
ranks <- as_tibble(colData(sce)) %>%
  arrange(sumY) %>%
  mutate(rank = 1:ncol(sce)) %>%
  mutate(sex = case_when(sumY > 300 ~ "M",
                            sumY < 100 ~ "F",
                            sumY < 300 & sumY > 100 ~ NA) )

# Flag embryos with conflicting sex results
# (i.e an individual where some cells would be flagged as F and other as M)
ambigous <- ranks %>%
  group_by(individual, sex) %>%
  dplyr::count() %>%
  na.omit() %>% ungroup() %>%
  filter(duplicated(individual)) %>%
  pull(individual)

ranks <- ranks %>%
  mutate(ambigous = ifelse(individual %in% ambigous, TRUE, FALSE))

# For non ambigous embryo, replace the NA in sex by the value of other cells
# from the same individual
F_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "F") %>%
  pull(individual)

M_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "M") %>%
  pull(individual)

ranks$sex[ranks$individual %in% F_individual] <- "F"
ranks$sex[ranks$individual %in% M_individual] <- "M"

individual_sex <- rbind(individual_sex, ranks)


## E6 embryos

day <- "E6"
sce <- embryo_sce_Petropoulos[, embryo_sce_Petropoulos$day == day]
sce$sumY <- colSums(assay(sce)[Y_linked, ])
ranks <- as_tibble(colData(sce)) %>%
  arrange(sumY) %>%
  mutate(rank = 1:ncol(sce)) %>%
  mutate(sex = case_when(sumY > 400 ~ "M",
                            sumY < 200 ~ "F",
                            sumY < 400 & sumY > 200 ~ NA))

# Flag embryos with conflicting sex results
# (i.e an individual where some cells would be flagged as F and other as M)
ambigous <- ranks %>%
  group_by(individual, sex) %>%
  dplyr::count() %>%
  na.omit() %>% ungroup() %>%
  filter(duplicated(individual)) %>%
  pull(individual)
ranks <- ranks %>%
  mutate(ambigous = ifelse(individual %in% ambigous, TRUE, FALSE))

# For non ambigous embryo, replace the NA in sex by the value of other cells
# from the same individual
F_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "F") %>%
  pull(individual)

M_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "M") %>%
  pull(individual)

ranks$sex[ranks$individual %in% F_individual] <- "F"
ranks$sex[ranks$individual %in% M_individual] <- "M"

individual_sex <- rbind(individual_sex, ranks)

## E7 embryos

day <- "E7"
sce <- embryo_sce_Petropoulos[, embryo_sce_Petropoulos$day == day]
sce$sumY <- colSums(assay(sce)[Y_linked, ])
ranks <- as_tibble(colData(sce)) %>%
  arrange(sumY) %>%
  mutate(rank = 1:ncol(sce)) %>%
  mutate(sex = case_when(sumY > 400 ~ "M",
                            sumY < 200 ~ "F",
                            sumY < 400 & sumY > 200 ~ NA))

# Flag embryos with conflicting sex results
# (i.e an individual where some cells would be flagged as F and other as M)
ambigous <- ranks %>%
  group_by(individual, sex) %>%
  dplyr::count() %>%
  na.omit() %>% ungroup() %>%
  filter(duplicated(individual)) %>%
  pull(individual)
ranks <- ranks %>%
  mutate(ambigous = ifelse(individual %in% ambigous, TRUE, FALSE))

# For non ambigous embryo, replace the NA in sex by the value of other cells
# from the same individual
F_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "F") %>%
  pull(individual)

M_individual <- ranks %>%
  filter(!ambigous) %>%
  dplyr::select(individual, sex) %>%
  filter(!is.na(sex)) %>%
  unique() %>%
  filter(sex == "M") %>%
  pull(individual)

ranks$sex[ranks$individual %in% F_individual] <- "F"
ranks$sex[ranks$individual %in% M_individual] <- "M"

individual_sex <- rbind(individual_sex, ranks)




## Summary of sexs

# summarise sex info in a table,
# flag individuals with conflicting sex determination
cell_sex <- individual_sex %>%
  dplyr::select(sample, individual, sex, ambigous)

coldata <- merge(colData(embryo_sce_Petropoulos), cell_sex)
rownames(coldata) <- coldata$sample
coldata <- coldata[rownames(colData(embryo_sce_Petropoulos)),]
colData(embryo_sce_Petropoulos) <- coldata


### ChrX ratios

X_linked <- transcripts_infos %>%
  filter(external_gene_name %in% rownames(embryo_sce_Petropoulos)) %>%
  filter(chromosome_name == "X") %>%
  pull(external_gene_name) %>%
  unique()

autosomal_genes <- transcripts_infos %>%
  filter(external_gene_name %in% rownames(embryo_sce_Petropoulos)) %>%
  filter(chromosome_name %in% c(1:22)) %>%
  pull(external_gene_name) %>%
  unique()

sce <- embryo_sce_Petropoulos

RPKM_sum_by_chr <- tibble(sample = colnames(sce),
                          chrY_RPKM_sum = colSums(assay(sce)[Y_linked,]),
                          chrX_RPKM_sum = colSums(assay(sce)[X_linked,]),
                          autosomal_RPKM_sum = colSums(assay(sce)[autosomal_genes,]),
                          chrY_RPKM_mean = colMeans(assay(sce)[Y_linked,]),
                          chrX_RPKM_mean = colMeans(assay(sce)[X_linked,]),
                          autosomal_RPKM_mean = colMeans(assay(sce)[autosomal_genes,])) %>%
  mutate(ratio_X_A = chrX_RPKM_mean / autosomal_RPKM_mean) %>%
  mutate(ratio_Y_A = chrY_RPKM_mean / autosomal_RPKM_mean) %>%
  left_join(as_tibble(colData(sce)))

RPKM_sum_by_chr$sex <- factor(RPKM_sum_by_chr$sex, levels = c("F", "M"))

ordered_samples <- RPKM_sum_by_chr %>%
  arrange(sex, day) %>%
  pull(sample)

# "E5.early.31" sample flagged as Xo (Turner) in Petropoulos
# "E6.15" sample flagged as embryo with biallelic expression of many X_linked
# genes in Petropoulos

bad_embryos <- c("E5.early.31", "E6.15")
embryo_sce_Petropoulos <-
  embryo_sce_Petropoulos[, !embryo_sce_Petropoulos$individual %in% bad_embryos]

## Complete colData
embryo_sce_Petropoulos$genotype <- embryo_sce_Petropoulos$sex
embryo_sce_Petropoulos$genotype[!is.na(embryo_sce_Petropoulos$genotype) &
                           embryo_sce_Petropoulos$sex == "F"] <- "XX"
embryo_sce_Petropoulos$genotype[!is.na(embryo_sce_Petropoulos$genotype) &
                           embryo_sce_Petropoulos$genotype == "M"] <- "XY"
embryo_sce_Petropoulos$day_stage[embryo_sce_Petropoulos$day == "E3"]  <-
  "Morula E3"
embryo_sce_Petropoulos$day_stage[embryo_sce_Petropoulos$day == "E4"] <-
  "Blastocyst E4"
embryo_sce_Petropoulos$day_stage[embryo_sce_Petropoulos$day == "E5"] <-
  "Blastocyst E5"
embryo_sce_Petropoulos$day_stage[embryo_sce_Petropoulos$day == "E6"] <-
  "Blastocyst E6"
embryo_sce_Petropoulos$day_stage[embryo_sce_Petropoulos$day == "E7"] <-
  "Blastocyst E7"
embryo_sce_Petropoulos$day_stage <- factor(embryo_sce_Petropoulos$day_stage,
                                    levels = c("Morula E3", "Blastocyst E4",
                                               "Blastocyst E5", "Blastocyst E6",
                                               "Blastocyst E7"))

save(embryo_sce_Petropoulos, file = "../../eh_data/embryo_sce_Petropoulos.rda",
     compress = "xz",
     compression_level = 9)

