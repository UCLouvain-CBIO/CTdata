## BiocManager::install("AzureStor")
library("AzureStor")

sas <- "taken_comes_here"
url <- "https://bioconductorhubs.blob.core.windows.net"
ep <- storage_endpoint(url, sas = sas)
container <- storage_container(ep, "staginghub")



## ------------------------------------------
## 2023-07-10

src <- c("/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CT_methylation_in_tissues.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CT_mean_methylation_in_tissues.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/TCGA_CT_methylation.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CT_genes.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CCLE_correlation_matrix.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/testis_sce.rda",
  "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/scRNAseq_HPA.rda")

dest <- c("CTdata/eh_data/v2/CT_methylation_in_tissues.rda",
         "CTdata/eh_data/v2/CT_mean_methylation_in_tissues.rda",
         "CTdata/eh_data/v2/TCGA_CT_methylation.rda",
         "CTdata/eh_data/v2/CT_genes.rda",
         "CTdata/eh_data/v2/CCLE_correlation_matrix.rda",
         "CTdata/eh_data/testis_sce.rda",
         "CTdata/eh_data/scRNAseq_HPA.rda")

storage_multiupload(container, src = src, dest = dest)

## ------------------------------------------
## 2023-10-11

src <- "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CT_methylation_in_tissues.rda"
dest <- "CTdata/eh_data/v3/CT_methylation_in_tissues.rda"

storage_multiupload(container, src = src, dest = dest)

## ------------------------------------------
## 2024-05

src <- c("/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/all_genes.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CCLE_correlation_matrix.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CCLE_data.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/CT_genes.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/DAC_treated_cells_multimapping.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/DAC_treated_cells.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/FGC_sce.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/GTEX_data.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/HPA_cell_type_specificities.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/mean_methylation_in_tissues.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/methylation_in_tissues.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/normal_tissues_multimapping_data.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/oocytes_sce.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/scRNAseq_HPA.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/TCGA_methylation.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/TCGA_TPM.rda",
         "/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/testis_sce.rda")

dest <- c("CTdata/eh_data/v2/CCLE_data.rda",
          "CTdata/eh_data/v2/DAC_treated_cells_multimapping.rda",
          "CTdata/eh_data/v2/DAC_treated_cells.rda",
          "CTdata/eh_data/v2/GTEX_data.rda",
          "CTdata/eh_data/v2/normal_tissues_multimapping_data.rda",
          "CTdata/eh_data/v2/scRNAseq_HPA.rda",
          "CTdata/eh_data/v2/TCGA_TPM.rda",
          "CTdata/eh_data/v2/testis_sce.rda",
          "CTdata/eh_data/v3/CT_genes.rda",
          "CTdata/eh_data/v3/CCLE_correlation_matrix.rda",
          "CTdata/eh_data/all_genes.rda",
          "CTdata/eh_data/FGC_sce.rda",
          "CTdata/eh_data/HPA_cell_type_specificities.rda",
          "CTdata/eh_data/mean_methylation_in_tissues.rda",
          "CTdata/eh_data/methylation_in_tissues.rda",
          "CTdata/eh_data/oocytes_sce.rda",
          "CTdata/eh_data/TCGA_methylation.rda")

# Removed CT_methylation_in_tissues.rda, CT_mean_methylation_in_tissues.rda,
# TCGA_CT_methylation.rda

storage_multiupload(container, src = src, dest = dest)
