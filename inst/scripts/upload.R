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
