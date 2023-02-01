## BiocManager::install("AzureStor")
library("AzureStor")

sas <- "taken_comes_here"
url <- "https://bioconductorhubs.blob.core.windows.net"
ep <- storage_endpoint(url, sas = sas)
container <- storage_container(ep, "staginghub")

src <- dir("/data/cbio/Packages/CTdata/eh_data",
           recursive=TRUE, full.names=TRUE)

files <- dir("/data/cbio/Packages/CTdata/eh_data", recursive=TRUE)
dest <- paste0("CTdata/eh_data/", files)

storage_multiupload(container, src = src, dest = dest)
