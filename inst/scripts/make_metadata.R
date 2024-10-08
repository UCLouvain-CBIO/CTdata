metadata <-
    data.frame(
        Title = c(
            "GTEX_data",
            "CCLE_data",
            "normal_tissues_multimapping_data",
            "DAC_treated_cells",
            "DAC_treated_cells_multimapping",
            "TCGA_TPM",
            "CT_methylation_in_tissues",
            "CT_mean_methylation_in_tissues",
            "TCGA_CT_methylation",
            "CT_genes",
            "CCLE_correlation_matrix",
            "testis_sce",
            "scRNAseq_HPA",
            "all_genes",
            "FGC_sce",
            "HPA_cell_type_specificities",
            "mean_methylation_in_tissues",
            "methylation_in_tissues",
            "oocytes_sce",
            "TCGA_methylation",
            "embryo_sce_Zhu",
            "embryo_sce_Petropoulos",
            "hESC_data",
            "methylation_in_hESC",
            "mean_methylation_in_hESC",
            "methylation_in_FGC",
            "mean_methylation_in_FGC",
            "methylation_in_embryo",
            "mean_methylation_in_embryo"),
        Description = c(
            "Gene expression data in normal tissues from GTEx database",
            "Gene expression data in cancer cell lines from CCLE",
            "Gene expression values in normal tissues with or without allowing multimapping",
            "Gene expression values in a set of cell lines treated or not with 5-Aza-2'-Deoxycytidine (DAC), a demethylating agent",
            "Gene expression values (multimapping allowed) in a set of cell lines treated or not with 5-Aza-2'-Deoxycytidine (DAC)",
            "Gene expression data in TCGA samples",
            "Methylation values of CpGs located within Cancer-Testis promoters in a set of normal tissues",
            "Mean methylation values of all CpGs located within Cancer-Testis (CT) promoters in a set of normal tissues",
            "Methylation of CT promoters in TCGA samples",
            "Cancer-Testis (CT) genes description",
            "Gene correlations in CCLE cancer cell lines",
            "Testis scRNAseq data",
            "Gene expression profiles in different human cell types based on scRNAseq data from The Human Protein Atlas",
            "All genes description, according to the analysis done for CT genes",
            "Fetal gonad scRNAseq",
            "Cell type specificities (from HPA)",
            "Mean methylation values of all CpGs located within all (CT) promoters in a set of normal tissues",
            "Methylation values of CpGs located within all promoters in a set of normal tissues",
            "Oocytes scRNAseq",
            "Methylation of all genes promoters in TCGA samples",
            "Gene expression values in early embryo (blastocyst)",
            "Gene expression values in early embryo (blastocyst and Morula)",
            "Gene expression data in human embryonic stem cells",
            "Methylation values of CpGs located within all promoters in human embryonic stem cells",
            "Mean methylation values of all CpGs located within all promoters in human embryonic stem cells",
            "Methylation values of CpGs located within all promoters in fetal germ cells",
            "Mean methylation values of all CpGs located within all promoters in fetal germ cells",
            "Methylation values of CpGs located within all promoters in early embryo",
            "Mean methylation values of all CpGs located within all promoters in early embryo"
            ),
        BiocVersion = c(
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.17",
            "3.17",
            "3.17",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19",
            "3.19"),
        Genome = c(rep("hg38", 25), rep ("hg19", 4)),
        SourceType = c(
            "tab",
            "CSV",
            "FASTQ",
            "FASTQ",
            "FASTQ",
            "TSV",
            "BED",
            "BED",
            "TXT",
            "TSV",
            "CSV",
            "FASTQ",
            "TSV",
            "TSV",
            "RDS",
            "TSV",
            "BED",
            "BED",
            "TXT",
            "TSV",
            "tab",
            "TXT",
            "FASTQ",
            "BED",
            "BED",
            "TXT",
            "TXT",
            "BED",
            "BED"),
        SourceUrl = c(
            "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
            "https://ndownloader.figshare.com/files/34989922, https://ndownloader.figshare.com/files/35020903",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://portal.gdc.cancer.gov/",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://portal.gdc.cancer.gov/",
            "https://zenodo.org/record/7537824/files/cancermine_collated.tsv?download=1",
            "https://ndownloader.figshare.com/files/34989922, https://ndownloader.figshare.com/files/35020903",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120508",
            "https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
            "https://zenodo.org/record/7537824/files/cancermine_collated.tsv?download=1",
            "https://cellxgene.cziscience.com/collections/661a402a-2a5a-4c71-9b05-b346c57bc451Data",
            "https://www.proteinatlas.org/",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154762",
            "https://portal.gdc.cancer.gov/",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81233",
            "https://www.ebi.ac.uk/biostudies/files/E-MTAB-3929/",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://www.encodeproject.org/",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107714",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107714",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81233",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81233"),
        SourceVersion = c(
            "8",
            "22Q2",
            "5",
            "5",
            "5",
            "36.0",
            "5",
            "5",
            "36.0",
            "49",
            "22Q2",
            "1",
            "1",
            "49",
            "1",
            "23",
            "5",
            "5",
            "1",
            "36.0",
            "1",
            "1",
            "5",
            "5",
            "5",
            "1",
            "1",
            "1",
            "1"),
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = "TRUE",
        DataProvider = c(
            "The Genotype-Tissue Expression",
            "Cancer Cell Line Encyclopedia",
            "ENCODE Project",
            "ENCODE Project",
            "ENCODE Project",
            "The Cancer Genome Atlas",
            "ENCODE Project",
            "ENCODE Project",
            "The Cancer Genome Atlas",
            "CancerMine",
            "Cancer Cell Line Encyclopedia",
            "The adult human testis transcriptional cell atlas, Guo et al. 2018",
            "The Human Protein Atlas",
            "CancerMine",
            "Single-cell roadmap of human gonadal development (Garcia-Alonso, Nature 2022)",
            "The Human Protein Atlas",
            "ENCODE Project",
            "ENCODE Project",
            "NCBI GEO : Decoding dynamic epigenetic landscapes in human oocytes
            using single-cell multi-omics sequencing (Yan et al. Cell Stem Cell 2021)",
            "The Cancer Genome Atlas",
            "Single Cell DNA Methylome Sequencing of Human Preimplantation Embryos (Zhu et al. 2018)",
            "Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in
            Human Preimplantation Embryos (Petropulous et al, 2014)",
            "ENCODE Project",
            "ENCODE Project",
            "ENCODE Project",
            "Data from Dissecting the epigenomic dynamics of human fetal germ cell
            development at single-cell resolution (Li et al. 2021)",
            "Data from Dissecting the epigenomic dynamics of human fetal germ cell
            development at single-cell resolution (Li et al. 2021)",
            "Single Cell DNA Methylome Sequencing of Human Preimplantation Embryos (Zhu et al. 2018)",
            "Single Cell DNA Methylome Sequencing of Human Preimplantation Embryos (Zhu et al. 2018)"),
        Maintainer = "Axelle Loriot <axelle.loriot@uclouvain.be>",
        RDataClass = c(
            "SummarizedExperiment",
            "SummarizedExperiment",
            "SummarizedExperiment",
            "SummarizedExperiment",
            "SummarizedExperiment",
            "SummarizedExperiment",
            "RangedSummarizedExperiment",
            "SummarizedExperiment",
            "RangedSummarizedExperiment",
            "tibble",
            "matrix",
            "SingleCellExperiment",
            "SingleCellExperiment",
            "tibble",
            "SingleCellExperiment",
            "SingleCellExperiment",
            "SummarizedExperiment",
            "RangedSummarizedExperiment",
            "SingleCellExperiment",
            "SummarizedExperiment",
            "SingleCellExperiment",
            "SingleCellExperiment",
            "SummarizedExperiment",
            "RangedSummarizedExperiment",
            "SummarizedExperiment",
            "RangedSummarizedExperiment",
            "RangedSummarizedExperiment",
            "RangedSummarizedExperiment",
            "RangedSummarizedExperiment"),
        DispatchClass = "Rda",
        RDataPath = c(
            "CTdata/eh_data/v2/GTEX_data.rda",
            "CTdata/eh_data/v2/CCLE_data.rda",
            "CTdata/eh_data/v2/normal_tissues_multimapping_data.rda",
            "CTdata/eh_data/v2/DAC_treated_cells.rda",
            "CTdata/eh_data/v2/DAC_treated_cells_multimapping.rda",
            "CTdata/eh_data/v2/TCGA_TPM.rda",
            "CTdata/eh_data/v2/CT_methylation_in_tissues.rda",
            "CTdata/eh_data/v2/CT_mean_methylation_in_tissues.rda",
            "CTdata/eh_data/v2/TCGA_CT_methylation.rda",
            "CTdata/eh_data/v3/CT_genes.rda",
            "CTdata/eh_data/v3/CCLE_correlation_matrix.rda",
            "CTdata/eh_data/v2/testis_sce.rda",
            "CTdata/eh_data/v2/scRNAseq_HPA.rda",
            "CTdata/eh_data/all_genes.rda",
            "CTdata/eh_data/FGC_sce.rda",
            "CTdata/eh_data/HPA_cell_type_specificities.rda",
            "CTdata/eh_data/mean_methylation_in_tissues.rda",
            "CTdata/eh_data/methylation_in_tissues.rda",
            "CTdata/eh_data/oocytes_sce.rda",
            "CTdata/eh_data/TCGA_methylation.rda",
            "CTdata/eh_data/embryo_sce_Zhu.rda",
            "CTdata/eh_data/embryo_sce_Petropoulos.rda",
            "CTdata/eh_data/hESC_data.rda",
            "CTdata/eh_data/methylation_in_hESC.rda",
            "CTdata/eh_data/mean_methylation_in_hESC.rda",
            "CTdata/eh_data/methylation_in_FGC.rda",
            "CTdata/eh_data/mean_methylation_in_FGC.rda",
            "CTdata/eh_data/methylation_in_embryo.rda",
            "CTdata/eh_data/mean_methylation_in_embryo.rda"),
        Tags = c("GeneExpression:NormalTissues",
                 "GeneExpression:CancerData:CellLines",
                 "GeneExpression:NormalTissues:Multimapping",
                 "GeneExpression:CellLines:DemethylatingAgent",
                 "GeneExpression:CellLines:DemethylatingAgent:Multimapping",
                 "GeneExpression:CancerData:TumorSamples",
                 "MethylationData:CancerTestis:NormalTissues",
                 "MethylationData:CancerTestis:NormalTissues",
                 "MethylationData:CancerTestis:CancerData:TumorSamples",
                 "CancerTestis:MethylationData:GeneExpression:CancerData:NormalTissues:TumorSamples",
                 "CancerTestis:GeneExpression:CancerData:CellLines",
                 "GeneExpression:scRNAseq:testis",
                 "GeneExpression:scRNAseq:normalCellTypes",
                 "MethylationData:GeneExpression:CancerData:NormalTissues:TumorSamples",
                 "GeneExpression:scRNAseq:FetalCells:Gonad",
                 "GeneActivation:scRNAseq:testis:normalCellTypes",
                 "MethylationData:NormalTissues",
                 "MethylationData:NormalTissues",
                 "GeneExpression:scRNAseq:Oocyts:Gonad",
                 "GeneExpression:CancerData:TumorSamples",
                 "GeneExpression:scRNAseq:embryo",
                 "GeneExpression:scRNAseq:embryo",
                 "GeneExpression:hESC",
                 "MethylationData:hESC",
                 "MethylationData:hESC",
                 "MethylationData:FetalGermCells",
                 "MethylationData:FetalGermCells",
                 "MethylationData:embryo",
                 "MethylationData:embryo")
    )



## Add default tags
metadata$Tags <- CTdata:::makeTags(metadata$Tags)

write.csv(metadata, file = "../extdata/metadata.csv", row.names = FALSE)

pathToPkg <- "~/dev/" ## change this accordingly

ExperimentHubData::makeExperimentHubMetadata(
                       pathToPackage = file.path(pathToPkg, "CTdata"),
                       fileName = "metadata.csv")
