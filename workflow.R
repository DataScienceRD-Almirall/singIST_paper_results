# 1. Append the library folder where packages are located
.libPaths(c("/Workspace/Users/amoruno@almirall.com/my_r_packages/v6/", .libPaths()))

# 2. Install singIST package from git repo
remotes::install_git(
    "https://github.com/DataScienceRD-Almirall/singIST.git",
    git = "external"
)
library(singIST)

# 3. Load human and disease model Seurat objects# Human Seurat
human <- readRDS('/Volumes/users_rd/user_amoruno/my_volume/-human-human_bangert.rds')
# Oxazolone and Imiquimod models Seurat
oxa_imq <- readRDS('/Volumes/users_rd/user_amoruno/my_volume/OXA_IMQ.rds')
oxa_imq <- Seurat::UpdateSeuratObject(oxa_imq)
# Ovalbumin model Seurat
ova <- readRDS('/Volumes/users_rd/user_amoruno/my_volume/OVA.rds')
ova <- Seurat::UpdateSeuratObject(ova)

# 4. Preprocess human Seurat object for its use in singIST: pseudobulk and select cell type clusters
# Description: Homogeneize cell type labels to its corresponding 
# Input: Annotated label cluster
# Output: Consensus cluster
# Homogenize cell type annotation clusters
  homogenize_cell_type <- function(celltype) {
    celltype <- ifelse(grepl("T", celltype, fixed = TRUE), "T-cell", celltype)
    celltype <- ifelse(grepl("KC", celltype, fixed = TRUE), "Keratinocytes", celltype)
    celltype <- ifelse(grepl("Melanocytes", celltype, fixed = TRUE), "Melanocytes", celltype)
    celltype <- ifelse(grepl("DC", celltype, fixed = TRUE), "Dendritic Cells", celltype)
    celltype <- ifelse(grepl("LC", celltype, fixed = TRUE), "Langerhans Cells", celltype)
    return(celltype)
  }
# Update cell type annotation
metadata <- human@meta.data
metadata$CELLTYPE_new <- sapply(metadata$CELLTYPE_new, homogenize_cell_type)
human@meta.data <- metadata
# Filter only HC and baseline of dupi treated AD patients (= Lesional)
human <- subset(human, subset = treatment_new %in% c("HC", "Baseline"))
# Pseudobulk log normalize data
human_pseudobulk <- Seurat::AggregateExpression(human, assays = "RNA", group.by = c("CELLTYPE_new", "Sample_id"), normalization.method = "LogNormalize", return.seurat = TRUE)$RNA$data
# Transpose matrix so rows are samples*cell type and columns genes
human_pseudobulk <- t(human_pseudobulk)
# Remove the cell types as we are not going to model
human_pseudobulk <- human_pseudobulk[!(sub("_.*$", "", rownames(human_pseudobulk)) %in% 
                                        c("MastC-other")), ]

# 5. Define superpathways to be analyzed
# Pathways to be analyzed according to Brunner et.al 2017
pathways <- c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION", "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES", "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
"BIOCARTA_INFLAM_PATHWAY", "BIOCARTA_TH1TH2_PATHWAY", "BIOCARTA_CYTOKINE_PATHWAY",
"BIOCARTA_DC_PATHWAY", "KEGG_JAK_STAT_SIGNALING_PATHWAY", "KEGG_ASTHMA", "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "PID_IL12_STAT4_PATHWAY",
"BIOCARTA_CD40_PATHWAY", "PID_IL4_2PATHWAY", "PID_IL23_PATHWAY", "PID_CXCR3_PATHWAY",
"KEGG_HEMATOPOIETIC_CELL_LINEAGE", "PID_IL2_STAT5_PATHWAY", "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY", "PID_CD8_TCR_DOWNSTREAM_PATHWAY", "REACTOME_SIGNALING_BY_INTERLEUKINS",
"KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
"REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"
# ARCHIVED: Archived Founder gene sets that are referenced by current Hallmarks
# This pathway is not included in analysis since it has been archived by msigdb and no longer updated nor curated
# "REACTOME_IMMUNE_SYSTEM"
)
# Data base source of the pathway
dbsource <- sub("_.*", "", pathways)
# Collection and subcollection is always the same for these pathways only Curated pathways "c2"
collection <- rep("c2", length(dbsource))
subcollection <- rep("CP", length(dbsource))
# Choose the same cell types for all pathways under analysis
celltypes <- c("Dendritic Cells", "Keratinocytes", "Langerhans Cells", "Melanocytes", "T-cell")
# Id of human samples
sample_id <- unique(sub(".*_", "", rownames(human_pseudobulk)))
# Identify class of each sample
sample_class <-  c(rep("AD",5), rep("HC", 4))
# Set the target and the base class
base_class <- "HC"
target_class <- "AD"

# 6. Initialize superpathway objects
# Load Gene Set Collection (gse) from MsigDB for human and HGNC gene symbol
gse <- msigdb::getMsigdb(org = "hs", id = "SYM", version = msigdb::getMsigdbVersions()[1])
# Pathway information
pathways_object <- lapply(seq(1, length(pathways)), function(i) new("pathway", standard_name = pathways[i], dbsource = dbsource[i], collection = collection[i], subcollection = subcollection[i]))
# Superpathway information
superpathway_gene_sets <- lapply(seq(1, length(pathways)), function(i) new("superpathway.gene.sets", pathway_info = pathways_object[[i]], celltypes = celltypes))
# As gene sets per cell type are all the same we retrieve the same gene set for all cell types with setRepeatGeneSets
superpathway_gene_sets <- lapply(seq(1, length(pathways)), function(i, gse_data = gse)setRepeatGeneSets(superpathway_gene_sets[[i]], gse = gse_data))

# 7. Define hyperparameters for singIST workflow
# We will use the same hyperparameters and CV strategy for all superpathways

# Choose the array of quantile combinations that will be tested in the Cross Validation. We choose a wide range of 
# quantiles to test, as a narrow selection might not find the globally optimal model. However, the more quantiles 
# we choose the more computation time it will take. 
quantile_comb_table <- lapply(seq(1, length(pathways)), function(i) as.matrix(RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.10), m = length(slot(superpathway_gene_sets[[i]], "celltypes")), TRUE)))
# All superpathways are for binary outcome 
outcome_type <- rep("binary", length(pathways))
# Number of PLS to test set to 3 for all superpathways. 
number_PLS <- rep(as.integer(3), length(pathways))
# Number of folds for CV set to 1, this will perform a LOOCV. Since the sample size is small (9) a LOOCV is more
# appropiate
folds_CV <- rep(as.integer(1), length(pathways))
# Since we are performing LOOCV the number of repetitions per fold can only be 1
repetition_CV <- rep(as.integer(1), length(pathways))

# 8. Initialize hyperparameter and superpathway input objects for singIST workflow
# Hyperparameters for asmbPLS-DA
hyperparameters <- lapply(seq(1, length(pathways)), function(i) new("hyperparameters", quantile_comb_table = quantile_comb_table[[i]], outcome_type = outcome_type[i], number_PLS = number_PLS[i], folds_CV = folds_CV[i], repetition_CV = repetition_CV[i]))
superpathway_input <- lapply(seq(1, length(pathways)), function(i) new("superpathway.input", superpathway_info = superpathway_gene_sets[[i]], hyperparameters_info = hyperparameters[[i]], pseudobulk_lognorm = human_pseudobulk, sample_id = sample_id,sample_class = sample_class, base_class = base_class, target_class = target_class))

# 9. (step1 singIST) fit and validate optimal asmbPLS-DA 
# enable parallellization with 9 workers
BiocParallel::register(BiocParallel::MulticoreParam(workers = 9,
                        exportglobals = FALSE, progressbar = TRUE),
                        default = TRUE)  
# CIP.GIP_significance_full = FALSE by default. Since we are not focusing the dicussion on all the superpathways, and storing all CIP.GIP 
# information might consume too much memory, we set it the default to CIP.GIP_significance_full = FALSE. We will compute later on the CIP/GIP values
# for the superpathways we are interested in a deeper look
superpathway_fit <- multiple_fitOptimal(superpathway_input, parallel = c(TRUE), npermut = c(10000), type = c("jackknife"), measure = c("F1"), nbObsPermut = c(9)) # Permute all samples (9) per permutation
# disable parallelization
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE) 

# 10. Define objects for the disease models
colnames(slot(ova, "meta.data"))[c(5, 8)] <- c("class", "celltype_cluster")
colnames(slot(oxa_imq, "meta.data"))[c(5, 11)] <- c("class", "celltype_cluster")
# Filter Seurat object by only cell types used in the singIST analysis
ova <- ova[, ova$celltype_cluster %in% c("T cells", "Keratinocytes", "Dendritic cells")]
# Filter OXA and IMQ objects separately
oxa <- oxa_imq[, oxa_imq$class %in% c("OXA", "ETOH")]
imq <- oxa_imq[, oxa_imq$class %in% c("IMQ", "VEH")]
# OVA
## Rename metadata group variable to `class` and cell type annotation to `celltype_cluster` as those are the variable names that the class `mapping.organism` identify
## Celltype_mapping slot left hand side of equality are levels of `celltype_cluster` while left hand side of equality
## are cell types as denoted in the `superpathway.gene.sets` object 
ova_org <- new("mapping.organism", target_class = "OVA", base_class = "SAL", organism = "Mus musculus",
                celltype_mapping = list(
                    "Dendritic Cells" = c("Dendritic cells"),
                    "Keratinocytes" = c("Keratinocytes"),
                    "Langerhans Cells" = c(), # There are < 50 Langerhans cells, hence we consider it as null
                    "Melanocytes" = c(),
                    "T-cell" = c("T cells")),
                counts = ova)
# OXA
oxa_org <- new("mapping.organism", target_class = "OXA", base_class = "ETOH", organism = "Mus musculus",
                celltype_mapping = list(
                    "Dendritic Cells" = c("cDC2", "cDC1", "migratory DCs"),
                    "Keratinocytes" = c("Keratinocytes"),
                    "Langerhans Cells" = c("LC"), 
                    "Melanocytes" = c(),
                    "T-cell" = c("DETC", "dγdT", "T")),
                counts = oxa)
# IMQ
imq_org <- new("mapping.organism", target_class = "IMQ", base_class = "VEH", organism = "Mus musculus",
                celltype_mapping = list(
                    "Dendritic Cells" = c("cDC2", "cDC1", "migratory DCs"),
                    "Keratinocytes" = c("Keratinocytes"),
                    "Langerhans Cells" = c("LC"), 
                    "Melanocytes" = c(),
                    "T-cell" = c("DETC", "dγdT", "T")),
                counts = imq)

# 11. Compute recapitulations
# Compute OXA recapitulations for all superpathways
recapitulations_oxa <- multiple_singISTrecapitulations(oxa_org, superpathway_fit)
# Compute IMQ recapitulations for all superpathways
recapitulations_imq <- multiple_singISTrecapitulations(imq_org, superpathway_fit)
# Compute OVA recapitulations for all superpathways
recapitulations_ova <- multiple_singISTrecapitulations(ova_org, superpathway_fit, logfc.treshold = 0.1)
# Adjust p-values of asmbPLS-DA global significance test with Benjamini-Hochberg
recapitulations_ova$superpathway$p_val <- stats::p.adjust(recapitulations_ova$superpathway$p_val, method = "BH", n = length(recapitulations_ova$superpathway$p_val))
recapitulations_oxa$superpathway$p_val <- stats::p.adjust(recapitulations_oxa$superpathway$p_val, method = "BH", n = length(recapitulations_oxa$superpathway$p_val))
recapitulations_imq$superpathway$p_val <- stats::p.adjust(recapitulations_imq$superpathway$p_val, method = "BH", n = length(recapitulations_imq$superpathway$p_val))
recapitulations_all <- render_multiple_outputs(objects = list(recapitulations_oxa, recapitulations_ova, recapitulations_imq))

# 12. Save optimal models and recapitulation objects
save(recapitulations_all, file = "/Volumes/users_rd/user_amoruno/my_volume/recapitulations_all.RData")
save(superpathway_fit, file = "/Volumes/users_rd/user_amoruno/my_volume/superpathway_fit.RData")
