# Load libraries
#Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT="TRUE")
remotes::install_github("DataScienceRD-Almirall/singIST")
library(singIST)
library(MASS)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(ggplot2)
library(dplyr)
library(broom)

simulate_pseudobulk <- function(p,
                                cell_types = c("Tcell","Bcell","Mono","NK","DC"),
                                n_case = 5,
                                n_ctrl = 5,
                                cor = TRUE,
                                effect_size = 20,
                                log_norm = TRUE) {
    # total samples
    N <- n_case + n_ctrl
    stopifnot(length(cell_types) >= 1)
    # get p real HGNC symbols
    all_genes <- AnnotationDbi::keys(org.Hs.eg.db, keytype="SYMBOL")
    stopifnot(p <= length(all_genes))
    gene_names <- sample(all_genes, p)
    # prepare effect vector
    if (length(effect_size) == 1) {
        eff <- rep(effect_size, p)
    } else if (length(effect_size) == p) {
        eff <- effect_size
    } else {
        stop("effect_size must be length 1 or length p")
    }
    # build covariance
    if (cor) {
        S <- matrix(runif(p^2, .6, .9), p, p)
        S <- (S + t(S)) / 2
        diag(S) <- 1
        S <- as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
    } else {
        S <- diag(p)
    }
    # labels
    labels <- factor(c(rep("case", n_case), rep("ctrl", n_ctrl)),
                     levels = c("case","ctrl"))
    case_idx <- which(labels == "case")
    
    # simulate each block
    blocks <- lapply(cell_types, function(ct, labels_ = labels) {
        # latent log-expression
        Z_ctrl <- MASS::mvrnorm(n = n_ctrl, mu = rep(10000, p), Sigma = S)
        Z_case <- MASS::mvrnorm(n = n_case, mu = eff, Sigma = S)
        Z <- rbind(Z_case, Z_ctrl)
        
        # re-log normalize to mimic log1p(CPM) truncate to 0 if negative
        if(log_norm){
            X <- abs(log2(ifelse(Z <= -1, 0, Z) + 1))
        }else{
            X <- abs(Z)
        }
        
        rownames(X) <- paste0(ct, "_", labels_, seq_len(N))
        colnames(X) <- gene_names
        X
    })
    names(blocks) <- cell_types
    
    list(blocks = do.call(rbind, blocks),
         labels = labels,
         genes  = gene_names)
}

# Usage examples:
# 1) constant shift of +3 on all genes:
#    sim1 <- simulate_pseudobulk(p = 100, effect_size = 3)
#
# 2) varying shifts per gene:
#    shifts <- runif(100, 0.5, 2)
#    sim2   <- simulate_pseudobulk(p = 100, effect_size = shifts)

# Simulate a list of fold-change data.frames for use as FC_list in biological_link_function
#
# cell_types  : character vector of cell‐type block names
# genes       : character vector of gene names to include in each block
# fc_min      : numeric, minimum of uniform distribution for avg_log2FC
# fc_max      : numeric, maximum of uniform distribution for avg_log2FC
# prop_signif : numeric in [0,1], proportion of genes with p_val_adj ≤ 0.05
#
# Returns a named list of data.frames, one per cell type, each with columns:
#   p_val, avg_log2FC, pct.1, pct.2, p_val_adj
simulate_FC_list <- function(cell_types,
                             genes,
                             fc_min = -1,
                             fc_max = 1,
                             prop_signif = 0.5) {
    stopifnot(is.character(cell_types),
              is.character(genes),
              is.numeric(fc_min), length(fc_min) == 1,
              is.numeric(fc_max), length(fc_max) == 1,
              fc_min <= fc_max,
              is.numeric(prop_signif),
              prop_signif >= 0, prop_signif <= 1)
    n_genes <- length(genes)
    n_signif <- floor(prop_signif * n_genes)
    FC_list <- setNames(vector("list", length(cell_types)), cell_types)
    for (ct in cell_types) {
        # draw avg_log2FC for each gene uniformly between fc_min and fc_max
        avg_log2FC <- runif(n_genes, min = fc_min, max = fc_max)
        
        # generate p-values and adjust to achieve prop_signif
        p_val     <- runif(n_genes, 0, 1)
        p_val_adj <- numeric(n_genes)
        sig_idx   <- sample.int(n_genes, n_signif)
        if (n_signif > 0) {
            p_val_adj[sig_idx] <- runif(n_signif, 0, 0.05)
        }
        if (n_signif < n_genes & n_signif > 0) {
            p_val_adj[-sig_idx] <- runif(n_genes - n_signif, 0.05, 1)
        }
        if(n_signif == 0){
            p_val_adj[1:n_genes] <- 1
        }
        
        # generate pct.1 and pct.2 (percent expressed in case/control)
        pct.1 <- runif(n_genes, 0, 1)
        pct.2 <- runif(n_genes, 0, 1)
        
        df <- data.frame(
            p_val      = p_val,
            avg_log2FC = avg_log2FC,
            pct.1      = pct.1,
            pct.2      = pct.2,
            p_val_adj  = p_val_adj,
            row.names  = genes,
            stringsAsFactors = FALSE
        )
        
        FC_list[[ct]] <- df
    }
    
    FC_list
}

# Adjust a Seurat object to have exactly new_genes, whatever they are
# - If new_gene exists in the original object, copy its values
# - If not, fill with zeros
# Clear reductions because they become invalid after changing features

adjust_seurat_features <- function(obj, new_genes) {
    require(Seurat)
    require(Matrix)
    
    # Old matrix
    counts_old <- obj@assays$RNA@counts
    old_feats  <- rownames(counts_old)
    cells      <- colnames(counts_old)
    
    # Prepares new matrix: rows = new_genes, columns = same cell types
    n_new <- length(new_genes)
    n_cells <- length(cells)
    counts_new <- Matrix(0, nrow = n_new, ncol = n_cells, sparse = TRUE)
    rownames(counts_new) <- new_genes
    colnames(counts_new) <- cells
    
    # For each new_gene: if it exists copy the row, otherwise keep it 0
    common <- intersect(new_genes, old_feats)
    if (length(common)>0) {
        idx_old <- match(common, old_feats)
        idx_new <- match(common, new_genes)
        counts_new[idx_new, ] <- counts_old[idx_old, ]
    }
    obj@assays$RNA@counts <- counts_new
    
    # Same logic for data (normalization) and scale.data if they exist
    if (!is.null(obj@assays$RNA@data)) {
        data_old <- obj@assays$RNA@data
        data_new <- Matrix(0, nrow = n_new, ncol = n_cells, sparse = TRUE)
        rownames(data_new) <- new_genes; colnames(data_new) <- cells
        common <- intersect(new_genes, rownames(data_old))
        if (length(common)>0) {
            data_new[match(common,new_genes), ] <-
                data_old[match(common,rownames(data_old)), ]
        }
        obj@assays$RNA@data <- data_new
    }
    if (!is.null(obj@assays$RNA@scale.data)) {
        scale_old <- obj@assays$RNA@scale.data
        scale_new <- matrix(0, nrow = n_new, ncol = n_cells)
        rownames(scale_new) <- new_genes; colnames(scale_new) <- cells
        common <- intersect(new_genes, rownames(scale_old))
        if (length(common)>0) {
            scale_new[match(common,new_genes), ] <-
                scale_old[match(common,rownames(scale_old)), ]
        }
        obj@assays$RNA@scale.data <- scale_new
    }
    
    # Meta.features: creates new files with NAs, copies the old ones wherever it exists
    meta_old <- obj@assays$RNA@meta.features
    meta_new <- data.frame(row.names = new_genes,
                           matrix(NA, nrow = n_new, ncol = ncol(meta_old)))
    colnames(meta_new) <- colnames(meta_old)
    common <- intersect(new_genes, rownames(meta_old))
    if (length(common)>0) {
        meta_new[common, ] <- meta_old[common, ]
    }
    obj@assays$RNA@meta.features <- meta_new
    
    # Removes reductions cause it no longer matches with the new gene set
    obj@reductions <- list()
    
    return(obj)
}

# assume helper functions are already sourced:
# simulate_pseudobulk(), simulate_FC_list(), adjust_seurat_features()

# reuse mapping_organism built earlier (with pbmc_small adjusted to 300 genes)

# pathway / superpathway.gene.sets / hyperparameters are the same as before
# so we can re‐use that code, only swapping out pseudobulk_lognorm and sample_class.

make_input <- function(pseudobulk, labels) {
    new("superpathway.input",
        superpathway_info = superpathway,
        hyperparameters_info= hyperparameters,
        pseudobulk_lognorm = pseudobulk,
        sample_id = as.character(seq_along(labels)),
        sample_class = as.character(labels),
        base_class = "ctrl",
        target_class = "case")
}
