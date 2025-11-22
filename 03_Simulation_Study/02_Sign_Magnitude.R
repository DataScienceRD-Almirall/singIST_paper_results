# ──────────────────────────────────────────────────────────────────────────────
# 1. Setup
# ──────────────────────────────────────────────────────────────────────────────
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT="TRUE")
if (!requireNamespace("singIST", quietly=TRUE)) {
    remotes::install_github("DataScienceRD-Almirall/singIST")
}

# Human simulated data objects
simulated_human <- simulate_pseudobulk(p = 300, effect_size = -20)
pseudobulk_human <- simulated_human$blocks
sample_id <- as.character(1:10)
sample_class <- as.character(simulated_human$labels)

# Fold Change
logFC <- simulate_FC_list(c("Tcell","Bcell","Mono","NK","DC"),
                          simulated_human$genes,
                          fc_min = 0.8, fc_max = 1.2, prop_signif = 1)

###############################################################################
##                         Initiliazing class objects                        ##
###############################################################################
# Pathway
pathway <- new("pathway",
               standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
               dbsource = "KEGG", collection = "c2", subcollection = "CP")
# Superpathway.gene.sets
superpathway <- new("superpathway.gene.sets", pathway_info = pathway,
                    celltypes = c("Tcell","Bcell","Mono","NK","DC"),
                    gene_sets_celltype = rep(list(simulated_human$genes), 5))
# Hyperparameters
hyperparameters <- new("hyperparameters",
                       quantile_comb_table = 
                           matrix(c(rep(0, 5), rep(1, 5)), ncol = 5, nrow = 2,
                                  byrow = TRUE),
                       outcome_type = "binary",
                       number_PLS = as.integer(2),
                       folds_CV = as.integer(1))
# Superpathway.input
superpathway_input <- new("superpathway.input",
                          superpathway_info = superpathway,
                          hyperparameters_info = hyperparameters,
                          pseudobulk_lognorm = pseudobulk_human,
                          sample_id = sample_id,
                          sample_class = sample_class,
                          base_class = "ctrl",
                          target_class = "case")
# mapping organism
counts <- SeuratObject::pbmc_small
counts <- adjust_seurat_features(counts, simulated_human$genes)
counts$class <- counts$groups
counts$celltype_cluster <- counts$RNA_snn_res.0.8
counts$donor <- rep("donor", 80)
mapping_organism <- new("mapping.organism",
                        organism = "Homo sapiens",
                        target_class = "g1",
                        base_class = "g2",
                        celltype_mapping =  list("Tcell" = 1,"Bcell" = 1,
                                                 "Mono" = 1,"NK" = 1,"DC" = 1),
                        counts = counts
)

# ──────────────────────────────────────────────────────────────────────────────
# 2. Simulate & fit for human effect = +20 and = −20
# ──────────────────────────────────────────────────────────────────────────────
# positive effect
sim_pos <- simulate_pseudobulk(p = 300, effect_size =  10000*2^0.2)
colnames(sim_pos$blocks) <- superpathway@gene_sets_celltype[[1]]
sim_pos$genes <- superpathway@gene_sets_celltype[[1]]
input_pos<- make_input(sim_pos$blocks, sim_pos$labels)
opt_pos <- fitOptimal(input_pos)

# negative effect
sim_neg <- simulate_pseudobulk(p = 300, effect_size = 10000*2^(-0.2))
colnames(sim_neg$blocks) <- superpathway@gene_sets_celltype[[1]]
sim_neg$genes <- superpathway@gene_sets_celltype[[1]]
input_neg<- make_input(sim_neg$blocks, sim_neg$labels)
opt_neg <- fitOptimal(input_neg)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Loop over fold‐change intervals
# ──────────────────────────────────────────────────────────────────────────────
fc_width <- 0.1
fc_mins <- seq(-1, 1 - fc_width, by = fc_width)
midpoints <- fc_mins + fc_width/2
n_steps <- length(fc_mins)

recap_pos <- numeric(n_steps)
recap_neg <- numeric(n_steps)

for (i in seq_len(n_steps)) {
    fc_min <- fc_mins[i]
    fc_max <- fc_min + fc_width
    # simulate a fully‐significant FC_list over [fc_min, fc_max]
    FC_list <- simulate_FC_list(
        cell_types  = sim_pos$blocks %>% rownames() %>% sub("_.*", "", .)
                        %>% unique(),
        genes = sim_pos$genes,
        fc_min  = fc_min,
        fc_max = fc_max,
        prop_signif = 1
    )
    # compute recapitulation for each pre‐fitted model
    rec_p <- singISTrecapitulations(mapping_organism, opt_pos,
                                    model_species = "hsapiens",
                                    FC_list = FC_list)
    rec_n <- singISTrecapitulations(mapping_organism, opt_neg,
                                    model_species = "hsapiens",
                                    FC_list = FC_list)
    # extract the superpathway recapitulation (%)
    recap_pos[i] <- rec_p$superpathway$recapitulation
    recap_neg[i] <- rec_n$superpathway$recapitulation
}

# assemble into a data.frame for plotting
df <- data.frame(
    midpoint = midpoints,
    recapit_pos20 = recap_pos,
    recapit_neg20 = recap_neg
)

# ──────────────────────────────────────────────────────────────────────────────
# 4. Plot
# ──────────────────────────────────────────────────────────────────────────────
ggplot(df, aes(x = midpoint)) +
    geom_line(aes(y = recapit_pos20, color = "Human log2FC = 0.20"), size = 1.2) +
    geom_point(aes(y = recapit_pos20, color = "Human log2FC = 0.20"), shape = 16, size = 2) +
    geom_line(aes(y = recapit_neg20, color = "Human log2FC = -0.20"), size = 1.2) +
    geom_point(aes(y = recapit_neg20, color = "Human log2FC = -0.20"), shape = 17, size = 2) +
    geom_hline(yintercept = c(0), color = "gray50") +
    scale_color_manual(name = "Human effect",
                       values = c("Human log2FC = 0.20" = "#1f78b4", "Human log2FC = -0.20" = "#e31a1c")) +
    labs(title    = "Superpathway Recapitulation vs. Disease model log2 Fold-Change Interval Midpoint",
         subtitle = "Comparing positive and negative human effects",
         x        = "Disease model log2 Fold-Change Midpoint",
         y        = "Superpathway Recapitulation (%)") +
    theme_minimal(base_size = 14) +
    scale_y_continuous(breaks = round(seq(-800, 800, by = 50),1)) + 
    scale_x_continuous(breaks = round(seq(-1, 1, by = 0.1),1)) +
    geom_vline(xintercept = 0.20, color = "black", linetype = "dashed") + 
    geom_vline(xintercept = -0.20, color = "black", linetype = "dashed") +
    geom_hline(yintercept = c(-100, 100), color = "black", linetype = "dashed") + 
    theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
        axis.title       = element_text(face = "bold"),
        legend.position  = "top",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )
