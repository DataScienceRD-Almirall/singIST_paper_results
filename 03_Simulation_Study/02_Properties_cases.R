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

###############################################################################
##                         Proportion of DEGs case                           ##
###############################################################################

# Prepare grid of proportions
props <- seq(0, 1, by = 0.01)
n_steps <- length(props)
recap_vals <- numeric(n_steps)

# Loop to compute recapitulation at each prop_signif
for (i in seq_along(props)) {
    FC_list <- simulate_FC_list(
        cell_types  = c("Tcell","Bcell","Mono","NK","DC"),
        genes       = sim_neg$genes,
        fc_min      = 0.19,
        fc_max      = 0.21,
        prop_signif = props[i]
    )
    rec      <- singISTrecapitulations(mapping_organism, opt_pos,
                                       model_species = "hsapiens",
                                       FC_list = FC_list)
    recap_vals[i] <- rec$superpathway$recapitulation
}

df2 <- data.frame(
    prop_signif    = props,
    recapitulation = recap_vals
)

# Plot
ggplot(df2, aes(x = prop_signif, y = recapitulation)) +
    geom_line(color = "#e31a1c", size = 1.2) +
    geom_point(color = "#e31a1c", shape = 17, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title    = "Effect of Proportion of Significant Genes on Superpathway Recapitulation",
         subtitle = "Human effect log2FC = 0.20; Fixed disease model log2FC interval Unif[0.19, 0.21]",
         x        = "Proportion of Genes with Adjusted p ≤ 0.05",
         y        = "Superpathway Recapitulation (%)") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    ) + 
    scale_y_continuous(n.breaks = 15)

###############################################################################
##                         Cell Type drop out case                           ##
###############################################################################

# Define cell types and fixed FC_list for [0.1, 0.5]
cell_types     <- c("Tcell","Bcell","Mono","NK","DC")
FC_list_fixed  <- simulate_FC_list(
    cell_types   = cell_types,
    genes        = sim_pos$genes,
    fc_min       = 0.19,
    fc_max       = 0.21,
    prop_signif  = 1.0
)

# Prepare to collect recapitulation values
n_cells    <- 1:length(cell_types)
recap_vals <- numeric(length(n_cells))

# Loop over number of mapped cell types
for (i in seq_along(n_cells)) {
    k <- n_cells[i]
    # Build new celltype_mapping: first k mapped, rest unmapped
    ct_map <- setNames(lapply(seq_along(cell_types), function(j) {
        if (j <= k) seq_len(1) else integer(0)
    }), cell_types)
    
    # Update mapping.organism object
    mapping_k <- mapping_organism
    mapping_k@celltype_mapping <- ct_map
    
    # Compute recapitulation
    rec <- singISTrecapitulations(
        mapping_k,
        opt_pos,
        model_species = "hsapiens",
        FC_list       = FC_list_fixed
    )
    recap_vals[i] <- rec$superpathway$recapitulation
}

# Assemble data.frame
df3 <- data.frame(
    n_mapped       = n_cells,
    recapitulation = recap_vals
)

# Create the plot
ggplot(df3, aes(x = n_mapped, y = recapitulation)) +
    geom_line(color = "#2c7fb8", size = 1.2) +
    geom_point(color = "#2c7fb8", size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
        title    = "Superpathway Recapitulation vs. Number of Mapped Cell Types",
        subtitle = "Fixed log2FC interval [0.19, 0.21]; Human effect log2FC = 0.20",
        x        = "Number of Cell Types Mapped",
        y        = "Superpathway Recapitulation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )

###############################################################################
##                         Orthology case                                    ##
###############################################################################

# Plot Recapitulation vs. “Observed” One‐to‐One Orthology (i.e. fraction of genes retained)
#
# We start from a full FC_list (100% orthologous) and randomly drop
# a given fraction of genes from every block’s data.frame. Then we
# compute recapitulation as that fraction goes from 0→1.

# 1. Prepare fixed FC_list for [0.1,0.5] with 100% of genes present
cell_types    <- c("Tcell","Bcell","Mono","NK","DC")
FC_list_full  <- simulate_FC_list(
    cell_types  = cell_types,
    genes       = sim_pos$genes,
    fc_min      = 0.19,
    fc_max      = 0.21,
    prop_signif = 1.0
)

# 2. Grid of “observed orthology” proportions and storage
props <- seq(0, 1, by = 0.05)
recap_obs <- numeric(length(props))

# 3. Loop: drop fraction of genes from each block
for (i in seq_along(props)) {
    prop_keep <- props[i]
    FC_list_i <- lapply(FC_list_full, function(df) {
        n_genes <- nrow(df)
        n_keep  <- floor(prop_keep * n_genes)
        keep_idx <- sample(n_genes, n_keep)
        df[keep_idx, , drop = FALSE]
    })
    
    # compute recapitulation with the reduced FC_list
    rec <- singISTrecapitulations(
        mapping_organism,
        opt_pos,
        model_species = "hsapiens",
        FC_list       = FC_list_i
    )
    recap_obs[i] <- rec$superpathway$recapitulation
}

# 4. Assemble data and plot
df_obs <- data.frame(
    pct_orthologous = props * 100,
    recapitulation  = recap_obs
)

ggplot(df_obs, aes(x = pct_orthologous, y = recapitulation)) +
    geom_line(color = "#4daf4a", size = 1.2) +
    geom_point(color = "#4daf4a", size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
        title    = "Superpathway Recapitulation vs. Observed One-to-One Orthology",
        subtitle = "Fixed FC interval [0.19,0.21]; Human effect log2FC = 0.20",
        x        = "Percent of Genes Retained (Observed one-to-one Orthology)",
        y        = "Superpathway Recapitulation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )

###############################################################################
##                         +100%/-100% recapitulation case                   ##
###############################################################################

# toy recapitulation function
fc_width <- 0.1
fc_mins <- seq(0, 10 - fc_width, by = fc_width)
midpoints <- fc_mins + fc_width/2
n_steps <- length(fc_mins)
recap <- numeric(n_steps)
recap_neg <- numeric(n_steps)
for (i in seq_len(n_steps)) {
    fc_min <- fc_mins[i]
    fc_max <- fc_min + fc_width
    # Human
    # positive effect
    sim <- simulate_pseudobulk(p = 300, effect_size = 10000*2^midpoints[i], cor = TRUE)
    colnames(sim$blocks) <- superpathway@gene_sets_celltype[[1]]
    sim$genes <- superpathway@gene_sets_celltype[[1]]
    inputs<- make_input(sim$blocks, sim$labels)
    opt_pos <- fitOptimal(inputs)
    #negative effect
    sim <- simulate_pseudobulk(p = 300, effect_size = 10000*2^(-midpoints[i]), cor = TRUE)
    colnames(sim$blocks) <- superpathway@gene_sets_celltype[[1]]
    sim$genes <- superpathway@gene_sets_celltype[[1]]
    inputs<- make_input(sim$blocks, sim$labels)
    opt_neg <- fitOptimal(inputs)
    # simulate a fully‐significant FC_list over [fc_min, fc_max]
    FC_list <- simulate_FC_list(
        cell_types  = sim$blocks %>% rownames() %>% sub("_.*", "", .)
        %>% unique(),
        genes = sim$genes,
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
    recap[i] <- rec_p$superpathway$recapitulation
    recap_neg[i] <- rec_n$superpathway$recapitulation
}

df <- data.frame(midpoints, recap)

# Fit lm models
mod_neg <- lm(recap_neg ~ midpoints, data = df)
mod_pos <- lm(recap ~ midpoints, data = df)
# get coefficient + 95% CI for each
sum_neg <- tidy(mod_neg, conf.int = TRUE)
sum_pos <- tidy(mod_pos, conf.int = TRUE)
# building label for each
lbl_neg <- with(sum_neg,
                sprintf(
                    "-100 model:\nIntercept = %.2f (%.2f to %.2f)\nSlope  = %.3f (%.3f to %.3f)",
                    estimate[1], conf.low[1], conf.high[1],
                    estimate[2], conf.low[2], conf.high[2]
                ))
lbl_pos <- with(sum_pos,
                sprintf(
                    "100 model:\nIntercept = %.2f (%.2f to %.2f)\nSlope  = %.3f (%.3f to %.3f)",
                    estimate[1], conf.low[1], conf.high[1],
                    estimate[2], conf.low[2], conf.high[2]
                ))
ann <- data.frame(
    midpoints = rep(min(midpoints) + 0.05*diff(range(midpoints)), 2),
    y = c(max(recap)+40, max(recap_neg)+40),
    label = c(lbl_pos, lbl_neg),
    color = c("steelblue","firebrick")
)

ggplot(df, aes(x=midpoints, y=recap)) +
    geom_point(aes(y=recap, color = "Human log2FC = Disease model log2FC")) +
    geom_point(aes(y=recap_neg, color = "Human log2FC = - Disease model log2FC")) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    labs(x = "Disease model log2FC", y = "Recapitulation (%)") +
    ylim(c(-200, 200)) +
    scale_color_manual(name = "Human effect",
                       values = c("Human log2FC = Disease model log2FC" = "steelblue", "Human log2FC = - Disease model log2FC" = "firebrick")) + 
    labs(title    = "Superpathway recapitulation when human and disease model log2FC perfectly align and oppose",
         subtitle = "",
         x        = "log2FC",
         y        = "Superpathway Recapitulation (%)") +
    theme_minimal(base_size = 14) +
    geom_text(data = ann,
              aes(x = midpoints, y = y, label = label, color = color),
              hjust = 0, size = 4.5,
              show.legend = FALSE) +
    theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
        axis.title       = element_text(face = "bold"),
        legend.position  = "top",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank()
    )


