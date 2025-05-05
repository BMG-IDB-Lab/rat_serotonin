source("R-scripts_paper/utils.R")
set.seed(22)
library(Seurat)
library(cacoa)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(magrittr)
library(clusterProfiler)

# Data ----
p28_sn_n <- readRDS("snP28_neurons.rds") 
p28_sn_n$sample_sex <- paste0(p28_sn_n$sample, "_", p28_sn_n$sex_predicted)
DefaultAssay(p28_sn_n) <- "RNA"

legend_pos <- c(0,1)

# Neuronal changes in P28 

# Cacoa object ----
## Meta data ----
meta <- as.data.frame(p28_sn_n@meta.data)

sample.groups.df <- meta %>% select(treatment, sample_sex) %>%
  group_by(treatment, sample_sex) %>%
  summarise(count = n())
sample.groups <- factor(sample.groups.df$treatment, levels = c("Ctrl", "5HT"))
names(sample.groups) <- sample.groups.df$sample_sex

cell.groups <- factor(meta[, "res.combined_det"]) 
names(cell.groups) <- rownames(meta)

sample.per.cell <- meta$sample_sex
names(sample.per.cell) <- rownames(meta)

## Create cao object ----
cao <- Cacoa$new(p28_sn_n,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level="Ctrl",
                 target.level="5HT",
                 graph.name="integrated_nn",
                 embedding=p28_sn_n@reductions$umap@cell.embeddings,
                 verbose = TRUE)


# Cluster-free analysis ----
## Compositional changes ----
cao$estimateCellDensity(method='graph')
cao$estimateDiffCellDensity(type='wilcox')

plot_grid(
  cao$plotEmbedding(color.by='cell.groups'),
  cao$plotDiffCellDensity(
    legend.position=legend_pos
  ),
  ncol=2
)

# Supplementary Figure 11A (full dataset)
cao$plotDiffCellDensity(
  size = 2, alpha = 0.5,
  legend.position=NULL
)
cao_compositional_rmCluster(cao, data = p28_sn_n, cluster = NULL)$plot

# Figure 3A (left, excluding DG cluster)
Idents(p28_sn_n) <- "res.combined_det_markers"
cao_compositional_rmCluster(cao, data = p28_sn_n, cluster = "17_Glut_Slc17a7.Nfib")$plot_zoom


## Expressional changes ----
cao$estimateClusterFreeExpressionShifts(
  gene.selection="expression", min.n.between=3, min.n.within=3
)
# Supplementary Figure 11B (both shifts and z-score)
cao$plotClusterFreeExpressionShifts(legend.position=legend_pos, font.size=3)

# Figure 3A (right, z-score with labeled top clusters)
Idents(p28_sn_n) <- "res.combined_det_markers"
plots <- cao_expression_top(cao, data = p28_sn_n,
                            top = 10,
                            cluster_rm = "51_Dop_Chol_Eef1b2.Cd59b")
plots$plot_zscore
plots$plot_barplot


### Expressional Programs -----
cao$estimateClusterFreeDE(
  n.top.genes=1000, min.expr.frac=0.01, adjust.pvalues=TRUE, smooth=TRUE,
  verbose=TRUE
)
cao$estimateGenePrograms(method="leiden", z.adj=TRUE, smooth=FALSE)
cao$plotGeneProgramScores(
  legend.position=legend_pos, plot.na=FALSE, size = 0.3,
  adj.list=theme(legend.key.width=unit(8, "pt"), legend.key.height=unit(12, "pt"))
)

plots_list <- cao_programs(cao)
# Figure 3B (left)
cowplot::plot_grid(plotlist = plots_list[c(1,2,3,4,6,7)], ncol = 2)
# Supplementary Figure 11C (all programs)
cowplot::plot_grid(plotlist = plots_list, ncol = 4)

# Supplementary Figure 11D
plots_list <- cao_programs_barplot(cao, p28_sn_n, anno = "res.combined_det_markers")
cowplot::plot_grid(plotlist = plots_list, ncol = 4)


### Annotation of expressional programs ----

pr_s <- cao[["test.results"]][["gene.programs"]]
pr <- pr_s[["genes.per.clust"]]
pr.id <- (sapply(pr, length) > 10) %>% which()
pr <- pr[pr.id]
rm(pr_s); gdata::mv("pr", "cl"); gc()
names(cl) <- as.character(as.numeric(names(cl)) + 1)

assay <- "RNA"
## universe 
genes_universe <- AnnotationDbi::mapIds(
  org.Rn.eg.db::org.Rn.eg.db,
  column = "ENTREZID",
  keytype = "SYMBOL",
  keys = rownames(p28_sn_n@assays[[assay]]@data)) # !
genes_universe <- genes_universe[!is.na(genes_universe)]


## symbol to entrez 
cl <- lapply(cl, function(x) {
  AnnotationDbi::mapIds(
    org.Rn.eg.db::org.Rn.eg.db,
    column = "ENTREZID",
    keytype = "SYMBOL",
    keys = x)
})

# start anno
df_all_ont <- data.frame()
for (go_ont in c("BP", "MF")) {
  
  df_all <- data.frame()
  for (cluster in names(cl)) {
    
    ck <- tryCatch(enrichGO(gene          = cl[[cluster]],
                            universe      = genes_universe,
                            OrgDb         = org.Rn.eg.db::org.Rn.eg.db,
                            ont           = go_ont, 
                            minGSSize = 5,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE), error=NULL)
    if (!is.null(ck)) {
      ck@result <- ck@result[ck@result$p.adjust<=0.05,]
      if (nrow(ck@result[ck@result$p.adjust<=0.05,]) != 0) {
        ck <- simplify(ck, cutoff=0.9, by="p.adjust", select_fun=min)
        
        terms <- (ck@result %>% as.data.frame() %>% 
                    mutate(Cluster = cluster) %>% 
                    filter(p.adjust <= 0.05) %>%
                    slice_min(order_by = p.adjust, n = 15))[,"Description",drop=T] %>% unique()
        terms_df <- ck@result %>% as.data.frame() %>% 
          mutate(Cluster = cluster) %>% 
          filter(Description %in% terms, p.adjust <= 0.05) %>% 
          arrange(p.adjust) %>% 
          dplyr::select(Cluster, ID, Description, GeneRatio, p.adjust, geneID) %>% 
          mutate(Ontology = go_ont)
        df_all <- rbind(df_all, terms_df)
      }
    }
  }

  df_all <- df_all %>% 
    tidyr::separate(col = GeneRatio, into = c("g", "g_pr"), sep = "/") %>% 
    mutate(g = as.numeric(g),
           g_pr = as.numeric(g_pr),
           ratio = g / g_pr)
  
  df_all_ont <- rbind(df_all_ont, df_all)
  
}


# Figure 3B (right)
### Top 2 by p.adjust 
df_plot <- df_all_ont %>% 
  filter(Cluster != 5) %>% 
  group_by(Cluster, Ontology) %>% 
  arrange(p.adjust, -g, .by_group = T) %>% 
  slice_head(n = 2)
df_plot$Description <- factor(df_plot$Description, levels = rev(unique(df_plot$Description)))

ggplot(df_plot, aes(y = Description, x = Cluster, color = p.adjust)) +
  geom_point(aes(size = ratio)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15, color = "black", 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black")) + 
  theme(strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(colour = 'black', face = "bold", size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  scale_y_discrete(labels = scales::label_wrap(50)) +
  viridis::scale_color_viridis(option = "inferno", begin = 0.1, end = 0.8, direction = -1) +
  theme(plot.margin = margin(10, 10, 10, 10))


# Supplementary Figure 11D (annotation of the 5th program)
### Top 2 by p.adjust 
df_plot <- df_all_ont %>% 
  filter(Cluster == 5) %>% 
  group_by(Cluster, Ontology) %>% 
  arrange(p.adjust, -g, .by_group = T) %>% 
  slice_head(n = 2)
df_plot$Description <- factor(df_plot$Description, levels = rev(unique(df_plot$Description)))

ggplot(df_plot, aes(y = Description, x = Cluster, color = p.adjust)) +
  geom_point(aes(size = ratio)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15, color = "black", 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black")) + 
  theme(strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(colour = 'black', face = "bold", size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  scale_y_discrete(labels = scales::label_wrap(20)) +
  viridis::scale_color_viridis(option = "inferno", begin = 0.1, end = 0.8, direction = -1) +
  theme(plot.margin = margin(10, 80, 10, 20),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left")
