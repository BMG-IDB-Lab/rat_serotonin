source("R-scripts_paper/utils.R")
set.seed(22)
library(Seurat)
library(cacoa)
library(ggplot2)
library(dplyr)
library(tidyr)
# library(cowplot)
# library(magrittr)
library(clusterProfiler)

# Data ----
e20_full <- readRDS("E20_whole_RNA.rds")
e20 <- readRDS("E20_whole_RNA_woZbtb18_Slc17a7_Slit3.Crym.rds")
DefaultAssay(e20) <- "RNA"

# Cacoa object ----
## Meta data ----
meta <- as.data.frame(e20@meta.data)

sample.groups.df <- meta %>% select(treatment, sample_sex) %>%
  group_by(treatment, sample_sex) %>%
  summarise(count = n())
sample.groups <- factor(sample.groups.df$treatment, levels = c("Ctrl", "5HT"))
names(sample.groups) <- sample.groups.df$sample_sex

cell.groups <- factor(meta[, "cluster_markers"]) 
names(cell.groups) <- rownames(meta)

sample.per.cell <- meta$sample_sex
names(sample.per.cell) <- rownames(meta)

## Create cao object ----
cao <- Cacoa$new(e20,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level="Ctrl",
                 target.level="5HT",
                 graph.name="integrated_nn",
                 embedding=e20@reductions$umap@cell.embeddings,
                 verbose = TRUE)
cao[["embedding"]] <- e20_full@reductions$umap.orig@cell.embeddings[Cells(e20), ]

# Cluster-free analysis ----
## Expressional Programs -----
cao$estimateClusterFreeDE(
  n.top.genes=1000, min.expr.frac=0.01, adjust.pvalues=TRUE, smooth=TRUE,
  verbose=TRUE
)
cao$estimateGenePrograms(method="leiden", z.adj=TRUE, smooth=FALSE)

plots_list <- cao_programs(cao)
# Figure 5H
cowplot::plot_grid(plotlist = plots_list[c(2,1,3,4,5,6)], ncol = 2)
# Supplementary Figure 21A (all programs)
cowplot::plot_grid(plotlist = plots_list[c(2,3,5,7,1,4,6)], ncol = 4)



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
  keys = rownames(e20@assays[[assay]]@data)) # !
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

# Figure 5I
### Top 1 by p.adjust 
df_plot <- df_all_ont %>% 
  filter(Cluster != 5) %>% 
  group_by(Cluster, Ontology) %>% 
  arrange(p.adjust, -g, .by_group = T) %>% 
  slice_head(n = 1)
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


