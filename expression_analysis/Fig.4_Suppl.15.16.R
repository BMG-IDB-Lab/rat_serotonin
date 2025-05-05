library(Seurat)
library(ggplot2)
library(dplyr)
library(cacoa)
library(cowplot)
library(magrittr)
library(gatom)
library(org.Rn.eg.db)
library(enrichplot)
library(clusterProfiler)

p28_sc <- readRDS("scP28_whole.rds")
p28_sn <- readRDS("snP28_whole.rds")
data <- readRDS("scP28_olig.rds")

# Whole P28 single-cell dataset ----
## Cacoa: Compositional changes ----
DefaultAssay(p28_sc) <- "RNA"
meta <- as.data.frame(p28_sc@meta.data)
sample.groups.df <- meta %>% select(treatment, sample) %>%
  group_by(treatment, sample) %>%
  summarise(count = n())
sample.groups <- factor(sample.groups.df$treatment, levels = c("Ctrl", "5HT"))
names(sample.groups) <- sample.groups.df$sample
cell.groups <- factor(meta$Cell_type) 
names(cell.groups) <- rownames(meta)
sample.per.cell <- meta$sample
names(sample.per.cell) <- rownames(meta)

cao <- Cacoa$new(p28_sc,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level="Ctrl",
                 target.level="5HT",
                 graph.name="integrated_nn",
                 embedding=p28_sc@reductions$umap@cell.embeddings,
                 verbose = TRUE)

## Cluster-based compositional analysis
# Figure 4A
cao$estimateCellLoadings()
cao$plotCellLoadings(show.pvals=TRUE)

## Cluster-free compositional changes 
# Supplementary Figure 5H: Cacoa
cao$estimateCellDensity(method='graph')
cao$estimateDiffCellDensity(type='wilcox')
plot_grid(
  cao$plotEmbedding(color.by='cell.groups'),
  cao$plotDiffCellDensity(
    legend.position=c(1,0)
  ),
  ncol=2
)



# P28 scRNA-Seq oligodendrocytes ----

DefaultAssay(data) <- "RNA"
meta <- as.data.frame(data@meta.data)
sample.groups.df <- meta %>% select(treatment, sample) %>%
  group_by(treatment, sample) %>%
  summarise(count = n())
sample.groups <- factor(sample.groups.df$treatment, levels = c("Ctrl", "5HT"))
names(sample.groups) <- sample.groups.df$sample
cell.groups <- factor(meta$Cell_type) 
names(cell.groups) <- rownames(meta)
sample.per.cell <- meta$sample
names(sample.per.cell) <- rownames(meta)

cao <- Cacoa$new(data,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level="Ctrl",
                 target.level="5HT",
                 graph.name="integrated_nn",
                 embedding=data@reductions$umap@cell.embeddings,
                 verbose = TRUE)

## Cluster-free compositional changes 
# Figure 4B
cao$estimateCellDensity(method='graph')
cao$estimateDiffCellDensity(type='wilcox')

cao_compositional_orderSize(cao)$p2


# Lipid pathways ----

# obtain lipid pathways list
ids.reactome.relation <- read.table("https://reactome.org/download/current/ReactomePathwaysRelation.txt")
colnames(ids.reactome.relation) <- c("id", "child")
org.db = org.Rn.eg.db::org.Rn.eg.db

id.Metabolism_lipid <- "R-RNO-556833"
colnames(ids.reactome.relation) <- c("from", "to")
g <- igraph::graph_from_data_frame(ids.reactome.relation, directed=TRUE)
ids <- id.Metabolism_lipid
res <- igraph::bfs(g, root=ids, mode="out", unreachable=FALSE, dist=TRUE)
ids.reactome.Metabolism_lipid <- names(which(res$dist > 0))


org.gatom.anno <- makeOrgGatomAnnotation(org.db = org.Rn.eg.db)
pathways <- org.gatom.anno$pathways
pathways_lipid <- pathways[which(gsub(":.*", "", names(pathways)) %in% ids.reactome.Metabolism_lipid)]
pathways_lipid_symbol <- lapply(pathways_lipid, function(x) {
  unname(
    AnnotationDbi::mapIds(
      org.Rn.eg.db::org.Rn.eg.db,
      column = "SYMBOL",
      keytype = "ENTREZID",
      keys = x)
  )
})

pathways <- Reduce(rbind, lapply(names(pathways_lipid_symbol), function (x) {
  data.frame(pathway = x,
             geneID = pathways_lipid_symbol[[x]])
}))

# de + enrichment
# olig de + intersection with list 

Idents(data) <- "treatment"

clusters <- c("mOD", "OPC/COP_1", "OPC/COP_2")
FC <- 0
pct <- 0

for (cluster in clusters) {
  
  cluster_name <- gsub("\\/", "\\.", cluster)
  
  cl_orig <- FindMarkers(object = subset(data, subset = Cell_type == cluster),
                            ident.1 = "5HT",
                            ident.2 = "Ctrl",
                            recorrect_umi = FALSE, 
                            assay = "SCT",
                            test.use = "MAST",
                            only.pos = FALSE,
                            logfc.threshold = 0,
                            min.pct = 0,
                            return.thresh = 1
  ) %>% tibble::rownames_to_column("gene") %>% 
    mutate(pct.diff = pct.1 - pct.2) %>% 
    filter(p_val_adj < 0.05, abs(avg_log2FC) > FC, pct.1 > pct | pct.2 > pct) %>% 
    mutate(cluster = ifelse(avg_log2FC > 0, "5HT", "Ctrl"))
  cl <- split(cl_orig, f = cl_orig$cluster)
  cl <- lapply(cl, function(x) x$gene)

  ck <- tryCatch(compareCluster(geneCluster = cl, fun = enricher,
                                TERM2GENE = pathways,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05), error=NULL)
  
  if (!is.null(ck)) {
    dotplot(ck, showCategory = 5, label_format = 50)
  }
}


# Cholesterol biosynthesis genes (and other lipid-related) ----
# Supplementary Figure 15-16
# Figure 3G

genes_cholesterol <- (pathways %>% filter(pathway == "R-RNO-191273: Cholesterol biosynthesis"))$geneID

# volcano plots

plots_volcano <- list()
plots_anno <- list()
olig_main <- list()

# process sc and sn sep
for (dt in c("sc", "sn")) {
  
  if (dt == "sc") {
    data_i <- p28_sc
  } else {
    data_i <- p28_sn
  }
  
  Idents(data_i) <- "treatment"
  Idents(data) <- "treatment"
  
  data_i$Cell_type_short2 <- gsub("_\\d+", "", data_i$Cell_type)
  data_i$Cell_type_short <- ifelse(data_i$Cell_type_short2 %in% c("mOD", "OPC/COP"), 
                                   "Olig", data_i$Cell_type_short2)
  
  Idents(data_i) <- "treatment"
  clusters <- (table(data_i$Cell_type_short, data_i$treatment) %>% as.data.frame() %>% 
                 tidyr::spread(key = "Var2", value = "Freq") %>% 
                 mutate(less20 = ifelse(`5HT` < 20 | Ctrl < 20, T, F),
                        Var1 = as.character(Var1)) %>% 
                 filter(less20 == FALSE) %>% dplyr::select(Var1) %>% 
                 as.vector())$Var1
  clusters <- c("whole", clusters, "mOD", "OPC/COP")
  if (dt == "sc") {clusters <- c(clusters, "sep_mOD", "sep_OPC/COP_1", "sep_OPC/COP_2")}
  
  plot_list1 <- list()
  plot_list2 <- list()
  for (cluster in clusters) {
    
    cluster_name <- gsub("\\/", "\\.", cluster); print(cluster_name)
    
    # DE 
    FC <- 0; pct <- 0
    
    if (cluster_name == "whole") {
      # Whole dataset
      cl_orig <- FindMarkers(object = data_i,
                             ident.1 = "5HT",
                             ident.2 = "Ctrl",
                             assay = "SCT",
                             test.use = "MAST",
                             only.pos = FALSE,
                             logfc.threshold = 0,
                             min.pct = 0,
                             return.thresh = 1
      ) %>% 
        mutate(pct.diff = pct.1 - pct.2) %>% 
        tibble::rownames_to_column("gene")
    } else if (cluster %in% c("mOD", "OPC/COP")) {
      # separate for mOD and OPC
      cl_orig <- FindMarkers(object = subset(data_i, subset = Cell_type_short2 == cluster),
                             ident.1 = "5HT",
                             ident.2 = "Ctrl",
                             recorrect_umi = FALSE, 
                             assay = "SCT",
                             test.use = "MAST",
                             only.pos = FALSE,
                             logfc.threshold = 0,
                             min.pct = 0,
                             return.thresh = 1
      ) %>% 
        mutate(pct.diff = pct.1 - pct.2) %>% 
        tibble::rownames_to_column("gene"); gc()
    } else if ((cluster %in% c("sep_mOD", "sep_OPC/COP_1", "sep_OPC/COP_2")) & dt == "sc") {
      # for P28 scRNA-Seq also perform DE on reintegrated oligodendrocytes object
      cl_orig <- FindMarkers(object = subset(data, subset = Cell_type == gsub("sep_", "", cluster)),
                             ident.1 = "5HT",
                             ident.2 = "Ctrl",
                             recorrect_umi = FALSE, 
                             assay = "SCT",
                             test.use = "MAST",
                             only.pos = FALSE,
                             logfc.threshold = 0.05,
                             min.pct = 0
                             # return.thresh = 1
      ) %>% 
        mutate(pct.diff = pct.1 - pct.2) %>% 
        tibble::rownames_to_column("gene"); gc()
    } else {
      # perform DE for each cell type (including whole Olig)
      cl_orig <- FindMarkers(object = subset(data_i, subset = Cell_type_short == cluster),
                             ident.1 = "5HT",
                             ident.2 = "Ctrl",
                             recorrect_umi = FALSE, 
                             assay = "SCT",
                             test.use = "MAST",
                             only.pos = FALSE,
                             logfc.threshold = 0,
                             min.pct = 0,
                             return.thresh = 1
      ) %>% 
        mutate(pct.diff = pct.1 - pct.2) %>% 
        tibble::rownames_to_column("gene"); gc()
    } 
    # filter DE
    cl_orig <- cl_orig %>% 
      filter(p_val_adj < 0.05, abs(avg_log2FC) > FC, pct.1 > pct | pct.2 > pct) %>% 
      mutate(cluster = ifelse(avg_log2FC > 0, "5HT", "Ctrl"))
    
    # skip cluster if there are no cholesterol genes
    if (length(intersect(cl_orig$gene, genes_cholesterol)) != 0) {
      # Volcano plots 
      p1 <- volcano_de_genes(markers = cl_orig,
                             p_val_adj_cutoff=0.05,
                             avg_log2FC_cutoff=FC,
                             pct_cutoff=pct, 
                             txt.size = 8,
                             genes_highlight = genes_cholesterol, 
                             up_label = "5HT", down_label = "Ctrl") +
        theme(text = element_text(size = 20)) +
        ggtitle(cluster_name)
      
      if (cluster == "sep_mOD") {
        p <- volcano_de_genes(markers = cl_orig,
                               p_val_adj_cutoff=0.05,
                               avg_log2FC_cutoff=FC,
                               pct_cutoff=pct, 
                               txt.size = 8,
                               genes_highlight = genes_cholesterol, 
                               label.left = T,
                               up_label = "5HT", down_label = "Ctrl") +
          theme(text = element_text(size = 20)) +
          ggtitle(cluster_name)
        olig_main[[cluster_name]] <- p
      }
      
      plot_list1[[cluster_name]] <- p1
    } else {
      print("... no genes ...")
    }
    
    # annotation with other pathways
    cl <- split(cl_orig, f = cl_orig$cluster)
    cl <- lapply(cl, function(x) x$gene)
    ck <- tryCatch(compareCluster(geneCluster = cl, fun = enricher,
                                  TERM2GENE = pathways,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05), error=NULL)
    if (!is.null(ck)) {
      p_anno <- dotplot(ck, showCategory = 5, label_format = 50) +
        ggtitle(cluster) +
        theme(axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20, face = "bold"),
              axis.title.x = element_blank(),
              legend.key.size = unit(ifelse(nrow(ck@compareClusterResult %>% filter(p.adjust < 0.05)) > 1,
                                            0.5, 0.15), "cm"),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 25))
      plot_list2[[cluster_name]] <- p_anno
    }
    
    
  }
  
  plots_volcano[[dt]] <- plot_list1
  
  plots_anno[[dt]] <- plot_list2
  
}

# Supplementary Figure 15A
# Cholesterol biosynthesis related differentially expressed genes in P28 scRNA-seq dataset
plot_grid(plotlist = plots_volcano[["sc"]], ncol = 3)
# For main figure genes on the left for olig:
olig_main[["sep_mOD"]]
# Supplementary Figure 15B
# Enriched lipid pathways in cell types of P28 scRNA-seq dataset
plot_grid(plotlist = plots_anno[["sc"]], ncol = 2)

# Supplementary Figure 16A
# Cholesterol biosynthesis related differentially expressed genes in P28 snRNA-seq dataset
plot_grid(plotlist = plots_volcano[["sn"]], ncol = 3)
# Supplementary Figure 16B
# Enriched lipid pathways in cell types of P28 snRNA-seq dataset
plot_grid(plotlist = plots_anno[["sn"]], ncol = 2)

