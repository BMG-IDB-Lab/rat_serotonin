source("R-scripts_paper/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(cacoa)
library(tidyr)
library(cowplot)

p28_sc <- readRDS("scP28_whole.rds")
DefaultAssay(p28_sc) <- "SCT"

# colors <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
#             RColorBrewer::brewer.pal(n = 8, "Accent")[-4],
#             RColorBrewer::brewer.pal(n = 8, "Set1")[-6], rainbow(31))[-c(16,19,22,23)]
# colors <- colors[-17]
# colors <- replace(colors, c(8, 9), colors[c(9, 8)])
# colors <- replace(colors, c(7, 14), colors[c(14, 7)])
# colors <- colors[c(1:6,14,9,8,10,11,12,7,15,16,19)]

colors <- c('#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', 
            '#B3B3B3', '#7FC97F', '#BEAED4', '#FDC086', '#386CB0', '#F0027F', 
            '#666666', '#E41A1C', '#4DAF4A', '#FF6300')

# General overview plots of P28 single-cell RNA-Seq dataset
# Supplementary Figure 5A-D
# A
DimPlot(p28_sc, group.by = "Cell_type", 
        reduction = "umap",
        split.by = "sample", ncol = 5,
        cols = colors)
# B
DimPlot(p28_sc, group.by = "Cell_type", 
        reduction = "umap",
        split.by = "sex_predicted", ncol = 2,
        cols = colors)
# C
DimPlot(p28_sc, group.by = "Cell_type", 
        reduction = "umap",
        split.by = "seq", ncol = 2,
        cols = colors)
# D
p28_sc$sex_tr <- paste0(p28_sc$sex_predicted, "_", p28_sc$treatment)
DimPlot(p28_sc, group.by = "Cell_type", 
        reduction = "umap",
        split.by = "sex_tr", ncol = 2,
        cols = colors)

# Main marker genes used for the cell annotation
# Supplementary Figure 5E
FeaturePlot(p28_sc, 
            c("Tmem130", "Plp1", "Vcan", "Agt", "Dnah3", "Crym", 
              "Rgs5", "Flt1", "P2ry12", "Cd163", "Dcn", "Cga"), 
            ncol = 4, pt.size = 0.1, order = T)

# Top differentially expressed of each cluster
# Supplementary Figure 5F
Idents(p28_sc) <- "Cell_type"
p28_sc <- PrepSCTFindMarkers(p28_sc)
markers <- FindAllMarkers(object = p28_sc,
                          assay='SCT',
                          only.pos = TRUE,
                          logfc.threshold = 0.25, min.pct = 0.1,
                          test.use = 'MAST')
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers$cluster <- factor(markers$cluster, levels=c('Neurons', 
                                                    'mOD',
                                                    'OPC/COP_1', 'OPC/COP_2',
                                                    'Astrocytes_1', 'Astrocytes_2', 
                                                    'Ependyma', 'Tanycytes',
                                                    'VEnC_1', 'VEnC_2', 'Mural_cells',
                                                    'Microglia', 'Macrophages', 
                                                    'Fibro_1', 'Fibro_2', "PTub")) 
# plot top 3 genes
cell_markers <- markers %>%
  filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T) %>% 
  slice_head(n = 3)
DoHeatmap(p28_sc, features = cell_markers$gene, 
          size = 3.5, group.colors = colors, angle = 60, 
          group.bar.height = 0.01, lines.width = 100,
          group.by = "Cell_type") +
  viridis::scale_fill_viridis(na.value = "white") +
  guides(color = "none") +
  theme(axis.text.y.left = element_text(size = 15))


# Differences in cellular composition: percanteges of cells from each sample
# Supplementary Figure 5G
prop <- proportionsPlot(obj = p28_sc, 
                        condition = "treatment", condition1 = "5HT", condition2 = "Ctrl",
                        samples = "sample", rev_group_color = F, 
                        group_color_4var = T, ngroups1 = c(3,6), ngroups2 = c(3,6),
                        sample_levels = (p28_sc@meta.data[, c("sample", "treatment", "seq")] %>% 
                                           as.data.frame() %>% distinct() %>% arrange(desc(treatment), desc(seq())))$sample,
                        cell_type = "Cell_type")
prop$cond_sample + theme(
  axis.text.x = element_text(size = 20),
  title = element_text(size = 25))


# Supplementary Figure 5H (see script for Fig.4)