source("R-scripts_paper/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)

e20 <- readRDS("E20_whole_RNA.rds")
DefaultAssay(e20) <- "RNA"
p28_sn <- readRDS("snP28_whole.rds")
DefaultAssay(p28_sn) <- "SCT"
p28_sn$Cell_type <- factor(p28_sn$Cell_type, levels = c('Neurons',
                                                        'mOD_1', 'mOD_2','mOD_3',
                                                        'OPC/COP_1', 'OPC/COP_2',
                                                        'Astrocytes',
                                                        'VEnC_Mural',
                                                        'MG_MF', "Lepto",
                                                        "Ependyma", "Tany",
                                                        "PTub"))
p28_sn_n <- readRDS("snP28_neurons.rds") 
DefaultAssay(p28_sn_n) <- "SCT"
p28_sc <- readRDS("scP28_whole.rds")
DefaultAssay(p28_sc) <- "SCT"

colors_e20 <- c("mediumpurple2", "tomato1",
                "#DF9200",  "#006FD2", "navajowhite3",
                "seagreen3", 
                "#019B0B", "#B3B3B3")
colors_p28 <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
                RColorBrewer::brewer.pal(n = 8, "Accent")[-4])
set.seed(123)
colors_p28_n <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
                  RColorBrewer::brewer.pal(n = 8, "Accent")[-4],
                  RColorBrewer::brewer.pal(n = 8, "Set1")[-6], rainbow(18), 
                  sample(rainbow(50), size = 50, replace = F)
)[-c(16,19,22)] %>% unique()


# Figure 2 ----

# Figure 2B
myFeaturePlot(e20, 
              features = c("Slc17a6", "Gad1", "Ptprz1", 
                           "Mki67", "Neurog2", "Ascl1"),
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"),
              ncol = 2,
              pt.size = 0.2,
              title.size = 17,
              title.face = "bold.italic",
              reduction = "umap.orig",
              keep.axis.titles = F,
              keep.scale.bar = F,
              keep.axis.numbers = F,
              keep.axis.ticks = F,
              keep.axis.lines = F,
              order = T)

# Figure 2C, see e20_scvelo.py
DimPlot(e20, group.by = "cell_types", 
        reduction = "umap.orig",
        split.by = "treatment", ncol = 1,
        cols = colors_e20) + NoLegend()

# Figure 2D
res <- "Cell_type"
DimPlot(p28_sn, reduction = "umap",
        group.by = res, cols = colors_p28,
        pt.size = 0.2)

# Figure 2E
DimPlot(p28_sn_n, reduction = "umap", group.by = "res.combined_det",
        cols = colors_p28_n, 
        label = T, label.size = 7, repel = T) + NoLegend()
myFeaturePlot(p28_sn_n, 
              features = c("Slc17a6", "Gad1", "Th", "Chat"),
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"),
              ncol = 2,
              pt.size = 0.2,
              title.size = 17,
              title.face = "bold.italic",
              reduction = "umap",
              keep.axis.titles = F,
              keep.scale.bar = F,
              keep.axis.numbers = F,
              keep.axis.ticks = F,
              keep.axis.lines = F,
              order = T)

# Figure 2F
# used list of neuropeptides combined with differentially expressed marker genes
# genes following the pattern"^LOC|^RGD|Gad1|Gad2|Slc17a6|^mt-|ND|COX|ATP" we excluded
genes <- c('Qrfp', 'Ghrh', 'Gal', 'Nts', 'Chat', 'Adcyap1', 'Tac3', 'Nr4a3', 
           'Penk', 'Crh', 'Pdyn', 'Hcrt', 'Sst', 'Cartpt', 'Oxt', 'Pde1c', 
           'Th', 'Ptprk', 'Pnoc', 'Satb2', 'Avp', 'Pomc', 'Tac1', 'Prlr', 
           'Homer1', 'Gpc5', 'Sulf1', 'Agrp', 'Npy', 'Tcf4', 'Rorb', 'Pmch', 
           'Adarb2', 'Reln', 'Ecel1', 'Grp', 'Vip', 'NMS', 'Rarb', 'Trh', 
           'Tmem163', 'Ache', 'Erbb4', 'Slc17a7', 'Cdh23', 'Gda', 'Syt2', 
           'Cdh4', 'Gabrg1', 'Calcr', 'Fign', 'Nrsn2', 'Il1rapl2', 'Bdnf', 
           'Scd2', 'Npas3', 'Nfib', 'Meis2', 'Scn9a', 'Cck', 'Tafa1')
DotPlot(p28_sn_n,
        features = genes,
        dot.min = 0.1,
        group.by = "res.combined_det",
        # cluster.idents = T,
        cols = c("blue", "red")) +
  coord_flip() + RotatedAxis() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face = "italic", size = 22),
        axis.text.x = element_text(size = 18, angle = 45),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(20,10,10,10)) +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

DotPlot(p28_sn_n,
        features = c("Slc17a6", "Gad1"),
        group.by = "res.combined_det",
        dot.min = 0.1,
        # cluster.idents = T,
        cols = c("blue", "red")) +
  coord_flip() + RotatedAxis() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(face = "italic", size = 22),
    axis.text.x = element_text(size = 18),
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    plot.margin = margin(20, 10, 10, 10),
    axis.ticks.x = element_blank()
  )

# Figure 2G
genes <- c("Avp",  "Agrp", "Pomc", 
           "Crh", "Trh", "Sst",
           "Gal", "Pmch", "Vip")
myFeaturePlot(p28_sn_n, features = genes,
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"),
              ncol = 3,
              pt.size = 0.8,
              title.size = 17,
              reduction = "umap",
              keep.axis.titles = F,
              keep.scale.bar = F,
              keep.axis.numbers = F,
              keep.axis.ticks = F,
              keep.axis.lines = F,
              title.face = "bold.italic",
              min.cutoff = "q25",
              order = T)



# Supplementary Figure 2 ----

# Supplementary Figure 2B-D
data_list = list(E20 = e20,
                 P28_sn = p28_sn,
                 P28_sc = p28_sc)

for (name in names(data_list)) {
  
  data <- data_list[[name]] 

  # RNA add log10
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  # SCT add log10
  data[['nCount_SCT_log10']] <- log10(data[['nCount_SCT']] + 1)
  data[['nFeature_SCT_log10']] <- log10(data[['nFeature_SCT']] + 1)
  
  ## data
  
  p <- VlnPlot(data, group.by = "sample", 
               features = c("nCount_RNA", "nFeature_RNA", "nCount_RNA_log10", "nFeature_RNA_log10",
                            "nCount_SCT", "nFeature_SCT", "nCount_SCT_log10", "nFeature_SCT_log10"),
               combine = F)
  cowplot::plot_grid(
    plotlist = lapply(p, function(x) x + theme(axis.title = element_blank()) + NoLegend()),
    ncol = 4)

  rm(data); gc()
}; rm(data_list); gc()



# Supplementary Figure 3 ----

# Supplementary Figure 3A
DimPlot(e20, group.by = "cell_types", 
        reduction = "umap.orig",
        cols = colors_e20) 

# Supplementary Figure 3B
library(tricycle)
data.m <- GetAssayData(object = e20, assay = "RNA", slot = "data")
projection.m <- tricycle:::.project_cycle_space(data.m, 
                                                gname.type = "SYMBOL", 
                                                species = "mouse")
e20$tricyclePosition <- tricycle:::.getTheta(projection.m,
                                             center.pc1 = 0,
                                             center.pc2 = 0)

e20[["tricycleEmbedding"]] <- CreateDimReducObject(
  embeddings = projection.m,
  key = "tricycleEmbedding_",
  assay = "RNA"
)
p <- FeaturePlot(e20, features = "tricyclePosition", pt.size = 2, reduction = "umap.orig") +
  scale_color_gradientn(limits = range(0, 2 * pi), 
                        breaks = seq(from = 0, to = 2 * pi, length.out = 500), 
                        colors = alpha(c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", 
                                         "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"), .6), 
                        guide = "none")
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
cowplot::plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))


# Supplementary Figure 3C - see e20_scvelo.py


# Supplementary Figure 3D-3H
myFeaturePlot(e20, features = c("Hes1", "Ascl1", "Dll1", "Dll3", "Rbpj", "Rbfox3",
                                "Dlx1", "Dlx5", "Neurog1", "Neurog2", "Neurod1", "Neurod2"),
              ncol = 6, order = T, pt.size = 0.5, 
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"),
              # min.cutoff = "q10",
              assay = "RNA", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = T, 
              title.size = 15,
              title.face = "bold.italic")

myFeaturePlot(e20, features = c("Pdyn", "Tac1", "Pomc", "Pnoc", "Nts", "Gal",
                                "Agrp", "Npy", "Sst", "Cartpt", "Pmch", "Rln1",
                                "Penk", "Trh", "Qrfp", "Hcrt", "Npvf", "Ghrh",
                                "Grp", "Crh", "Bdnf", "Adcyap1", "Cck", "Vip"),
              ncol = 6, order = T, pt.size = 0.5, 
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"),
              assay = "RNA", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = T, 
              title.size = 15,
              title.face = "bold.italic")

myFeaturePlot(e20, features = c("Oxt", "Avp", "Sox2", "Sox9", "Sox10", "Slc6a3"),
              ncol = 6, order = T, pt.size = 0.5,
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"), 
              assay = "RNA", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = T, 
              title.size = 15,
              title.face = "bold.italic")


# Supplementary Figure 4 ----

# Supplementary Figure 4B
genes <- c("Maoa", "Maob", "Ddc", "Slc6a4")
myFeaturePlot(e20, 
              features = genes,
              ncol = 2, order = T, pt.size = 0.8, 
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"), 
              assay = "RNA", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = T, 
              title.size = 18,
              reduction = "umap.orig",
              title.face = "bold.italic", 
              legend.key.size = unit(0.5, "cm"))

# Supplementary Figure 4C-D
genes <- c(paste0("Htr",
                  c("1a", "1b", "1d", "1f", "2a", "2b",
                    "2c", "3a", "3b", "4", "5a", "5b", 
                    "6", "7")),
           "Tph1", "Tph2")
myFeaturePlot(e20, 
              features = genes,
              ncol = 4, order = T, pt.size = 0.8, 
              custom.cols = c("darkblue", "cyan", "yellow",  "red", "red", "darkred", "coral4"), 
              assay = "RNA", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = T, 
              title.size = 18,
              reduction = "umap.orig",
              title.face = "bold.italic", 
              legend.key.size = unit(0.5, "cm"))


# Supplementary Figure 6 ----

# Supplementary Figure 6A
id_new <- p28_sn_n$res.combined_det_markers
cell_names <- names(id_new); id_new <- as.character(id_new)
id_new[grep("51_Dop_Chol", id_new)] <- "51_Mixed_Eef1b2.Cd59b"
names(id_new) <- cell_names
id_new <- factor(id_new, levels = stringr::str_sort(unique(as.character(id_new)), numeric = T))
p28_sn_n <- AddMetaData(p28_sn_n, metadata = id_new, col.name = "res.combined_det_markers_mix")

DimPlot(p28_sn_n, reduction = "umap", group.by = "res.combined_det_markers_mix",
        cols = colors_p28_n, label = F)
DimPlot(p28_sn_n, reduction = "umap", group.by = "res.combined_det",
        cols = colors_p28_n, label = T, repel = T) + NoLegend()


# Supplementary Figure 6B
# used list of neuropeptides combined with differentially expressed marker genes
# genes following the pattern"^LOC|^RGD|Gad1|Gad2|Slc17a6|^mt-|ND|COX|ATP" we excluded
genes <- c('Qrfp', 'Ghrh', 'Gal', 'Nts', 'Chat', 'Adcyap1', 'Tac3', 'Nr4a3', 
           'Penk', 'Crh', 'Pdyn', 'Hcrt', 'Sst', 'Cartpt', 'Oxt', 'Pde1c', 
           'Th', 'Ptprk', 'Pnoc', 'Satb2', 'Avp', 'Pomc', 'Tac1', 'Prlr', 
           'Homer1', 'Gpc5', 'Sulf1', 'Agrp', 'Npy', 'Tcf4', 'Rorb', 'Pmch', 
           'Adarb2', 'Reln', 'Ecel1', 'Grp', 'Vip', 'NMS', 'Rarb', 'Trh', 
           'Tmem163', 'Ache', 'Erbb4', 'Slc17a7', 'Cdh23', 'Gda', 'Syt2', 
           'Cdh4', 'Gabrg1', 'Calcr', 'Fign', 'Nrsn2', 'Il1rapl2', 'Bdnf', 
           'Scd2', 'Npas3', 'Nfib', 'Meis2', 'Scn9a', 'Cck', 'Tafa1')
myFeaturePlot(p28_sn_n, 
              features = sort(genes),
              ncol = 7, order = T, pt.size = 0.8, 
              assay = "SCT", slot = "data", 
              keep.axis.titles = F, keep.axis.numbers = F, 
              keep.axis.ticks = F, keep.axis.lines = F, 
              keep.scale.bar = F, 
              title.size = 18,
              title.face = "bold.italic")
