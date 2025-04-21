"""
This script reproduce analysis of connections between proximal and distal tCREs 
for hypothalamic samples of Rattus Norvegicus at day 28 of pastnatal development (P28 stage).
Using Cicero we reconstructed connections between tCREs for control and 5HTP-stimulated samples separately 
and compared the reconstructed conections between conditions for NFI-family transcription factors.
Using ChromVar we analyze motif activity of NFI TFs.
Download .zip archive from 10.5281/zenodo.15237638, unpack it and make folder <scafe_cicero_chromvar> as working directory.
"""

setwd("~/scafe_cicero_chromvar/")

set.seed(1234)
library(monocle3)
#library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cicero)
library(Signac)
library(EnsDb.Rnorvegicus.v79)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Rnorvegicus.UCSC.rn7)
library(chromVAR)

# 1 Upload P28 integrated SeuratObject containig peak x cell count matrix after running SCAFE --------------
p28.scafe = readRDS('./p28.scafe.Rds')

# split seurat object by stim 
DefaultAssay(p28.scafe) = 'RNA'
unique(p28.scafe@meta.data$stim)
serotonin.scafe.list = list()
for(st in unique(p28.scafe@meta.data$stim)){
  print(st)
  v = subset(p28.scafe, subset = stim == st)
  print(dim(v))
  print(unique(v@meta.data$stim))
  serotonin.scafe.list = base::append(serotonin.scafe.list, v)
}
names(serotonin.scafe.list) = unique(p28.scafe@meta.data$stim)

# 2. Convert to cell_data_set object ------------------------
cds.list <- lapply(serotonin.scafe.list, function(x) as.cell_data_set(x))
# since it misses the gene_short_name column, let's add it
for(i in 1:length(cds.list)){
  fData(cds.list[[i]])$gene_short_name <- rownames(cds.list[[i]])
}
lapply(cds.list, function(x) fData(x))

# Cluster cells (using clustering info from seurat's UMAP)
# let's use the clustering information have
recreate.partitions = lapply(cds.list, function(x) 
  recreate.partition = c(rep(1,length(x@colData@rownames))))
for(i in 1:length(cds.list)){
  names(recreate.partitions[[i]]) = cds.list[[i]]@colData@rownames
  recreate.partitions[[i]] = as.factor(recreate.partitions[[i]])}
for(i in 1:length(cds.list)){
  cds.list[[i]]@clusters$UMAP$partitions <- recreate.partitions[[i]]}
# Assign the cluster info 
lapply(serotonin.scafe.list, function(x) x@active.ident)
list_clusters = lapply(serotonin.scafe.list, function(x) x@active.ident)
for(i in 1:length(cds.list)){
  cds.list[[i]]@clusters$UMAP$clusters <- list_clusters[[i]]}
# Assign UMAP coordinate - cell embeddings
for(i in 1:length(cds.list)){
  cds.list[[i]]@int_colData@listData$reducedDims$UMAP  = serotonin.scafe.list[[i]]@reductions$umap@cell.embeddings}

# 3. Make cicero cds -----------------------------------------
umap_coords = lapply(cds.list, function(x) reducedDims(x)$UMAP)
cicero_cds = list()
for(i in 1:length(cds.list)){
  cicero_cds[[i]] =  make_cicero_cds(cds.list[[i]], reduced_coordinates = umap_coords[[i]])
  names(cicero_cds)[i] = names(cds.list)[i]}

# 4. Run cicero -----------------------------------------
mratbn7.2 = read.table('./rat.bn7.2.genome.txt')

conns <- lapply(cicero_cds, function(x) run_cicero(x, mratbn7.2, sample_num = 100))
conns = lapply(conns, na.omit)
conns = lapply(conns, function(x) dplyr::filter(x, coaccess > 0.05))

# 5. Annotate peaks -----------------------------------------
all.tCREs = read.csv("./p28_all.tCREs.csv")
all.tCREs$Peak1 = all.tCREs$peak
# Заменить все значения NA на "Unknown" в колонке "SYMBOL"
all.tCREs <- all.tCREs %>% mutate(SYMBOL = ifelse(is.na(SYMBOL), "Unknown", SYMBOL))
# add annotation to Peak in conns
# for conns
for(i in 1:length(conns)){
  conns[[i]]$Peak1_anno = ifelse(conns[[i]]$Peak1 %in% 
                                   all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
  conns[[i]]$Peak2_anno = ifelse(conns[[i]]$Peak2 %in% 
                                   all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
}
# select enhancer-promoter pairs
conns_6 = lapply(conns, function(x) x[x$Peak1_anno != x$Peak2_anno,])

# 6. find specific and common connections between cntrl and 5ht -----------------------------------------
conns_6$Ctrl$in_5ht <- compare_connections(conns_6$Ctrl, conns_6$`5HT`)
conns_6$`5HT`$in_cntrl <- compare_connections(conns_6$`5HT`, conns_6$Ctrl)

# ... select unique connection (Promoter - Enhancer)
conns_6.unique = lapply(conns_6, function(x) x[x$Peak1_anno == 'Promoter',])

# add gene symbol
conns_6.unique = lapply(conns_6.unique, function(x) merge(x, all.tCREs[,c(26,29)], by="Peak1", all.x=TRUE))

# 7. Visualization archplots (Supplementary figure 23B )------------------------------------------------------------------
# download and unzip annotation
temp <- tempfile()
download.file("http://ftp.ensembl.org/pub/release-109/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.109.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

conns_6.unique$Ctrl$conns_color = "#56B1F7"
conns_6.unique$Ctrl$conns_color[conns_6.unique$Ctrl$in_5ht] = "lightgray"

conns_6.unique$`5HT`$conns_color = "#fb7669"
conns_6.unique$`5HT`$conns_color[conns_6.unique$`5HT`$in_cntrl] = "lightgray"

# NFIX
plot_connections(conns_6.unique$Ctrl, 
                 "chr19", 23355498-200000, 23448265+100000,
                 viewpoint = "chr19_23354187_23356823",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 comparison_peak_color = "#fb7669",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.0,
                 comparison_coaccess_cutoff = 0.0,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr19", 23355498-200000, 23448265+100000,
                 viewpoint = "chr19_23354187_23356823",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$Ctrl,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.01,
                 comparison_coaccess_cutoff = 0.01,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

# NFIA
plot_connections(conns_6.unique$Ctrl, 
                 "chr5", 112439724-230000, 112443239+300000,
                 viewpoint = "chr5_112439724_112443239",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 comparison_peak_color = "#fb7669",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = 1.5,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr5", 112439724-230000, 112443239+300000,
                 viewpoint = "chr5_112439724_112443239",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$Ctrl,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = 1.5,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

# NFIB
plot_connections(conns_6.unique$Ctrl, 
                 "chr5", 96970365-270000, 96976509+230000,
                 viewpoint = "chr5_96970365_96976509",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 comparison_peak_color = "#fb7669",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.05,
                 comparison_coaccess_cutoff = 0.05,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr5", 96970365-270000, 96976509+230000,
                 viewpoint = "chr5_96970365_96976509",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$Ctrl,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.05,
                 comparison_coaccess_cutoff = 0.05,
                 connection_width = 3, 
                 collapseTranscripts = "longest")


# NFI motifs analysis unsing ChromVar (Supplementary figures 23C) ---------------------------------------------------------------
p28.signac = readRDS("./whole_assay/p28.signac.Rds")
levels(Idents(p28.signac))
DimPlot(p28.signac)
View(p28.signac@meta.data)
unique(p28.signac$stim)
p28.signac$stim[p28.signac$stim == "5HT"] = "5HTP"
p28.signac$stim = factor(p28.signac$stim, levels = c("Ctrl", "5HTP"))

library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(
  object = p28.signac,
  features = c("MA0670.1", "MA1643.1", "MA0671.1"), #A _ B _X 
  col = SpatialColors(n = 100)[51:100],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  split.by = "stim",
  label = F) & theme_void() &
  theme(
    legend.position = "none",       # Убираем легенду
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Заголовок по центру
  )  

# violin plots
srat.subset = subset(p28.signac, subset = celltype %in% "Tany")
srat.subset

# NFIA
# Get data for Violin plot
plot_data <- VlnPlot(srat.subset, 
                     assay = NULL, 
                     features = "MA0670.1", #NFIA
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data
plot_data$ident = "Tanycytes"
plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP"))
# Draw Violin plot with boxplot
library(ggpubr)

nfia_chrom = ggplot(plot_data, aes(x = group, y = MA0670.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  labs(title = "", y = "", x = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  # Убираем текст на оси X
        axis.ticks.x = element_blank(),  # Убираем отметки на оси X
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  # Указать соответствующие группы

# NFIB
# Получаем данные для Violin plot
plot_data <- VlnPlot(srat.subset, 
                     assay = NULL, 
                     features = "MA1643.1", #NFIB
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data
plot_data$ident = "Tanycytes"
plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP"))

# Создаем Violin plot с соответствующими boxplots
nfib_chrom = ggplot(plot_data, aes(x = group, y = MA1643.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  labs(title = "", y = "ChromVar score", x = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  # Убираем текст на оси X
        axis.ticks.x = element_blank(),  # Убираем отметки на оси X
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  # Указать соответствующие группы

# NFIX
# Получаем данные для Violin plot
plot_data <- VlnPlot(srat.subset, 
                     assay = NULL, 
                     features = "MA0671.1",  # NFIX
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data
plot_data$ident = "Tanycytes"
plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP"))
# Создаем Violin plot с соответствующими boxplots
nfix_chrom = ggplot(plot_data, aes(x = group, y = MA0671.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  labs(title = "", y = "", x = "Tanycytes") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  # Убираем текст на оси X
        axis.ticks.x = element_blank(),  # Убираем отметки на оси X
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  # Указать соответствующие группы



