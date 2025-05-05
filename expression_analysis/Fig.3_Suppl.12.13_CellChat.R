source("R-scripts_paper/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(CellChat)
library(grid)
library(ComplexHeatmap)
library(RColorBrewer)

# Data ----
p28_sn_n <- readRDS("snP28_neurons.rds") 
p28_sn_n$sample_sex <- paste0(p28_sn_n$sample, "_", p28_sn_n$sex_predicted)
DefaultAssay(p28_sn_n) <- "SCT"

p28_sn <- readRDS("snP28_whole.rds")
DefaultAssay(p28_sn) <- "RNA"


# Analysis of neuronal communication ----

res_select <- "res.combined_det_markers"

x <- .Random.seed
set.seed(123)
colors_orig <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
                 RColorBrewer::brewer.pal(n = 8, "Accent")[-4],
                 RColorBrewer::brewer.pal(n = 8, "Set1")[-6], rainbow(18), sample(rainbow(50), 
                                                                                  size = 50, replace = F)
)[-c(16,19,22)] %>% unique()
.Random.seed <- x
colors_n_here <- colors_orig[1:59]

# exclude small cluster and cluster with mixed dopamine/choline markers
cl_rm <- (table(p28_sn_n@meta.data[, c("res.combined_det_markers", "treatment")]) %>% 
            as.data.frame() %>% 
            mutate(treatment = ifelse(treatment == "5HT", "HT", "Ctrl")) %>% 
            tidyr::spread(key = "treatment", value = "Freq") %>% 
            mutate(less10 = ifelse(HT < 10 | Ctrl < 10, T, F)) %>% 
            filter(less10 == TRUE))$res.combined_det_markers %>% as.character()
cl_rm <- c(cl_rm, "51_Dop_Chol_Eef1b2.Cd59b")
cl_rm
p28_sn_n <- subset(p28_sn_n, res.combined_det_markers %in% cl_rm, invert = T)

## CellChat analysis separately for control/5HTP samples ----
data_list <- SplitObject(p28_sn_n, split.by = "treatment")
res <- res_select
for (i in 1:length(data_list)) {
  Idents(data_list[[i]]) <- res
  print(DimPlot(data_list[[i]], label = T) + ggtitle(names(data_list[i])))
}

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB

cellchat_list <- list()
for (i in 1:length(data_list)) {
  data.input <- GetAssayData(data_list[[i]], assay = "SCT", slot = "data")
  meta <- data_list[[i]]@meta.data[, c(res, "sample"), drop = F]
  colnames(meta) <- c("labels", "samples")
  meta$samples <- factor(meta$samples, levels = unique(meta$sample))

  meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),
                                                          unique(meta$labels)))

  cellchat_i <- createCellChat(object = data.input, meta = meta, group.by = "labels")

  cellchat_i <- addMeta(cellchat_i, meta = meta)
  cellchat_i <- setIdent(cellchat_i, ident.use = "labels") # set "labels" as default cell identity
  groupSize <- as.numeric(table(cellchat_i@idents))

  # set the used database in the object
  cellchat_i@DB <- CellChatDB.use

  cellchat_i <- subsetData(cellchat_i) # This step is necessary even if using the whole database

  cellchat_i <- identifyOverExpressedGenes(cellchat_i,
                                           thresh.pc = 0.25, 
                                           thresh.fc = 0.25)
  cellchat_i <- identifyOverExpressedInteractions(cellchat_i)

  cellchat_i <- computeCommunProb(cellchat_i)
  cellchat_i <- filterCommunication(cellchat_i, 
                                    min.cells = 10)

  cellchat_i <- computeCommunProbPathway(cellchat_i)

  cellchat_i <- aggregateNet(cellchat_i)

  cellchat_i <- netAnalysis_computeCentrality(cellchat_i, slot.name = "netP")

  cellchat_list[[i]] <- cellchat_i

  names(cellchat_list)[i] <- names(data_list)[i]
}

cellchat_list <- list(cellchat_list[["Ctrl"]],
                      cellchat_list[["5HT"]])
names(cellchat_list) <- c("Ctrl", "5HT")
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))


## Communication flow ----

# Figure 3D
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, 
                  do.flip = F,
                  color.use = c("#56B1F7", "#fb7669")) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 14))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, 
                  do.flip = F,
                  color.use = c("#56B1F7", "#fb7669")) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 14))
(gg1 + theme(plot.margin = margin(60,2,2,2))) / (gg2 + theme(plot.margin = margin(5,2,2,2)))



## Hypocretin signaling ----

cc1 <- aggregateNet_my(cellchat_list[["5HT"]], signaling = "HCRT", 
                       return.object = T, remove.isolate = F)
cc2 <- aggregateNet_my(cellchat_list[["Ctrl"]], signaling = "HCRT", 
                       return.object = T, remove.isolate = F)

# Figure 3E (left)
cellchat_p <- mergeCellChat(list(cc2, cc1), add.names = c("Ctrl", "5HT"))
cellchat_p <- renameClusters(cellchat_p, keep_parts = 1)
names(colors_n_here) <- c(0:(length(colors_n_here)-1))
cl <- "50"
netVisual_diffInteraction_my(cellchat_p, weight.scale = T, measure = "weight",
                             vertex.label.cex = 1.2,
                             margin = 0.1, edge.curved = F,
                             label.edge = F, 
                             color.edge = rev(c("#56B1F7", "#fb7669")),
                             color.use = colors_n_here, 
                             title.name = "",
                             sources.use = cl, remove.isolate = T)

# Figure 3E (right)
cellchat_p <- mergeCellChat(list(cc2, cc1), add.names = c("Ctrl", "5HT"))
names(colors_n_here) <- stringr::str_sort(levels(cellchat_p@meta$labels), numeric = T)
cl <- "50_Glut_Hcrt"
netVisual_heatmap_my(cellchat_p, sources.use = cl,
                     font.size = 14, font.size.title = 18,
                     measure = "weight",
                     color.use = colors_n_here,
                     remove.isolate = T,
                     cluster.rows = T, cluster.cols = T)


## Overall in/outgoing neuronal signaling ----

### Heatmaps ----

#### Ingoing ----

in_ <- cbind(
  data.frame(Ctrl = rowSums(sapply(cellchat@netP[["Ctrl"]][["centr"]], function (x) x$indeg))),
  data.frame(HT = rowSums(sapply(cellchat@netP[["5HT"]][["centr"]], function (x) x$indeg)))) %>% 
  tibble::rownames_to_column("cluster") %>% 
  mutate(abs_diff = abs(HT - Ctrl),
         diff = HT - Ctrl,
         diff_ratio = ifelse(HT > Ctrl, HT / Ctrl, -(Ctrl / HT))) %>% 
  tidyr::gather(key = "condition", value = "prob", -c("cluster", diff, abs_diff, diff_ratio)) %>% 
  mutate(condition_x = ifelse(condition == "Ctrl", -1, 1),
         condition = ifelse(condition == "HT", "5HTP", "Ctrl"))

in_df1 <- as.data.frame(sapply(cellchat@netP[["5HT"]][["centr"]], function (x) x$indeg)) 
in_df2 <- as.data.frame(sapply(cellchat@netP[["Ctrl"]][["centr"]], function (x) x$indeg))
diff_order <- (in_ %>% arrange(-diff))$cluster %>% unique() %>% as.character()

# Get the union of column names from both dataframes
all_cols <- union(names(in_df1), names(in_df2))

# Process in_df1
missing_cols1 <- setdiff(all_cols, names(in_df1))
in_df1[missing_cols1] <- 0
in_df1 <- in_df1[diff_order,all_cols]

# Process in_df2
missing_cols2 <- setdiff(all_cols, names(in_df2))
in_df2[missing_cols2] <- 0
in_df2 <- in_df2[diff_order,all_cols]

in_df <- in_df1 - in_df2


#### Outgoing ----

out <- cbind(
  data.frame(Ctrl = rowSums(sapply(cellchat@netP[["Ctrl"]][["centr"]], function (x) x$outdeg))),
  data.frame(HT = rowSums(sapply(cellchat@netP[["5HT"]][["centr"]], function (x) x$outdeg)))) %>% 
  tibble::rownames_to_column("cluster") %>% 
  mutate(abs_diff = abs(HT - Ctrl),
         diff = HT - Ctrl,
         diff_ratio = ifelse(HT > Ctrl, HT / Ctrl, -(Ctrl / HT))) %>% 
  tidyr::gather(key = "condition", value = "prob", -c("cluster", diff, abs_diff, diff_ratio)) %>% 
  mutate(condition_x = ifelse(condition == "Ctrl", -1, 1),
         condition = ifelse(condition == "HT", "5HTP", "Ctrl"))

out_df1 <- as.data.frame(sapply(cellchat@netP[["5HT"]][["centr"]], function (x) x$outdeg)) 
out_df2 <- as.data.frame(sapply(cellchat@netP[["Ctrl"]][["centr"]], function (x) x$outdeg))
diff_order <- (out %>% arrange(-diff))$cluster %>% unique() %>% as.character()

# Get the union of column names from both dataframes
all_cols <- union(names(out_df1), names(out_df2))

# Process out_df1
missing_cols1 <- setdiff(all_cols, names(out_df1))
out_df1[missing_cols1] <- 0
out_df1 <- out_df1[diff_order,all_cols]

# Process out_df2
missing_cols2 <- setdiff(all_cols, names(out_df2))
out_df2[missing_cols2] <- 0
out_df2 <- out_df2[diff_order,all_cols]

out_df <- out_df1 - out_df2


#### Plot both ----

min_both <- min(min(in_df), min(out_df))
max_both <- max(max(in_df), max(out_df))

#### Supplementary Figure 12D (Ingoing)
top_diff_abs <- unique(slice_max(in_, order_by = abs_diff, n = 30)$cluster) %>% 
  as.character()
ht1 <- ComplexHeatmap::Heatmap(in_df[top_diff_abs, !colSums(in_df[top_diff_abs,]) == 0],
                               cluster_rows = T,
                               col = colorRamp3(c(min_both, 0, max_both), 
                                                c('#2166ac', "#f7f7f7", '#b2182b')),
                               na_col = "white",
                               column_names_rot = 45,
                               row_names_side = "left",row_names_rot = 0,
                               row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7)
); ht1

#### Supplementary Figure 12C (Outgoing)
top_diff_abs <- unique(slice_max(out, order_by = abs_diff, n = 30)$cluster) %>% 
  as.character()
ht1 <- ComplexHeatmap::Heatmap(out_df[top_diff_abs, !colSums(out_df[top_diff_abs,]) == 0],
                               cluster_rows = T,
                               col = colorRamp3(c(min_both, 0, max_both), 
                                                c('#2166ac', "#f7f7f7", '#b2182b')),
                               na_col = "white",
                               column_names_rot = 45,
                               row_names_side = "left",row_names_rot = 0,
                               row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7)
); ht1


### Dot plots ----

colors <- scPalette(61)

#### Ingoing ----
# Supplementary Figure 12B

top_diff <- unique(slice_max(in_, order_by = abs_diff, n = 30)$cluster) %>% 
  as.character()
in_$cluster <- factor(in_$cluster,
                      levels = stringr::str_sort(unique(in_$cluster), numeric = T))
in_$alpha <- ifelse(in_$cluster %in% top_diff, 1, 0.5)


ggplot(in_, aes(x = condition_x, y = prob, 
                color = cluster, fill = cluster, alpha = alpha, 
                label = cluster)) +
  geom_line(aes(group = cluster), linewidth = 1, show.legend = FALSE) +
  geom_point(data = in_[order(in_$abs_diff), ], size = 5) +
  guides(color = "none", fill = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(
    color = "none",  # No legend for color
    fill = "none",   # No legend for fill
    alpha = guide_legend(override.aes = list(shape = 16, size = 5, linetype = 0)) 
  ) +
  ggrepel::geom_text_repel(
    data = subset(in_, condition_x == -1 & cluster %in% top_diff),
    size = 4,
    aes(x = condition_x - 0.5, y = prob),  # Adjust 0.5 to control leftward distance
    hjust = 1,  # Right-align text to ensure labels' right edges align
    max.overlaps = 100, min.segment.length = 100,
    direction = "y",
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = subset(in_, condition_x == 1 & cluster %in% top_diff),
    size = 4,
    aes(x = condition_x + 0.5, y = prob),  
    hjust = 0,  
    max.overlaps = 100, min.segment.length = 100,
    direction = "y",
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  expand_limits(x = c(-5, 5))+
  theme(panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_alpha_continuous(
    name = "Cluster Impact",  # Legend title
    breaks = c(0.5, 1),      # Only show alpha values 0.5 and 1
    labels = c("Others", "Most Affected")  # Custom labels
  ) +
  # Add "Ctrl" and "5HTP" labels at the bottom
  annotate("text", x = -1, y = 0 - 0.05, 
           label = "Ctrl", size = 7, vjust = 1) +
  annotate("text", x = 1, y = 0 - 0.05, 
           label = "5HTP", size = 7, vjust = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = Inf, 
           color = "black", linewidth = 0.5) +
  theme(legend.position = "bottom") +
  ggtitle("Incoming signalling") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))


#### Outgoing ----
# Supplementary Figure 12A

top_diff <- unique(slice_max(out, order_by = abs_diff, n = 30)$cluster) %>% 
  as.character()

out$cluster <- factor(out$cluster,
                      levels = stringr::str_sort(unique(out$cluster), numeric = T))
out$alpha <- ifelse(out$cluster %in% top_diff, 1, 0.5)


ggplot(out, aes(x = condition_x, y = prob, 
                color = cluster, fill = cluster, alpha = alpha, 
                label = cluster)) +
  geom_line(aes(group = cluster), linewidth = 1, show.legend = FALSE) +
  geom_point(data = out[order(out$abs_diff), ], size = 5) +
  guides(color = "none", fill = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(
    color = "none",  # No legend for color
    fill = "none",   # No legend for fill
    alpha = guide_legend(override.aes = list(shape = 16, size = 5, linetype = 0)) 
  ) +
  ggrepel::geom_text_repel(
    data = subset(out, condition_x == -1 & cluster %in% top_diff),
    size = 4,
    aes(x = condition_x - 0.5, y = prob),  # Adjust 0.5 to control leftward distance
    hjust = 1,  # Right-align text to ensure labels' right edges align
    max.overlaps = 100, min.segment.length = 100,
    direction = "y",
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = subset(out, condition_x == 1 & cluster %in% top_diff),
    size = 4,
    aes(x = condition_x + 0.5, y = prob),
    hjust = 0,  
    max.overlaps = 100, min.segment.length = 100,
    direction = "y",
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  expand_limits(x = c(-5, 5))+
  theme(panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_alpha_continuous(
    name = "Cluster Impact",  # Legend title
    breaks = c(0.5, 1),      # Only show alpha values 0.5 and 1
    labels = c("Others", "Most Affected")  # Custom labels
  ) +
  # Add "Ctrl" and "5HTP" labels at the bottom
  annotate("text", x = -1, y = 0 - 0.05, 
           label = "Ctrl", size = 7, vjust = 1) +
  annotate("text", x = 1, y = 0 - 0.05, 
           label = "5HTP", size = 7, vjust = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = Inf, 
           color = "black", linewidth = 0.5) +
  theme(legend.position = "bottom") +
  ggtitle("Outgoing signalling") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))



## Nrxn-Nlgn signaling ----

### Nrxn3 ----
# Supplementary Figure 12E

LR.pair <- data.frame(interaction_name = "NRXN3_NLGN1") 
cc1 <- aggregateNet_my(cellchat_list[["5HT"]], pairLR.use = LR.pair, 
                       return.object = T, remove.isolate = F)
cc2 <- aggregateNet_my(cellchat_list[["Ctrl"]], pairLR.use = LR.pair, 
                       return.object = T, remove.isolate = F)
cellchat_p <- mergeCellChat(list(cc2, cc1), add.names = c("Ctrl", "5HT"))
cellchat_p <- renameClusters(cellchat_p, keep_parts = 2)

w <- cellchat_p@net[["5HT"]][["weight"]] - cellchat_p@net[["Ctrl"]][["weight"]] %>% 
  as.matrix()
color.heatmap.use <- colorRamp3(c(min(w), 0, max(w)), c('#2166ac',"#f7f7f7",'#b2182b'))
ht <- ComplexHeatmap::Heatmap(w,
                              col = color.heatmap.use,
                              na_col = "white",
                              name = "Communication Prob. (relative)",
                              cluster_rows = F,
                              cluster_columns = F,
                              row_names_side = "left",  
                              column_names_side = "bottom",
                              column_title = "NRXN3-NLGN1",
                              column_names_gp = grid::gpar(fontsize = 7),
                              row_names_gp = grid::gpar(fontsize = 7),
                              heatmap_legend_param = list(
                                title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 8)
                              )); ht


### Nrxn1 ----
# Supplementary Figure 12F

LR.pair <- data.frame(interaction_name = "NRXN1_NLGN1") 
cc1 <- aggregateNet_my(cellchat_list[["5HT"]], pairLR.use = LR.pair, 
                       return.object = T, remove.isolate = F)
cc2 <- aggregateNet_my(cellchat_list[["Ctrl"]], pairLR.use = LR.pair, 
                       return.object = T, remove.isolate = F)
cellchat_p <- mergeCellChat(list(cc2, cc1), add.names = c("Ctrl", "5HT"))
cellchat_p <- renameClusters(cellchat_p, keep_parts = 2)

w <- cellchat_p@net[["5HT"]][["weight"]] - cellchat_p@net[["Ctrl"]][["weight"]] %>% 
  as.matrix()
color.heatmap.use <- colorRamp3(c(min(w), 0, max(w)), c('#2166ac',"#f7f7f7",'#b2182b'))
ht <- ComplexHeatmap::Heatmap(w,
                              col = color.heatmap.use,
                              na_col = "white",
                              name = "Communication Prob. (relative)",
                              cluster_rows = F,
                              cluster_columns = F,
                              row_names_side = "left",  
                              column_names_side = "bottom",
                              column_title = "NRXN1-NLGN1",
                              column_names_gp = grid::gpar(fontsize = 7),
                              row_names_gp = grid::gpar(fontsize = 7),
                              heatmap_legend_param = list(
                                title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 8)
                              )); ht


## NTS signaling ----
# Supplementary Figure 12G

cc1 <- aggregateNet_my(cellchat_list[["5HT"]], signaling = "NTS", 
                       return.object = T, remove.isolate = F)
cc2 <- aggregateNet_my(cellchat_list[["Ctrl"]], signaling = "NTS", 
                       return.object = T, remove.isolate = F)
cellchat_p <- mergeCellChat(list(cc2, cc1), add.names = c("Ctrl", "5HT"))

netVisual_heatmap_my(cellchat_p, signaling = "NTS",
                     remove.isolate = T,
                     cluster.rows = T, cluster.cols = T)



# Analysis of communication between major cell types ----

x <- .Random.seed
set.seed(123)
colors <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
            RColorBrewer::brewer.pal(n = 8, "Accent")[-4],
            RColorBrewer::brewer.pal(n = 8, "Set1")[-6], sample(rainbow(31)))[-c(16,19,22,23)]
.Random.seed <- x
colors_n <- c(colors[1],"palegreen2", "darkseagreen4", "lightblue3",
              colors[-c(1,8,9,10,13)])


p28_sn <- subset(p28_sn, Cell_type %in% c("MG_MF", "VEnC_Mural", "Lepto", "PTub"), 
                 invert = T)
DimPlot(p28_sn, group.by = "Cell_type")

# add info on broad neuronal cell types and filter cells

cells <- names(p28_sn_n$res.combined_det_markers)
neuro_cl <- setNames(gsub("(\\d+_.*)_.*", "\\1", p28_sn_n$res.combined_det_markers),
                     cells)
cells <- names(neuro_cl)
neuro_cl <- gsub("\\d+_", "", neuro_cl); names(neuro_cl) <- cells
cells <- names(p28_sn$Cell_type)
all <- as.character(p28_sn$Cell_type); names(all) <- cells
all <- all[which(names(all) %not in% names(neuro_cl))]
new_ct <- c(all, neuro_cl)
p28_sn <- AddMetaData(p28_sn, new_ct, "Cell_type")
DimPlot(p28_sn, group.by = "Cell_type", label = T)

Idents(p28_sn) <- "Cell_type"
p28_sn <- subset(p28_sn, Cell_type == "Dop_Chol", invert = T)
p28_sn$Cell_type <- factor(p28_sn$Cell_type,
                           levels = c('GABA', 'Glut', 'Dop', 'Chol',
                                      'mOD_1', 'mOD_2','mOD_3',
                                      'OPC/COP_1', 'OPC/COP_2',
                                      'Astrocytes',
                                      'Ependyma', 'Tany'))
DimPlot(p28_sn, group.by = "Cell_type", label = T)

## CellChat ----

data_list <- SplitObject(data, split.by = "treatment")
res <- "Cell_type"
for (i in 1:length(data_list)) {
  data_list[[i]] <- NormalizeData(data_list[[i]], assay = "RNA")
  Idents(data_list[[i]]) <- res
  print(DimPlot(data_list[[i]], label = T) + ggtitle(names(data_list[i])))
}

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB

cellchat_list <- list()
for (i in 1:length(data_list)) {
  data.input <- GetAssayData(data_list[[i]], assay = "RNA", slot = "data")
  meta <- data_list[[i]]@meta.data[, c(res, "sample"), drop = F]
  colnames(meta) <- c("labels", "samples")

  cellchat_i <- createCellChat(object = data.input, meta = meta, group.by = "labels")

  cellchat_i <- addMeta(cellchat_i, meta = meta)
  cellchat_i <- setIdent(cellchat_i, ident.use = "labels") # set "labels" as default cell identity
  groupSize <- as.numeric(table(cellchat_i@idents))

  # set the used database in the object
  cellchat_i@DB <- CellChatDB.use

  cellchat_i <- subsetData(cellchat_i) # This step is necessary even if using the whole database

  cellchat_i <- identifyOverExpressedGenes(cellchat_i)
  cellchat_i <- identifyOverExpressedInteractions(cellchat_i)

  cellchat_i <- computeCommunProb(cellchat_i)
  cellchat_i <- filterCommunication(cellchat_i, 
                                    min.cells = 10)

  cellchat_i <- computeCommunProbPathway(cellchat_i)

  cellchat_i <- aggregateNet(cellchat_i)

  cellchat_i <- netAnalysis_computeCentrality(cellchat_i, slot.name = "netP")

  cellchat_list[[i]] <- cellchat_i

  names(cellchat_list)[i] <- names(data_list)[i]
}

cellchat_list <- list(cellchat_list[["Ctrl"]],
                      cellchat_list[["5HT"]])
names(cellchat_list) <- c("Ctrl", "5HT")
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))


# perform de for LR pairs filtering and more robust results
pos.dataset = "5HT"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset, features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1,
                                       thresh.fc = 0.05,thresh.p = 0.05,
                                       group.DE.combined = FALSE)
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "5HT",
                              ligand.logFC = 0.05, receptor.logFC = 0.05)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Ctrl",
                                ligand.logFC = -0.05, receptor.logFC = -0.05)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]


### Compare overall communication flow ----

# Figure 3F (left)
netVisual_diffInteraction(cellchat, weight.scale = T, 
                          color.edge = c("#fb7669", "#56B1F7"),
                          color.use = colors_n,
                          edge.width.max = 17,
                          title.name = "")

# Supplementary Figure 13B
colors_n_here <- colors_n[1:length(cellchat@meta$labels %>% levels())]
names(colors_n_here) <- levels(cellchat@meta$labels)
gg1 <- netVisual_heatmap(cellchat, color.heatmap = rev(c("#fb7669", "#56B1F7")), 
                         color.use = colors_n_here)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         color.heatmap = rev(c("#fb7669", "#56B1F7")),
                         color.use = colors_n_here)
gg1 + gg2


### Tanycytes and Neurons communication ----

# Main affected pathways 
# Figure 3F (right)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,
        do.flip = T, 
        sources.use = c("Tany", "Ependyma"),
        targets.use = c("GABA", "Glut", "Dop", "Chol"),
        color.use = c("#56B1F7", "#fb7669"))

# Closer look at LR pairs 
# Supplementary Figure 13C
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                 sources.use = c("Tany", "Ependyma"),
                 targets.use = c('GABA', 'Glut', 'Dop', 'Chol'),
                 comparison = c(1, 2),  angle.x = 90, max.dataset = which(names(cellchat_list) == "5HT"),
                 remove.isolate = T,
                 color.text =  c("#56B1F7", "#fb7669"),
                 title.name = paste0("Up-regulated signaling in ", pos.dataset))

# Nrxn signaling LR pairs
# Figure 3G (left)
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up %>% 
                   filter(grepl("NRXN", interaction_name)), 
                 sources.use = c("Tany", "Ependyma"),
                 targets.use = c('GABA', 'Glut', 'Dop', 'Chol'),
                 comparison = c(1,2),  angle.x = 90, max.dataset = which(names(cellchat_list) == "5HT"),
                 remove.isolate = T,
                 color.text =  c("#56B1F7", "#fb7669"),
                 title.name = paste0("Up-regulated signaling in ", pos.dataset)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(40,2,2,20))


### OPC/COP and Neurons communication ----

# Main affected pathways 
# Supplementary Figure 13D
gg1 <- rankNet_my(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,
                  do.flip = F,
                  sources.use = "OPC/COP_2", targets.use = c("GABA", "Glut"),
                  color.use = c("#56B1F7", "#fb7669"))
gg2 <- rankNet_my(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,
                  do.flip = F,
                  sources.use = "OPC/COP_2", targets.use = c("GABA", "Glut"),
                  color.use = c("#56B1F7", "#fb7669"))
gg1 / gg2


# Closer look at LR pairs 
# Supplementary Figure 13E
netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                 sources.use = "OPC/COP_2",
                 targets.use = c("GABA", "Glut"),
                 comparison = c(1,2),  angle.x = 90, max.dataset = which(names(cellchat_list) == "Ctrl"),
                 remove.isolate = T,
                 color.text =  c("#56B1F7", "#fb7669"),
                 title.name = paste0("Down-regulated signaling in ", pos.dataset))

# Cntn signaling LR pairs
# Figure 3G (right)
netVisual_bubble(cellchat, pairLR.use = pairLR.use.down %>% 
                   filter(grepl("CNTN", interaction_name)), 
                 sources.use = "OPC/COP_2",
                 targets.use = c("GABA", "Glut"),
                 comparison = c(1,2),  angle.x = 90, max.dataset = which(names(cellchat_list) == "Ctrl"),
                 remove.isolate = T,
                 color.text =  c("#56B1F7", "#fb7669"),
                 title.name = paste0("Down-regulated signaling in ", pos.dataset)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(40,2,2,20))
