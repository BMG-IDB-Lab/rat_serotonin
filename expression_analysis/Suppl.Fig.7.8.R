source("R-scripts_paper/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)
`%not in%` <- Negate(`%in%`)

colors <- c(RColorBrewer::brewer.pal(n = 8, "Set2")[-7],
            RColorBrewer::brewer.pal(n = 8, "Accent")[-4],
            RColorBrewer::brewer.pal(n = 8, "Set1")[-6], rainbow(31), rainbow(80))[-c(8,9,16,19,22,23)]

# Data ----
# data with rat gene names changed to corresponding mouse genes where possible
p28_sn <- readRDS("snP28_whole_mm38.rds")
p28_sn <- NormalizeData(p28_sn, assay = "RNA")
p28_sn <- ScaleData(p28_sn, assay = "RNA")

p28_sn_n <- readRDS("snP28_neurons_mm38.rds") 
p28_sn_n <- NormalizeData(p28_sn_n, assay = "RNA")
p28_sn_n <- ScaleData(p28_sn_n, assay = "RNA")

# Mapping of rat cells to other datasets ----
# preprocessed data from open sources
# choose one option
ref <- readRDS("influenzaA.rds")
ref <- readRDS("dentate_gyrus.rds")

DefaultAssay(p28_sn) <- "RNA"
DefaultAssay(ref) <- "integrated"

# transfer labels
anchors <- FindTransferAnchors(reference = ref, query = p28_sn, dims = 1:30,
                               reference.reduction = "pca", 
                               normalization.method = "SCT", 
                               query.assay = "RNA", 
                               reference.assay = "integrated")
# Influenza = cluster_id / Dentate gyrus = cluster_name
predictions <- TransferData(anchorset = anchors, refdata = ref$cluster_id, dims = 1:30)
predictions <- predictions %>% 
  mutate(predicted.id_mod = ifelse(prediction.score.max >= 0.5, predicted.id, "Unassigned"))

p28_sn <- AddMetaData(p28_sn, metadata = predictions)

# Influenza A
p28_sn$predicted.id_mod2 <- ifelse(grepl("Gaba_|Glut", p28_sn$predicted.id_mod),
                                                  gsub("_\\d+", "", p28_sn$predicted.id_mod),
                                                  p28_sn$predicted.id_mod)
cl_un <- which((unique(p28_sn$predicted.id_mod2) %>% sort()) == "Unassigned")
DimPlot(p28_sn, 
        group.by = "predicted.id_mod2", 
        reduction = red,
        pt.size = 0.4,
        cols = swap_elements(colors, 7, cl_un))

# Dentate Gyrus
cl_un <- which((unique(p28_sn$predicted.id_mod) %>% sort()) == "Unassigned")
DimPlot(p28_sn, 
        group.by = "predicted.id_mod", 
        reduction = red,
        pt.size = 0.4,
        cols = swap_elements(colors, 7, cl_un))



# HypoMap ----
gdata::mv("p28_sn_n", "query")
query <- mapscvi::map_new_seurat_hypoMap(
  query, 
  suffix="whole",
  label_col = "C286_named",
  max_epochs=20)

# set.seed(785065)
# colors_hypomap <- unname(Polychrome::createPalette(300,  colors[c(1:7)]))
# 
# # Define the target color in hex
# target_color <- colors[7]
# # Convert the target color to RGB
# target_rgb <- col2rgb(target_color)
# # Define a similarity threshold (adjust as needed)
# threshold <- 50
# # Function to calculate Euclidean distance between two colors
# color_distance <- function(hex1, hex2) {
#   rgb1 <- col2rgb(hex1)
#   rgb2 <- col2rgb(hex2)
#   sqrt(sum((rgb1 - rgb2)^2))
# }
# # Calculate distances and filter out similar colors
# filtered_colors <- colors_hypomap[sapply(color_hypomap, function(hex) {
#   color_distance(hex, target_color) > threshold
# })]
# # Print the filtered colors
# print(filtered_colors)
# saveRDS(filtered_colors, file = "colors_hypomap.rds")
colors_hypomap <- readRDS("colors_hypomap.rds")

new_meta <- query@meta.data[, c("predicted", "prediction_probability")] %>% 
  mutate(confidence = ifelse(prediction_probability >= 0.5, T, F),
         new_id = ifelse(confidence == TRUE, predicted, "Unassigned"))
query <- AddMetaData(query, setNames(factor(new_meta$new_id,
                                            levels = rev(stringr::str_sort(unique(new_meta$new_id), numeric = T))), 
                                     rownames(new_meta)), "predicted_0.5")

cl_un <- which((unique(query$predicted_0.5) %>% sort()) == "Unassigned")
colors_hypomap[cl_un] <- colors[7]
DimPlot(query, reduction = red, group.by = "predicted_0.5",
        cols = colors_hypomap,
        order = T,
        pt.size = 1.5) +
  guides(color=guide_legend(nrow=35,byrow=TRUE,override.aes = list(size=4)))

# only if cell number > 90
keep <- (table(query$predicted_0.5) %>% as.data.frame() %>% filter(Freq >= 90))$Var1 %>% 
  as.character()
plot <- DimPlot(query, reduction = red, group.by = "predicted_0.5",
                cols = colors_hypomap,
                order = T,
                label = F,
                pt.size = 1.5) + NoLegend()
LabelClusters(plot, id = "predicted_0.5", clusters = keep,
              repel = T,
              size = 5,
              fontface = "bold")
