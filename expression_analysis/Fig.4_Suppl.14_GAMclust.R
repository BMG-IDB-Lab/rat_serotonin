library(GAMclust)
library(gatom)
library(mwcsr)
library(fgsea)
library(data.table)
library(Seurat)

set.seed(42)

reduction <- "umap"

work.dir <- "work.dir"
dir.create(work.dir, showWarnings = F, recursive = T)
setwd(work.dir)

net_db <- "new"
topology <- "metabolites"
set.cor <- 0.4
start.base <- 0.5

# network ------------
if (net_db == "KEGG") {
  # KEGG:
  network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
  metabolites.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))
} else if (net_db == "Rhea") {
  # Rhea:
  network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rhea.rds"))
  metabolites.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.rhea.db.rds"))
} else if (net_db == "new") {
  network <- readRDS("network.combined.rds")
  metabolites.annotation <- readRDS("met.combined.db.rds")
}

# network anno ------------------
# makeOrgGatomAnnotation(org.db = org.Rn.eg.db)
network.annotation <- readRDS("org.Rn.eg.gatom.anno.rds")

# met to filter -----------------
met.to.filter <- data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID

# cplex -------------------------
cplex.dir <- "CPLEX_Studio1271"
solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

# stats -------------------------
stats.dir <- paste0(work.dir, "/stats")
dir.create(stats.dir, showWarnings = F, recursive = T)

# obj ---------------------------
seurat_object <- readRDS("scP28_olig.rds")
E <- as.matrix(Seurat::GetAssayData(object = seurat_object,
                                    assay = "SCT",
                                    slot = "scale.data"))

# prep --------------------------
E.prep <- prepareData(E = E,
                      gene.id.type = "Symbol",
                      use.PCA = TRUE,
                      use.PCA.n = 50,
                      network.annotation = network.annotation)

network.prep <- prepareNetwork(E = E.prep,
                               network = network,
                               met.to.filter = met.to.filter,
                               topology = topology,
                               network.annotation = network.annotation)

# pre clustering ----------------
cur.centers <- preClustering(E.prep = E.prep,
                             network.prep = network.prep,
                             initial.number.of.clusters = 32,
                             network.annotation = network.annotation)
pheatmap::pheatmap(
  GAMclust:::normalize.rows(cur.centers),
  cluster_rows=F, cluster_cols=F,
  show_rownames=F, show_colnames=T)

# clustering ----------------
results <- gamClustering(E.prep = E.prep,
                         network.prep = network.prep,
                         cur.centers = cur.centers,
                         
                         start.base = start.base,
                         base.dec = 0.05,
                         max.module.size = 50,
                         
                         cor.threshold = set.cor, 
                         p.adj.val.threshold = 1e-5,
                         
                         batch.solver = seq_batch_solver(solver),
                         work.dir = work.dir,
                         
                         show.intermediate.clustering = TRUE,
                         verbose = TRUE,
                         collect.stats = TRUE)

# results ----------------
## obtain graphs for each module ----
# Figure 4E, Supplementary Figure 14A
getGraphs(modules = results$modules,
          network.annotation = network.annotation,
          metabolites.annotation = metabolites.annotation,
          seed.for.layout = 42,
          work.dir = work.dir)
# further processing in Cytoscape
# get list of genes for each module
m.gene.list <- getGeneTables(modules = results$modules,
                             nets = results$nets,
                             patterns = results$patterns.pos,
                             gene.exprs = E.prep,
                             network.annotation = network.annotation,
                             work.dir = work.dir)

## visualize expression patterns of each module ----
# Supplementary Figure 14B
library(Seurat)
library(fgsea)

Seurat::DefaultAssay(seurat_object) <- "SCT"
geseca_list <- list()
for(i in seq_along(m.gene.list)){
  
  # print individuals:
  print(ifelse("SCT" %in% names(seurat_object@assays),
               "SCT", "RNA"))
  Seurat::DefaultAssay(seurat_object) <- ifelse("SCT" %in% names(seurat_object@assays),
                                                "SCT", "RNA")
  Seurat::FeaturePlot(seurat_object, 
                      slot = "data",
                      reduction = reduction,
                      features = m.gene.list[[i]],
                      ncol = 5,
                      order = T,
                      combine = T)
  
  ggplot2::ggsave(file = sprintf("%s/m.%s.genes.png", work.dir, i), 
                  width = 15, 
                  height = round(length(m.gene.list[[i]]) / 5)*2 + 1)
  
  # print pattern:
  print(ifelse("integrated" %in% names(seurat_object@assays),
               "integrated", DefaultAssay(seurat_object)))
  Seurat::DefaultAssay(seurat_object) <- ifelse("integrated" %in% names(seurat_object@assays),
                                                "integrated", Seurat::DefaultAssay(seurat_object))
  p <- fgsea::plotCoregulationProfileReduction(m.gene.list[[i]], seurat_object, 
                                               order = T,
                                               reduction = reduction) +
    ggplot2::ggtitle(paste0("Module_", i))
  
  ggplot2::ggsave(p,
                  file = sprintf("%s/m.%s.pattern.png", work.dir, i), 
                  width = 4.5, 
                  height = 4)
  
  geseca_list[[i]] <- p
  
} 

cowplot::plot_grid(plotlist = geseca_list,
                   ncol = 3)
ggplot2::ggsave(file = sprintf("%s/m.all.pattern.png", work.dir), 
                bg = "white",
                dpi = 400,
                width = 3 * 4, 
                height = 4)


## get annotations ----
# Figure 4K
getAnnotationTables(network.annotation = network.annotation,
                    nets = results$nets,
                    work.dir = work.dir)

# modification of getAnnotationHeatmap function to obtain dot plots instead of heatmap
# files with pathways
m_files <- list.files(work.dir, "m\\.[0-9]+\\.pathways\\.tsv", full.names = T)
m_files <- m_files[order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files)))]
# files with genes
m_files_genes <- list.files(work.dir, "m\\.[0-9]+\\.genes\\.tsv", full.names = T)
m_files_genes <- m_files_genes[order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files_genes)))]
m_files_genes <- m_files_genes[sort(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files)))]
# order
m_files <- m_files[c(1,3,2)]
m_files_genes <- m_files_genes[c(1,3,2)]

set <- list()
padj.threshold <- Inf
threshold <- 0.05
for (i in 1:length(m_files)) {
  GAMclust:::messagef("Processing module %s", as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files[i])))
  file <- GAMclust:::read.tsv(m_files[i])
  moduleSize <- nrow(GAMclust:::read.tsv(m_files_genes[i]))
  GAMclust:::messagef(sprintf("Module size: %s", moduleSize))
  if (length(file$overlap) == 0) {
    paths <- file[, c("pathway", "overlap")]
    paths$percentOfGenesInThePthw <- c(0)
    paths$PATHNAME <- c("No pathways")
  } else {
    paths <- file[which(file$padj < padj.threshold), ]
    paths <- paths[, c("pathway", "overlap")]
    paths <- paths[order(paths$overlap, decreasing = T), ]
    paths$percentOfGenesInThePthw <- paths$overlap / moduleSize
    paths <- paths[which(paths$percentOfGenesInThePthw > threshold), ]
  }
  paths <- mutate(paths, 
                  pathway = gsub(".*: ", "", pathway))
  set <- c(set, list(paths))
}

pathNames <-c()
for(i in seq_along(set)){
  pathNames <- c(pathNames, set[[i]]$pathway)
}
pathNames <- unique(pathNames)
set_df <- as.data.frame(matrix(NA, nrow = length(set), ncol=length(pathNames)))
colnames(set_df) <- pathNames
for(i in seq_along(set)){
  set_df[i, set[[i]]$pathway] <- set[[i]]$percentOfGenesInThePthw
}
set_df[is.na(set_df)] <- 0
rownames(set_df) <- c(paste("Module", as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files))))

plot_df <- set_df %>% 
  tibble::rownames_to_column("Module") %>% 
  tidyr::gather(key = "Pathway", value = "Percent", -Module) %>% 
  mutate(Module = factor(Module, levels = paste0("Module ", c(1,3,2)),
                         labels = c(1,3,2)),
         Pathway = factor(Pathway, levels = rev(colnames(set_df)),
                          labels = gsub(", .*", "", rev(colnames(set_df)))))
plot_df[plot_df == 0] <- NA

colors4dotplot <- c("white", 
                    "purple4", "#5E4FA2", "#486ecf", "#3288BD", "#66adc2", 
                    "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#fad56b", 
                    "#f7bd16", "#FDAE61", "#fa912a", "#F46D43", "#D53E4F", 
                    "#9E0142", "#8B0001", "#830001", "#720000", "#690000")
ggplot(plot_df, aes(x = Module, y = Pathway, color = Percent)) +
  geom_point(size = 10) +
  scale_color_gradientn(colours = colors4dotplot, limits = c(0,1), na.value = NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black")) +
  scale_y_discrete(labels = scales::label_wrap(40))
