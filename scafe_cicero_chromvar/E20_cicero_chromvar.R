"""
This script reproduce analysis of connections between proximal and distal tCREs 
for hypothalamic samples of Rattus Norvegicus at day 20 of embryo development (E20 stage).
Using Cicero we reconstructed connections between proximal and distal tCREs for control and 5HTP-stimulated samples separately 
and compared the reconstructed conections between conditions.
Download .zip archive from 10.5281/zenodo.15237638, unpack it and make folder <scafe_cicero_chromvar> as working directory.
"""
setwd("~/scafe_cicero_chromvar/")

set.seed(1234)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cicero)
library(dplyr)
library(data.table)

library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# 1 Upload E20 integrated SeuratObject containig peak x cell count matrix after running SCAFE --------------

# upload E20 SeuratObject
e20.scafe = readRDS("~/e20.scafe.Rds")
Idents(e20.scafe) = e20.scafe@meta.data$stim
DefaultAssay(e20.scafe) = "RNA"

# split by stim
unique(e20.scafe@meta.data$stim)
e20.scafe.list = list()
for(st in unique(e20.scafe@meta.data$stim)){
  print(st)
  v = subset(e20.scafe, subset = stim == st)
  print(dim(v))
  print(unique(v@meta.data$stim))
  e20.scafe.list = base::append(e20.scafe.list, v)
}
names(e20.scafe.list) = unique(e20.scafe@meta.data$stim)
e20.scafe.list
sapply(e20.scafe.list, DefaultAssay)
rm(v)

# 2 Convert to cell_data_set object  ------------------------------------------
cds.list <- lapply(e20.scafe.list, function(x) as.cell_data_set(x))
saveRDS(cds.list, "./serotonin6.cds.stim.Rds")

# since it misses the gene_short_name column, let's add it
for(i in 1:length(cds.list)){
  fData(cds.list[[i]])$gene_short_name <- rownames(cds.list[[i]])
}
lapply(cds.list, function(x) fData(x))
# to get counts
lapply(cds.list, counts)

# 2.1 Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have
recreate.partitions = lapply(cds.list, function(x) 
  recreate.partition = c(rep(1,length(x@colData@rownames))))
for(i in 1:length(cds.list)){
  names(recreate.partitions[[i]]) = cds.list[[i]]@colData@rownames
  recreate.partitions[[i]] = as.factor(recreate.partitions[[i]])}
for(i in 1:length(cds.list)){
  cds.list[[i]]@clusters$UMAP$partitions <- recreate.partitions[[i]]}
# Assign the cluster info 
lapply(e20.scafe.list, function(x) x@active.ident)
list_clusters = lapply(e20.scafe.list, function(x) x@active.ident)
for(i in 1:length(cds.list)){
  cds.list[[i]]@clusters$UMAP$clusters <- list_clusters[[i]]}
# Assign UMAP coordinate - cell embeddings
for(i in 1:length(cds.list)){
  cds.list[[i]]@int_colData@listData$reducedDims$UMAP  = e20.scafe.list[[i]]@reductions$umap@cell.embeddings}


# 3 make cicero cds. ----------------------------------------
umap_coords = lapply(cds.list, function(x) reducedDims(x)$UMAP)
cicero_cds = list()

for(i in 1:length(cds.list)){
  cicero_cds[[i]] =  make_cicero_cds(cds.list[[i]], reduced_coordinates = umap_coords[[i]])
  names(cicero_cds)[i] = names(cds.list)[i]}


# 4. Run cicero and getting pairwise connection between al peaks ----------------------------------------
mratbn7.2 = read.table('./rat.bn7.2.genome.txt')

conns <- lapply(cicero_cds, function(x) run_cicero(x, mratbn7.2, sample_num = 100))

# 5. Filtering connections -----------------------------
conns = lapply(conns, na.omit)
conns = lapply(conns, function(x) filter(x, coaccess > 0.1))

# 6. Annotate peaks -----------------
all.tCREs = read.csv("./all.tCREs.csv") # this table correspond to all peaks getting from runnig <scafe.workflow.cm.aggregate> 
                                        # and annotated as promoter, distal intergenic or intron regions by chipseeker.
all.tCREs$Peak1 = all.tCREs$peak
# Replace all NA values with "Unknown" in the "SYMBOL" column
all.tCREs <- all.tCREs %>% mutate(SYMBOL = ifelse(is.na(SYMBOL), "Unknown", SYMBOL))
# add annotation to Peak in conns
for(i in 1:length(conns)){
  conns[[i]]$Peak1_anno = ifelse(conns[[i]]$Peak1 %in% 
                                   all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
  conns[[i]]$Peak2_anno = ifelse(conns[[i]]$Peak2 %in% 
                                   all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
}

# select enhancer-promoter pairs
conns_6 = lapply(conns, function(x) x[x$Peak1_anno != x$Peak2_anno,])

# 7. Find unique and common connections between cntrl and 5ht
conns_6$CNTRL$in_5ht <- compare_connections(conns_6$CNTRL, conns_6$`5HT`)
conns_6$`5HT`$in_cntrl <- compare_connections(conns_6$`5HT`, conns_6$CNTRL)

# select unique connection (Promoter - Enhancer)
conns_6.unique = lapply(conns_6, function(x) x[x$Peak1_anno == 'Promoter',])

# add gene symbol
conns_6.unique = lapply(conns_6.unique, function(x) merge(x, all.tCREs[,c(26,29)], by="Peak1", all.x=TRUE))

# 8. Make dataframe for hierarhical edge bundling plot --------------------------------------------------

# add column with type of connection
conns_6.unique$CNTRL$connection = ifelse(conns_6.unique$CNTRL$in_5ht == T, 'both','cntrl')
conns_6.unique$`5HT`$connection = ifelse(conns_6.unique$`5HT`$in_cntrl == T, 'both','5ht')
# combine two datasets
all(conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$Peak1 %in% 
      conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$Peak1) # FALSE

all(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$Peak1 %in% 
      conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$Peak1) # FAlSE
# compare connection via pairs of peaks
conns_6.unique$CNTRL$pairs = paste0(conns_6.unique$CNTRL$Peak1, "_", conns_6.unique$CNTRL$Peak2)
conns_6.unique$`5HT`$pairs = paste0(conns_6.unique$`5HT`$Peak1, "_", conns_6.unique$`5HT`$Peak2)

all(conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs %in% 
      conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs) # FALSE

all(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs %in% 
      conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs) # FAlSE

length(conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs) # 11400
length(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs) # 10564

length(conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == F,]$pairs) # 11532
length(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == F,]$pairs) # 4071

sum(conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs %in% 
      conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs) # 9979

sum(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs %in% 
      conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs) # 9979

length(unique(c(conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl == T,]$pairs, 
                conns_6.unique$CNTRL[conns_6.unique$CNTRL$in_5ht == T,]$pairs)))


# combine full CNTRL dataset and unique 5HT interaction
connect = bind_rows(conns_6.unique$CNTRL, conns_6.unique$`5HT`[conns_6.unique$`5HT`$in_cntrl==F,])
connect = connect[,c(1:6,8)]
sum(connect$connection == "cntrl") # 11532
sum(connect$connection == "both")  # 11400
sum(connect$connection == "5ht")   # 4071

length(unique(connect[connect$connection == "cntrl",]$SYMBOL)) # 4700 cntrl unique genes
length(unique(connect[connect$connection == "5ht",]$SYMBOL))   # 2030 5ht unique genes
length(unique(connect[connect$connection == "both",]$SYMBOL))  # 3820 both unique genes

# create a dataframe with connection between leaves (individuals)
head(connect)
connect = connect[,c(2,6,3,7)]
connect$Peak2 = as.character(connect$Peak2)
str(connect)
colnames(connect) = c('from', 'to', 'coaccess', 'connection')
connect = connect[order(connect$connection, connect$from, decreasing=T),] # ordering by 

# 9 Draw edge plot (Figure 5B) --------------------------------------------
# Libraries
library(ggraph)
library(igraph)
library(tidyverse)

# make groups of vertices 
cntrl.genes = unique(connect[connect$connection == 'cntrl',]$to) # have additional connection in cntrl
ht.genes = unique(connect[connect$connection == '5ht',]$to)      # have additional connection in 5ht
cntrl_5ht.genes = cntrl.genes[cntrl.genes %in% ht.genes]             # have different connection between cntrl and 5ht

cntrl.genes = cntrl.genes[!(cntrl.genes %in% cntrl_5ht.genes)]       # remove cntrl_5gt.genes 
ht.genes = ht.genes[!(ht.genes %in% cntrl_5ht.genes)]                # remove cntrl_5gt.genes 

common.genes = unique(connect[connect$connection == 'both',]$to)
common.genes = common.genes[!(common.genes %in% c(cntrl.genes,ht.genes,cntrl_5ht.genes))]  # select only common
sum(common.genes %in% c(cntrl.genes,ht.genes,cntrl_5ht.genes))

enhancers = as.character(unique(connect$from)) # 11307

#  create a data frame giving the hierarchical structure of your individuals
d1 <- data.frame(from="Origin", to=c('Enhancer', 'both', 'cntrl', '5ht', 'cntrl_5ht'))
d2 <- data.frame(from=c(rep("Enhancer", length(enhancers)),
                        rep("both", length(common.genes)), 
                        rep("cntrl", length(cntrl.genes)),
                        rep("5ht", length(ht.genes)),
                        rep("cntrl_5ht", length(cntrl_5ht.genes))),
                 to=c(enhancers,
                      common.genes,
                      cntrl.genes,
                      ht.genes,
                      cntrl_5ht.genes)
)
hierarchy <- rbind(d1, d2)

# create a vertices data.frame. One line per object of our hierarchy
vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))) 
vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]


# Create a graph object
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)
col = match(connect$connection, vertices$name)

# Advanced graph
library(RColorBrewer)
# nodes color 5ht -> both -> cntrl -> cntrl_5ht -> enh 

p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + theme_void()

p + geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group), size=1) +
  #scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)[6:10]) +
  scale_colour_manual(values= c('#fb7669', "grey", '#56B1F7',"#33A02C", "#e10c72" )) +
  geom_conn_bundle(data = get_con(from = from, to = to,col = connect$connection), alpha=0.1, aes(color=col), tension = .9)

# cntrl
sum(connect$connection == "cntrl") #11532
p + geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group), size=1) +
  scale_colour_manual(values= c('#fb7669', "grey", '#56B1F7',"#33A02C", "#e10c72" ),
                      name = "Nodes", labels = c("5HT a.c.", "CNTRL a.c.", "Common c.", "5HT-CNTRL a.c.", "Enhancers")) +
  geom_conn_bundle(data = get_con(from = from[1:11532], to = to[1:11532],
                                  col = connect$connection[1:11532]), alpha=0.05, colour='#56B1F7', tension = .9)
# both
sum(connect$connection == "both") #11400
p + geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group), size=1) +
  scale_colour_manual(values= c('#fb7669', "grey", '#56B1F7',"#33A02C", "#e10c72" ),
                      name = "Nodes", labels = c("5HT a.c.", "CNTRL a.c.", "Common c.", "5HT-CNTRL a.c.", "Enhancers")) +
  geom_conn_bundle(data = get_con(from = from[11533:22932], to = to[11533:22932],
                                  col = connect$connection[11533:22932]), alpha=0.05,  colour='lightgrey', tension = .9)
# 5ht
sum(connect$connection == "5ht") # 4071
p + geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group), size=1) +
  scale_colour_manual(values= c('#fb7669', "grey", '#56B1F7',"#33A02C", "#e10c72" ),
                      name = "Nodes", labels = c("5HT a.c.", "CNTRL a.c.", "Common c.", "5HT-CNTRL a.c.", "Enhancers")) +
  geom_conn_bundle(data = get_con(from = from[22933:27004], to = to[22933:27004],
                                  col = connect$connection[22933:27004]), alpha=0.05,  colour='#fb7669', tension = .9)

# ...8.1 Visualization archplots (Figure 6B) ------------------------------------------------------------------
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

"Nfia" %in% cntrl.genes
"Nfix" %in% cntrl_5ht.genes
"Nfib" %in% cntrl_5ht.genes

conns_6.unique$CNTRL$conns_color = "#56B1F7"
conns_6.unique$CNTRL$conns_color[conns_6.unique$CNTRL$in_5ht] = "lightgray"

conns_6.unique$`5HT`$conns_color = "#fb7669"
conns_6.unique$`5HT`$conns_color[conns_6.unique$`5HT`$in_cntrl] = "lightgray"

# NFIX
plot_connections(conns_6.unique$CNTRL, 
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
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 1.5, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr19", 23355498-200000, 23448265+100000,
                 viewpoint = "chr19_23354187_23356823",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 1.5, 
                 collapseTranscripts = "longest")

# NFIA
plot_connections(conns_6.unique$CNTRL, 
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
                 connection_width = 1.5, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr5", 112439724-230000, 112443239+300000,
                 viewpoint = "chr5_112439724_112443239",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = 1.5,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 1.5, 
                 collapseTranscripts = "longest")

# NFIB
plot_connections(conns_6.unique$CNTRL, 
                 "chr5", 96970365-270000, 96976509+230000,
                 #viewpoint = "chr5_96970365_96976509",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 comparison_peak_color = "#fb7669",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 1.5, 
                 collapseTranscripts = "longest")+NoAxes()

plot_connections(conns_6.unique$`5HT`, 
                 "chr5", 96970365-270000, 96976509+230000,
                 #viewpoint = "chr5_96970365_96976509",
                 connection_color = "conns_color",
                 comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.0,
                 comparison_coaccess_cutoff = 0.0,
                 connection_width = 1.5*2, 
                 collapseTranscripts = "longest")

# supplementary figure 20 "Explanation of edge plots" ------------------------------------------
conns_6.unique = lapply(conns_6.unique, drop_na)
sapply(conns_6.unique, function(x) sum(is.na(x$SYMBOL)))

# Common connections
common.genes
conns_6.unique$CNTRL[conns_6.unique$CNTRL$SYMBOL == "Slitrk1",]
conns_6.unique$`5HT`[conns_6.unique$`5HT`$SYMBOL == "Slitrk1",]

plot_connections(conns_6.unique$CNTRL, 
                 "chr15", 85688283-10000, 85728537+10000,
                 viewpoint = "chr15_85725319_85728537",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 #comparison_peak_color = "#fb7669",
                 #comparison_connection_color = "lightgray",
                 #comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 #comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr15", 85688283-10000, 85728537+10000,
                 viewpoint = "chr15_85725319_85728537",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")


# Ctrl - specific connections
cntrl.genes
conns_6.unique$CNTRL[conns_6.unique$CNTRL$SYMBOL == "H1f6",]
conns_6.unique$`5HT`[conns_6.unique$`5HT`$SYMBOL == "H1f6",]

plot_connections(conns_6.unique$CNTRL, 
                 "chr17", 41391587-10000, 41575970+10000,
                 viewpoint = "chr17_41428008_41428509",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 #comparison_peak_color = "#fb7669",
                 #comparison_connection_color = "lightgray",
                 #comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 #comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr17", 41391587-10000, 41575970+10000,
                 viewpoint = "chr17_41428008_41428509",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

# 5HTP - specific connections
ht.genes
conns_6.unique$CNTRL[conns_6.unique$CNTRL$SYMBOL == "Ctnna1",]
conns_6.unique$`5HT`[conns_6.unique$`5HT`$SYMBOL == "Ctnna1",]

plot_connections(conns_6.unique$CNTRL, 
                 "chr18", 26381039-10000, 26729484+50000,
                 viewpoint = "chr18_26727561_26729733",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 #comparison_peak_color = "#fb7669",
                 #comparison_connection_color = "lightgray",
                 #comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 #comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr18", 26381039-10000, 26729484+50000,
                 viewpoint = "chr18_26727561_26729733",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")


# Dual - specific connections
cntrl_5ht.genes
conns_6.unique$CNTRL[conns_6.unique$CNTRL$SYMBOL == "Snhg11",]
conns_6.unique$`5HT`[conns_6.unique$`5HT`$SYMBOL == "Snhg11",]

plot_connections(conns_6.unique$CNTRL, 
                 "chr3", 146964892-10000, 147224895+10000,
                 viewpoint = "chr3_147030413_147031777",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$`5HT`,
                 peak_color = "#56B1F7",
                 #comparison_peak_color = "#fb7669",
                 #comparison_connection_color = "lightgray",
                 #comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 #comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")

plot_connections(conns_6.unique$`5HT`, 
                 "chr3", 146964892-10000, 147224895+10000,
                 viewpoint = "chr3_147030413_147031777",
                 connection_color = "conns_color",
                 #comparison_track = conns_6.unique$CNTRL,
                 peak_color = "#fb7669",
                 comparison_peak_color = "#56B1F7",
                 comparison_connection_color = "lightgray",
                 comparison_connection_width = .75,
                 gene_model_color = "#2DD881",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.1,
                 comparison_coaccess_cutoff = 0.1,
                 connection_width = 3, 
                 collapseTranscripts = "longest")




# analysis of co-transcribed regulatory modules activity among general cell types

### 1 Create full cds object  ---------- 
cds.full = as.cell_data_set(e20.scafe)

# since it misses the gene_short_name column, let's add it 
fData(cds.full)$gene_short_name = rownames(cds.full)

# ..1.2. Cluster cells (using clustering info from seurat's UMAP)
# let's use the clustering information have
recreate.partition = c(rep(1,length(cds.full@colData@rownames)))
names(recreate.partition) = cds.full@colData@rownames
recreate.partition = as.factor(recreate.partition)
cds.full@clusters$UMAP$partitions <- recreate.partition

# Assign the cluster info 
Idents(e20.scafe) = e20.scafe@meta.data$general_celltype
clusters = e20.scafe@active.ident
cds.full@clusters$UMAP$clusters <- clusters

# Assign UMAP coordinate - cell embeddings
cds.full@int_colData@listData$reducedDims$UMAP  = e20.scafe@reductions$umap@cell.embeddings

# ..1.3. make cicero cds. 
umap_coords = reducedDims(cds.full)$UMAP
cicero_cds.full = make_cicero_cds(cds.full, reduced_coordinates = umap_coords)

# ..1.4 run cicero
conns_full <- run_cicero(cicero_cds.full, genomic_coords = mratbn7.2, sample_num = 100)

conns_full$Peak1_anno = ifelse(conns_full$Peak1 %in% all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
conns_full$Peak2_anno = ifelse(conns_full$Peak2 %in% all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')

conns_full = conns_full[conns_full$Peak1_anno != conns_full$Peak2_anno,]

ccans <- generate_ccans(conns_full, coaccess_cutoff_override = 0.1)
links = ConnectionsToLinks(conns = conns_full, ccans = ccans)

e20.signac = readRDS("./e20.signac.Rds")

Links(e20.signac) = links

# 3.1 set names of subsetted object
Idents(e20.scafe) = e20.scafe@meta.data$stim_gct
stim_gct = sort(levels(Idents(e20.scafe)), decreasing = T)

serotonin.stim_idents.list = list()
for(st_gct in stim_gct){
  print(st_gct)
  v = subset(e20.scafe, subset = stim_gct == st_gct)
  print(dim(v))
  serotonin.stim_idents.list  = base::append(serotonin.stim_idents.list, v)}
names(serotonin.stim_idents.list) = stim_gct
serotonin.stim_idents.list

all(sort(sapply(serotonin.stim_idents.list, function(x) dim(x)[2])) ==
      sort(table(e20.scafe@meta.data$stim_gct)))

# upload cds_full and add new cluster metadata
e20.meta = read.csv("~/serotonin6/srat-genes/e20_meta.csv")
e20.meta = e20.meta[e20.meta$cell %in% colnames(cds),]
e20.meta <- e20.meta[match(colnames(cds), e20.meta$cell), ]
all(e20.meta$cell == colnames(cds.full))


cds.full@colData["gct"] = e20.meta$cell_types
cds.full@colData["stim_gct"] = paste0(cds.full@colData$stim, "_", cds.full@colData$gct)

# remove DG and Irx cells
# Фильтруем клетки, исключая те, где значения равны "DG" или "Irx1.Irx2.Lh9" в колонке "cell_type"
cds.full <- cds.full[, !cds.full@colData$gct %in% c("DG", "Irx1.Irx2.Lh9")]
cds

# create list of cell_data_set object
stim_gct = unique((cds.full@colData$stim_gct))

cds.list = list()
for(st_gct in stim_gct){
  print(st_gct)
  # Выбираем клетки с аннотацией "CNTRL_IN" в колонке stim_gct
  v <- cds.full[, cds.full@colData$stim_gct == st_gct]
  
  print(dim(v))
  cds.list[[st_gct]] <- v
}
cds.list

names(cds.list) = c("CNTRL_IN", "CNTRL_EX", "CNTRL_BRIDGE_1", "CNTRL_GLIA", "CNTRL_BRIDGE_2",
                    "5HT_EX", "5HT_GLIA", "5HT_IN", "5HT_BRIDGE_1", "5HT_BRIDGE_2")
cds.list = cds.list[order(names(cds.list), decreasing = T)]

# 3.2 create cicero_cds
# optimaze k for kNN by k = √N/2
k0.5.values = round((sqrt(sapply(cds.list, function(x) dim(x)[2]))/2))
#
cicero_cds = list()

for(i in 1:length(cds.list)){
  cicero_cds[[i]] =  make_cicero_cds(cds.list[[i]], reduced_coordinates = reducedDims(cds.list[[i]])$UMAP, k = k0.5.values[i])
  names(cicero_cds)[i] = names(cds.list)[i]}


# 3.3 run cicero
# estimate distance parameter
paramaters_cicero_cds = lapply(cicero_cds, function(x) estimate_distance_parameter(cds = x, genomic_coords = mratbn7.2))
paramaters_cicero_cds_mean = sapply(paramaters_cicero_cds, mean)

models = parallel::mcmapply(function(x, y) 
  generate_cicero_models(cds = x, 
                         distance_parameter = y, 
                         genomic_coords = mratbn7.2), 
  cicero_cds, 
  paramaters_cicero_cds_mean, 
  mc.cores = 12, SIMPLIFY = FALSE)

conns.st_gct = mclapply(models, function(x) assemble_connections(x,),mc.cores = 12)


for(i in 1:length(conns.st_gct)){
  conns.st_gct[[i]]$Peak1_anno = ifelse(conns.st_gct[[i]]$Peak1 %in% 
                                          all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
  conns.st_gct[[i]]$Peak2_anno = ifelse(conns.st_gct[[i]]$Peak2 %in% 
                                          all.tCREs[all.tCREs$annotation == "Promoter", ]$peak,'Promoter', 'Enhancer')
}

conns.st_gct = lapply(conns.st_gct, function(x) x[x$Peak1_anno != x$Peak2_anno,])

# split e20.signac by st_gct 
e20.meta_signac = e20.meta[e20.meta$cell %in% colnames(e20.signac),]
e20.meta_signac <- e20.meta[match(colnames(e20.signac), e20.meta$cell), ]

e20.signac@meta.data["cell_type_new"] = e20.meta_signac$cell_types
e20.signac[["stim_gct"]] = paste0(e20.signac$stim, "_", e20.signac$cell_type_new)

# umap logo 
DimPlot(
  e20.signac, 
  group.by = "cell_type_new", 
  cols = c("#dd9500", "#0077d5", "#b5b5b5", "#ff6450", "#009f0b", "#9f78f3", "#b5b5b5", "#42cc7c"),
  pt.size = 1
) + 
  NoAxes() + 
  labs(title = NULL)  

Idents(e20.signac) = e20.signac@meta.data$stim_gct
levels(e20.signac)
stim_gct = levels(Idents(e20.signac))

signac.stim_idents.list = list()
for(st_gct in stim_gct){
  print(st_gct)
  v = subset(e20.signac, subset = stim_gct == st_gct)
  print(dim(v))
  signac.stim_idents.list  = base::append(signac.stim_idents.list, v)}
names(signac.stim_idents.list) = stim_gct
signac.stim_idents.list

all(sort(sapply(signac.stim_idents.list, function(x) dim(x)[2])) ==
      sort(table(e20.signac@meta.data$stim_gct)))

names(signac.stim_idents.list)
signac.stim_idents.list = signac.stim_idents.list[c(1,2,3,4,6,7,9,10,11,12,13,14)]
names(signac.stim_idents.list)
# add conns ccans to signac subsets
names(conns.st_gct)

links_list = lapply(conns.st_gct, function(x) ConnectionsToLinks(conns = x, ccans = ccans))


all(names(links_list) == names(signac.stim_idents.list))
all(names(links_list) == names(conns.st_gct))
all(names(signac.stim_idents.list) == names(conns.st_gct))


for(i in 1:length(links_list)) {
  Links(signac.stim_idents.list[[i]]) = links_list[[i]]
}


# count co-activity score among new general cell types beetwen cntrl and 5ht
names(signac.stim_idents.list)
signac.stim_idents.list = signac.stim_idents.list[c(1,6,5,4,3,2,9,10,12,8,11,7)]
names(signac.stim_idents.list)
signac.stim_idents.list = signac.stim_idents.list[c(1:7, 9,8, 10:12)]
signac.control = signac.stim_idents.list[1:6]
signac.ht = signac.stim_idents.list[7:12]
general_celltypes = c("Inh. neu.", "Ascl1-axis", "Progen.", "Glia", "Neurog-axis", "Ex. neu.")

Links(signac.stim_idents.list$CNTRL_Inhibitory_Neurons)
Links(serotonin.signac)
library(ggpubr)
library(reshape2)
num = 1
df_new = data.frame()

# Formula for counting co-activity score via log weighted mean
weighted_mean <- function(mean_value, sample_size) {
  return(mean_value * log(sample_size + 1))
}

for(gct in 1:6){
  df <- data.frame()
  a <- data.frame(Links(signac.control[[gct]]))
  b <- data.frame(Links(signac.ht[[gct]]))
  print(length(intersect(unique(a$group),unique(b$group)))/length(union(unique(a$group),unique(b$group))))
  for (i in unique(ccans$CCAN)) {
    if (i %in% unique(a$group)) {
      temp_1 <- weighted_mean(mean(a[a$group %in% i,]$score), length(a[a$group %in% i,]$score))
    } else {
      temp_1 <- 0
    }
    if (i %in% unique(b$group)) {
      temp_2 <- weighted_mean(mean(b[b$group %in% i,]$score), length(b[b$group %in% i,]$score))
    } else {
      temp_2 <- 0
    }
    df <- rbind(df,data.frame(group=i,control=temp_1,ht=temp_2))
  }
  #  rownames(df) <- df$group
  #  df <- df[,-1]
  
  my_comparison <- list(c("control","ht"))
  p <- ggpubr::ggpaired(df, cond1 = "control", cond2 = "ht",
                        fill = "condition", palette = "jco",line.size = 0.1,point.size = 0.2)+
    stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "wilcox.test")
  assign(paste0("p",num),p) 
  df_new <- rbind(df_new,
                  data.frame(melt(df,measure.vars = c("control", "ht")),cluster=general_celltypes[gct]))
  num <- num + 1
}

p1 | p2 | p3 | p4 | p5 | p6

# final var
df_new$log.value = -log10(df_new$value)

# Figure 5C
ggboxplot(df_new, x = "cluster", y = "value", fill = "variable", size=NULL, outlier.shape = 20)  +
  stat_compare_means(aes(group=variable), label = "p.signif", method = "wilcox.test", vjust = -1.5) +
  labs(y = "co-activity score", x = "") +
  scale_fill_manual(values = c("#56B1F7", "#fb7669"), labels = c("Ctrl", "5HT"),name = "") +
  theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) +
  theme(panel.border = element_blank()) +
  theme_set(theme_minimal()) +
  theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) +
  theme(panel.border = element_blank()) +
  theme(axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain", size = 14,colour = "black", angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(face = "plain", size = 14, colour = "black"))

p.adjust(c(0.0440, 1.6e-06, 0.0076, 0.0122, 4.7e-11, 0.0416), method = "BH")

# Analysis of motif activity using ChroVar ------------------
# NFI/CTF family

View(e20.signac@meta.data)
e20.signac$stim[e20.signac$stim == "CNTRL"] = "Ctrl"
e20.signac$stim[e20.signac$stim == "5HT"] = "5HTP"

e20.signac$stim = factor(e20.signac$stim, levels = c("Ctrl", "5HTP"))

# Figure 6C left (UMAP)
FeaturePlot(
  object = e20.signac,
  features = c("MA0670.1", "MA1643.1", "MA0671.1"),
  cols = SpatialColors(n = 100)[51:100],  
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  split.by = "stim",
  label = FALSE
) & 
  theme_void() & 
  theme(
    legend.position = "none",      
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  
  ) 


MotifPlot(object = e20.signac, motifs = c("MA0670.1"),assay = 'peaks') & theme_void()
MotifPlot(object = e20.signac, motifs = c("MA1643.1"),assay = 'peaks') & theme_void()
MotifPlot(object = e20.signac, motifs = c("MA0671.1"),assay = 'peaks') & theme_void()

levels(Idents(e20.signac))

# Violin plots for Figure 6C


signac.subset = subset(e20.signac, subset = cluster_type %in% c("Bridge1_Neurod1.Neurod2", "Tanycytes"))
signac.subset
DefaultAssay(signac.subset)
signac.subset$cluster_type[signac.subset$cluster_type == "Bridge1_Neurod1.Neurod2"] = "Neurog-axis"
Idents(signac.subset) = signac.subset$cluster_type
# NFIA
plot_data <- VlnPlot(signac.subset, 
                     assay = NULL, 
                     features = "MA0670.1", #NFIA
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data

plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP",
                                                      "Neurog-axis.Ctrl", "Neurog-axis.5HTP"))
# Violin plot + boxplots
library(ggpubr)

nfia_chrom = ggplot(plot_data, aes(x = group, y = MA0670.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  scale_y_continuous(breaks = c(-2, 0, 2)) + 
  labs(title = "", y = "", x = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size=16)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP"),
                                        c("Neurog-axis.Ctrl", "Neurog-axis.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  

# NFIB

plot_data <- VlnPlot(signac.subset, 
                     assay = NULL, 
                     features = "MA1643.1", #NFIB
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data
plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP",
                                                      "Neurog-axis.Ctrl", "Neurog-axis.5HTP"))

# Violin plot + boxplots
nfib_chrom = ggplot(plot_data, aes(x = group, y = MA1643.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +  
  labs(title = "", y = "ChromVar score", x = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size=16)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP"),
                                        c("Neurog-axis.Ctrl", "Neurog-axis.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  

# NFIX
plot_data <- VlnPlot(signac.subset, 
                     assay = NULL, 
                     features = "MA0671.1",  # NFIX
                     split.by = "stim",
                     cols = c("#56B1F7", "#fb7669"),
                     pt.size = 0.1, 
                     combine = FALSE)[[1]]$data
plot_data$group <- interaction(plot_data$ident, plot_data$split)
levels(plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Tanycytes.Ctrl", "Tanycytes.5HTP",
                                                      "Neurog-axis.Ctrl", "Neurog-axis.5HTP"))
# Violin plot + boxplots
nfix_chrom = ggplot(plot_data, aes(x = group, y = MA0671.1, fill = split)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +  
  labs(title = "", y = "", x = "") +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),  
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 0),
        axis.text.x = element_text(size=0, angle = 45, vjust = 1, hjust = 1)) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("Tanycytes.Ctrl", "Tanycytes.5HTP"),
                                        c("Neurog-axis.Ctrl", "Neurog-axis.5HTP")), 
                     method = "wilcox.test", 
                     label = "p.signif", bracket.size = 1, size=0)  

# Find differentially active motifs between stim for Neurog-axis cells -----------
DefaultAssay(e20.signac) <- 'chromvar'
table(e20.signac$stim_gct)
# bridge1 DAMs cntrl 5ht
bridge1.dams <- FindMarkers(
  object = e20.signac,
  ident.1 = '5HT_Bridge1_Neurod1.Neurod2',
  ident.2 = 'CNTRL_Bridge1_Neurod1.Neurod2',
  test.use = "LR",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  latent.vars = "nCount_peaks"
)

bridge1.dams.motifs = ConvertMotifID(object = serotonin.signac,
                                     id = rownames(bridge1.dams),
                                     assay = 'peaks')
bridge1.dams$motif = bridge1.dams.motifs

bridge1.dams.s = bridge1.dams[bridge1.dams$p_val_adj<0.001 & abs(bridge1.dams$avg_diff) > 0.4,]

# TFs motifs that regulate RXF genes
upstream.rfx.regulators = c("SP2","E2F4","KLF16","SP8","SP3","EGR3","ESR1","Creb5","ZNF740","MZF1(var.2)","ATF7","SOX21","MZF1","Tcfl5","KLF5","SP1","EGR1","TFAP2C(var.3)","JDP2(var.2))")
# cntrl specific motifs in RFX regulators
bridge1.dams.s %>% dplyr::filter(motif %in% upstream.rfx.regulators & avg_diff < 0) %>% pull(motif)
# 5ht specific motifs in RFX regulators
bridge1.dams.s %>% dplyr::filter(motif %in% upstream.rfx.regulators & avg_diff > 0) %>% pull(motif)

# Figure 6D
ggplot(data = bridge1.dams, aes(x = avg_diff, y = -log10(p_val_adj))) +
  geom_point(aes(color = 
                   ifelse(avg_diff < 0.4 & p_val_adj > 0.001, "gray",
                          ifelse(avg_diff >= 0.4 & p_val_adj < 0.001, "#fb7669",
                                 ifelse(avg_diff <= -0.4 & p_val_adj < 0.001, "#56B1F7", "gray")))),
             alpha = 1, 
             size = 2.5) +  # Прозрачность и размер точек
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "brown") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "brown") +
  geom_text(data = subset(bridge1.dams, p_val_adj < 0.001 & abs(avg_diff) > 0.4),
            aes(label = motif), 
            size = 4, 
            hjust = ifelse(bridge1.dams[bridge1.dams$p_val_adj < 0.001 
                                        & abs(bridge1.dams$avg_diff) > 0.4,]$avg_diff < 0, 1.2, -0.2),
            check_overlap = T,
            color = "black") +  # Подписи только для значимых мотивов
  xlim(-1.5, 1.5) +
  scale_color_identity() +  # Цвет точек
  theme_minimal() +  # Минималистичная тема
  labs(title = "Neurod-axis DAMs", x = "avg_diff", y = "-log10(p_val_adj)")


# Find differentially active motifs between stim for Ascl1-axis cells -----------

bridge2.dams <- FindMarkers(
  object = e20.signac,
  ident.1 = '5HT_Bridge2_Ascl1',
  ident.2 = 'CNTRL_Bridge2_Ascl1',
  test.use = "LR",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  latent.vars = "nCount_peaks"
)

bridge2.dams.motifs = ConvertMotifID(object = e20.signac,
                                     id = rownames(bridge2.dams),
                                     assay = 'peaks')
bridge2.dams$motif = bridge2.dams.motifs

bridge2.dams.s = bridge2.dams[bridge2.dams$p_val_adj<0.001 & abs(bridge2.dams$avg_diff) > 0.4,]

# Volcano plot (Figure 6D)
ggplot(data = bridge2.dams, aes(x = avg_diff, y = -log10(p_val_adj))) +
  geom_point(aes(color = 
                   ifelse(avg_diff < 0.4 & p_val_adj > 0.001, "gray",
                          ifelse(avg_diff >= 0.4 & p_val_adj < 0.001, "red",
                                 ifelse(avg_diff <= -0.4 & p_val_adj < 0.001, "blue", "gray")))),
             alpha = 0.5, 
             size = 2) +  
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "brown") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "brown") +
  geom_text(data = subset(bridge2.dams, p_val_adj < 0.001 & abs(avg_diff) > 0.4),
            aes(label = motif), 
            size = 2, 
            hjust = ifelse(bridge2.dams[bridge2.dams$p_val_adj < 0.001 
                                        & abs(bridge2.dams$avg_diff) > 0.4,]$avg_diff < 0, 1.2, -0.2),
            check_overlap = T,
            color = "black") +  
  xlim(-1.5, 1.5) +
  scale_color_identity() +  
  theme_minimal() + 
  labs(title = "Volcano Plot", x = "avg_diff", y = "-log10(p_val_adj)")


# DOT PLOT FOR INTERSTING MOTIFS --------------------------------------------
# add column with general cell types
ct_before = serotonin.signac@meta.data$cluster_type

rules <- data.frame(
  celltype = c("GABA_Meis2.Bcl11b",            "Glut_Peg10.Calcr",          "Bridge1_Neurod1.Neurod2", "GABA_Dlx1",
               "Glut_Lhx9",                    "GABA_Asic4.Gata3.Tal1",     "OPC",                     "GABA_Bcl11b.Arx.Lhx6",    
               "Tanycytes",                    "Irx1.Irx2.Lh9",             "Astrocytes",              "GABA_Bcl11b.Adarb.Lhx6", 
               "GABA_Meis2.Isl1.Colq.Slc44a5", "Progenitors",               "Glut_Tfap2d",             "Bridge2_Ascl1",
               "Glut_Tbr1.Meis2.Ebf3",         "Ependyma",                  "Glut_Lmx1a.Ebf3",         "DG",
               "GABA_Peg10",                   "Glut_Foxb1.Sim1.Nxph4",     "Glut_Nr5a1.Pdyn.Fezf1",   "GABA_Sim1"),
  
  general_celltype = c("Inh. neu.",            "Ex. neu.",                  "Neurod-axis",             "Inh. neu.", 
                       "Ex. neu.",             "Inh. neu.",                 "OPC",                     "Inh. neu.", 
                       "Tanycytes",            "Irx1.Irx2.Lh9",             "Astrocytes",              "Inh. neu.", 
                       "Inh. neu.",            "Progenitors",               "Ex. neu.",                "Ascl1-axis", 
                       "Ex. neu.",             "Ependyma",                  "Ex. neu.",                "DG", 
                       "Inh. neu.",            "Ex. neu.",                  "Ex. neu.",                 "Inh. neu.")
)

e20.signac@meta.data[["celltypes2"]] <- 
  rules$general_celltype[match(e20.signac@meta.data$cluster_type,rules$celltype)]

ct_after = e20.signac@meta.data$celltypes2
all(ct_before == ct_after) # TRUE
unique_values <- unique(subset(e20.signac@meta.data, select=c(cluster_type, celltypes2)))
Idents(e20.signac) = e20.signac$celltypes2

levels(Idents(serotonin.signac))[c(1:5, 7:10)]

serotonin.signac.subset <- subset(
  e20.signac,
  subset = (celltypes2 %in% levels(Idents(e20.signac))[c(1:5, 7:10)])
)

levels(Idents(serotonin.signac.subset))

Idents(serotonin.signac.subset) = serotonin.signac.subset$celltypes2
Idents(serotonin.signac.subset) <- factor(Idents(serotonin.signac.subset), 
                                          levels = c("Inh. neu.", "Ascl1-axis", "Progenitors", "OPC", 
                                                     "Astrocytes", "Ependyma", "Tanycytes", "Neurod-axis", "Ex. neu."))
ko_enrich = c("MA0057.1", "MA0660.1", "MA0768.1", "MA0700.2", "MA0798.2", "MA0509.2")

View(serotonin.signac.subset@meta.data)
serotonin.signac.subset$stim[serotonin.signac.subset$stim == "CNTRL"] <- "Ctrl"
serotonin.signac.subset$stim[serotonin.signac.subset$stim == "5HT"] <- "5HTP"

serotonin.signac.subset@meta.data[["celltypes2_stim"]] = paste0(serotonin.signac.subset$celltypes2, "_", 
                                                                serotonin.signac.subset$stim)
Idents(serotonin.signac.subset) = serotonin.signac.subset$celltypes2_stim
levels(Idents(serotonin.signac.subset))
Idents(serotonin.signac.subset) <- factor(Idents(serotonin.signac.subset), 
                                          levels = c("Inh. neu._Ctrl", "Inh. neu._5HTP",
                                                     "Ascl1-axis_Ctrl", "Ascl1-axis_5HTP",
                                                     "Progenitors_Ctrl", "Progenitors_5HTP",
                                                     "OPC_Ctrl", "OPC_5HTP",
                                                     "Astrocytes_Ctrl", "Astrocytes_5HTP",
                                                     "Ependyma_Ctrl", "Ependyma_5HTP",
                                                     "Tanycytes_Ctrl",  "Tanycytes_5HTP",
                                                     "Neurod-axis_Ctrl", "Neurod-axis_5HTP",
                                                     "Ex. neu._Ctrl", "Ex. neu._5HTP" ))


dp = DotPlot(serotonin.signac.subset, features = ko_enrich, 
             cols = c("blue", "white", "red"), 
             dot.scale = 10) + 
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200))) + RotatedAxis() + coord_flip()

dpdf = dp$data

rownames(dpdf) = 1:nrow(dpdf)
str(dpdf)
dpdf$id  = as.character(dpdf$id)
dpdf$celltype = sapply(dpdf$id, function(x) strsplit(x, "_")[[1]][1])
dpdf$stim = sapply(dpdf$id, function(x) strsplit(x, "_")[[1]][2])

gene_order <- rev(c("MA0057.1", "MA0660.1", "MA0768.1", "MA0700.2", "MA0798.2", "MA0509.2"))
celltype_order <- c("Inh. neu.", "Ascl1-axis", "Progenitors", "OPC", 
                    "Astrocytes", "Ependyma", "Tanycytes", "Neurod-axis", "Ex. neu.")

dpdf <- dpdf %>%
  mutate(
    features.plot = factor(features.plot, levels = gene_order),  
    celltype = factor(celltype, levels = celltype_order)         
  ) %>%
  arrange(features.plot, stim) %>%  
  mutate(
    features.stim = factor(interaction(features.plot, stim, sep = "_"), 
                           levels = unique(interaction(features.plot, stim, sep = "_"))),
    

    features_label = rep(c("5HTP", "Ctrl"), length.out = nrow(dpdf))
  )


#  DotPlot (Figure 6K)
ggplot(dpdf, aes(x = celltype, y = features.stim, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() + 

  scale_size(name = "Cell expression ratio", range = c(3, 10)) +  

  scale_color_gradientn(name = "Chromvar Z-score", colors = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200))) +

  theme_classic() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=16),  
    axis.text.y = element_text(size=16), 
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),  
    legend.title = element_text(angle = 90, vjust = 1.5)
  ) +
  scale_y_discrete(labels = dpdf$features_label)



