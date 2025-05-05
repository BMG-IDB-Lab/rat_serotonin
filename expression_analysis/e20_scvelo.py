import scvelo as scv
import loompy
import scanpy as sc
import scipy
import numpy as np
import pandas as pd
import re
import os

# settings
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo', figsize=[5,4])  # for beautified visualization

# combine looms for all samples
loom = os.listdir('results_loom')
loom = list(filter(lambda x: 'merged' not in x, loom))
loom = [os.path.join('results_loom', file) for file in loom]
out = 'results_loom/merged.loom'
loompy.combine(loom, out, key="Accession")

# data
data_dir = "data"
loom_dir = "results_loom"
# loom
output_filename = f"{loom_dir}/merged.loom"
adata_loom = scv.read(output_filename)
# counts and metadata from seurat
adata_s = sc.read_h5ad(f"{data_dir}/whole.integrated_5HT_E20_int.h5ad")
# raw to X
tempAdata = adata_s.raw.to_adata() 
tempAdata.var_names = adata_s.var_names 
adata_s.raw = tempAdata 
adata_s = adata_s.raw.to_adata()
# rename cells for consistency
sample_to_number = {
    'Rat_1_CR7default' : '1',
    'Rat_2_CR7_fc24k_new' : '2',
    'Rat_3_CR7_fc85k_new' : '3',
    'Rat_4_CR7_fc12k_new' : '4',
    'Rat_5_CR7_fc105k_new' : '5',
    'Rat_6_CR7default' : '6'}
adata_loom.obs_names = list(map(lambda x: re.sub("x", "-1" , re.sub(".*:", "", x)) + "_" + sample_to_number[re.sub(":.*", "", x)], 
         adata_loom.obs_names))
# use umap.orig instead of umap
adata_s.obsm["X_umap"] = adata_s.obsm["X_umap.orig"]
# assign colors to clusters
adata_s.uns["cell_types_colors"] = [
    "#DF9200", "#006FD2", "#B3B3B3",
    "#FF6347", "#019B0B", "#9F79EE",
    "#CDB38B", "#43CD80"
]
# merge loom and adata
adata = scv.utils.merge(adata_s, adata_loom)

# preprocessing
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=3000, enforce=True)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs=15, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# scvelo statistics
scv.pl.proportions(adata, groupby="cell_types")
# run with stochastic model
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata, n_jobs = 32)
# plot results
scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell_types', 
                                legend_loc='none', size=20, alpha=0.5,
                                save="stream_alpha.0.5.png",
                                dpi = 600)
                                
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', colorbar=False)
