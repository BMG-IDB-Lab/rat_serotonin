import seaborn as sns
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import numba
import pynndescent
import matplotlib.patches as mpatches
import glob


#Ploting one slice of MERFISH data
def plot_section_one(xx, yy, cc = None, val = None, fig_width = 10, fig_height = 8, cmap = None):
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)
    if cmap is not None:
        plt.scatter(xx, yy, s=15, c=val, marker='.', cmap=cmap)
    elif cc is not None:
        scatter=plt.scatter(xx, yy, s=15, color=cc)
    ax.set_ylim(11, 0)
    ax.set_xlim(0, 11)
    ax.axis('equal')
    plt.locator_params(axis='both', nbins=50)
    plt.grid()
    return fig, ax

#Ploting all slices of MERFISH data
def ploting_all_sections(meta_plot, label_column, color_column, name='no_name'):
    secs=[]
    for i in np.unique(meta_plot['brain_section_label']):
        if ('HY' in meta_plot[meta_plot['brain_section_label'] == i].division.values):
            secs.append(i)
    
    
    for i in range(len(secs)):
        ntexp_sec = meta_plot[meta_plot['brain_section_label'] == secs[i]]
        fig, ax,  = plot_section_one(ntexp_sec['x'], ntexp_sec['y'], cc=ntexp_sec[color_column])
        fig.set_size_inches(15,15)
        df1 = ntexp_sec
        a=df1.groupby([label_column,color_column]).size().reset_index()
        names = a[label_column]
        colors = a[color_column]
        recs = []
    
        for j in range(0, len(colors)):
            recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=colors[j]))
        plt.legend(recs, names, bbox_to_anchor= [1,1.02], ncol=1)
        plt.title(secs[i])
        plt.locator_params(axis='both', nbins=50)
        plt.grid()
        plt.savefig(f'/home/z.starinnov/spatial/slides2/{secs[i]}.png', pad_inches=0.2)
    from PIL import Image
    if name!='no_name':
        images = [
        Image.open(f"/home/z.starinnov/spatial/slides2/{f}.png") for f in secs]
        
        pdf_path = f"/home/z.starinnov/spatial/pdfs/{name}.pdf"
            
        images[0].save(
            pdf_path, "PDF" ,resolution=100.0, save_all=True, append_images=images[1:]
        )
    files = glob.glob('/home/z.starinnov/spatial/slides2/*')
    for f in files:
        os.remove(f)
#Plotting all slices of MERFISH data with levels of expression of certain genes
def ploting_all_sections_exp(meta_plot, gene_column, name='no_name'):
    secs=[]
    for i in np.unique(meta_plot['brain_section_label']):
        if ('HY' in meta_plot[meta_plot['brain_section_label'] == i].division.values):
            secs.append(i)
    for i in range(len(secs)):
        ntexp_sec = meta_plot[meta_plot['brain_section_label'] == secs[i]]
        fig, ax,  = plot_section_one(ntexp_sec['x'], ntexp_sec['y'], val = gene_column, cmap = 'magma')
        plt.savefig(f'/home/z.starinnov/spatial/slides2/{secs[i]}.png', pad_inches=0.2)
    from PIL import Image
    if name!='no_name':
        images = [
        Image.open(f"/home/z.starinnov/spatial/slides2/{f}.png") for f in secs]
        
        pdf_path = f"/home/z.starinnov/spatial/pdfs/{name}.pdf"
            
        images[0].save(
            pdf_path, "PDF" ,resolution=100.0, save_all=True, append_images=images[1:]
        )
#Used for sorting names of clusters, purely for my own sanity
def get_number(data):
    return int(data.split('_')[0])

#Function for weighing predictive properties of each neigbour of labeled cell
@numba.njit
def weighted_prediction(weights, ref_cats):
    """Get highest weight category."""
    N = len(weights)
    predictions = np.zeros((N,), dtype=ref_cats.dtype)
    uncertainty = np.zeros((N,))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                uncertainty[i] = max(1 - best_prob, 0)

    return predictions, uncertainty