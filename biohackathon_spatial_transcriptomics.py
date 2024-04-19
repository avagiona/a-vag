
import scanpy as sc
import squidpy as sq
#import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import umap

seed = 42
#import data

adata = sc.read_10x_h5(filename='cell_feature_matrix.h5')
#adata
df = pd.read_csv('cells.csv')
df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()

#quality control metrics
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# Value 300 from the tutorial percent_top parameter is showing "Positions outside range of features.", I just kept until 200 and it worked
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=(50, 100, 200), inplace=True)

fig, axs = plt.subplots(2, 1, figsize=(15, 6))
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    bins=60,
    ax=axs[1],
)

sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)
#adata

adata.layers["counts"] = adata.X.copy()
#adata
sc.pp.normalize_total(adata, inplace=True)    # Using default parameters, computes CPM normalization
#adata
sc.pp.log1p(adata)                            # Computes log(X+1)
#adata

# Scanpy doesn't support the change of the output data as the native UMAP
# https://umap-learn.readthedocs.io/en/latest/embedding_space.html
hyperbolic_mapper_me = umap.UMAP(output_metric='hyperboloid', random_state=seed).fit(adata.X)
hyperbolic_trans_me = hyperbolic_mapper_me.transform(adata.X)
plt.scatter(hyperbolic_mapper_me.embedding_.T[0],
            hyperbolic_mapper_me.embedding_.T[1], #c=adata.obs.leiden_num,
            cmap='Spectral')
sc.pp.pca(adata)                          # Performs the PCA. With the default parameters, we'll have svd solver arpack, 50 clusters or 1 - minimum dimension size of selected representation
#adata
#adata.varm['PCs'].shape

sc.pp.neighbors(adata)                        # Computes a neighborhood graph of observation
sc.tl.umap(adata)
sc.tl.leiden(adata)



# The cluster numbers must be categorized to be used in the hyperbolic plot
adata.obs["leiden_num"] = adata.obs["leiden"].cat.codes

sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
    size=8, 
)
sc.set_figure_params(figsize=(15,15))
sc.pl.umap(adata, color=["leiden"], save='_classes.png', size=8, palette="Spectral")

adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
sq.pl.spatial_scatter(adata, library_id="spatial", shape=None,color=["leiden"], size=8, palette="Spectral", save="spatial_brain.png", figsize=(15,15))
sq.pl.spatial_scatter(adata, library_id="spatial", shape=None,color=["leiden"], palette="Spectral", save="legend.png", figsize=(15,30))

#NEW CLUSTERS BASED ON HYPER CLUSTERS 
# The cluster numbers must be categorized to be used in the hyperbolic plot
df = pd.read_csv('FINAL.csv')

sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden", 
    ],
    wspace=0.4,
    size=8, 
)
sc.set_figure_params(figsize=(15,15))
sc.pl.umap(adata, color=["leiden"], save='_classes.png', size=8, palette="Spectral")

adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
sq.pl.spatial_scatter(adata, library_id="spatial", shape=None,color=["leiden"], size=8, palette="Spectral", save="spatial_brain.png", figsize=(15,15))
sq.pl.spatial_scatter(adata, library_id="spatial", shape=None,color=["leiden"], palette="Spectral", save="legend.png", figsize=(15,30))















