import scanpy as sc


def plot_umap(data, labels=None, show=False):
    adata = sc.AnnData(data)
    adata.obs['labels'] = labels

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='label', show=show, frameon=False)
