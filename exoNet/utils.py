import anndata
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder


def remove_sparsity(adata: sc.AnnData):
    if issparse(adata.X):
        adata.X = adata.X.A
    return adata


def train_test_split_adata(adata, label_key=None, train_frac=0.85, stratify=True):
    adata = remove_sparsity(adata)

    if label_key and stratify:
        labels = adata.obs[label_key].values
        x_train, x_valid, y_train, y_valid = train_test_split(adata.X, labels, test_size=1. - train_frac,
                                                              stratify=labels)

        train_adata = anndata.AnnData(X=x_train)
        train_adata.obs = y_train
        train_adata.var = adata.var.copy(deep=True)
        train_adata.var_names = adata.var_names

        valid_adata = anndata.AnnData(X=x_valid)
        valid_adata.obs = y_valid
        valid_adata.var = adata.var.copy(deep=True)
        valid_adata.var_names = adata.var_names
    else:
        x_train, x_valid = train_test_split(adata.X, test_size=1. - train_frac)

        train_adata = anndata.AnnData(X=x_train)
        train_adata.var = adata.var.copy(deep=True)
        train_adata.var_names = adata.var_names

        valid_adata = anndata.AnnData(X=x_valid)
        valid_adata.var = adata.var.copy(deep=True)
        valid_adata.var_names = adata.var_names

    return train_adata, valid_adata


def label_encoder(adata, label_encoder, label_key='condition'):
    labels = np.zeros(adata.shape[0])
    if isinstance(label_encoder, dict):
        for condition, label in label_encoder.items():
            labels[adata.obs[label_key] == condition] = label
    elif isinstance(label_encoder, LabelEncoder):
        labels = label_encoder.transform(adata.obs[label_key].values)
    else:
        label_encoder = LabelEncoder()
        labels = label_encoder.fit_transform(adata.obs[label_key].values)
    return labels.reshape(-1, 1), label_encoder
