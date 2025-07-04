import phenograph # pip
import leidenalg # phenograph
import harmonypy
import importlib
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import loompy
import sklearn
import warnings
import os
import subprocess
import struct
from datetime import datetime
import umap
import copy
import scipy
import gc
assert importlib.metadata.version("harmonypy") in ["0.0.5", "0.0.6"]
import typing
from collections import Counter

def adata_obs_to_string(adata):
    tmp_adata = adata.copy()
    for ii in tmp_adata.obs.columns.values:
        tmp_adata.obs[ii] = tmp_adata.obs[ii].astype('str')
    return tmp_adata

def adata_var_to_string(adata):
    tmp_adata = adata.copy()
    for ii in tmp_adata.var.columns.values:
        tmp_adata.var[ii] = tmp_adata.var[ii].astype('str')
    return tmp_adata

def adata_na_format_correct(adata):
    tmp_adata = adata.copy()
    for ii in tmp_adata.obs.columns.values:
        this_metadata = tmp_adata.obs[ii].copy().astype('str')
        this_metadata[this_metadata.str.lower().isin(['na', 'nan'])] = pd.NA
        tmp_adata.obs[ii] = this_metadata
    return tmp_adata

def read_10x_mtx(mtx_path, barcode_path, feature_path, feature_column_name):
    obs_df = pd.read_table(barcode_path, sep="\t", header=None, index_col=0)
    var_df = pd.read_table(feature_path, sep="\t", header=None, names=feature_column_name)
    var_df.index=var_df['gene_ids']
    return(ad.AnnData(X=scipy.sparse.csr_matrix(scipy.io.mmread(mtx_path).astype(np.float32).transpose()), obs=obs_df, var=var_df))

def AssertSameDictAdataVar(dict):
    this_var = None
    for index, this_adata in enumerate(dict.values()):
        if index == 0:
            this_var = this_adata.var.copy()
        else:
            assert np.all(this_var == this_adata.var)
    print("AssertSameDictAdataVar: Pass")
    return True
