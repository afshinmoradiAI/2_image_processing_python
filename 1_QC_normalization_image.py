import argparse
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import warnings
from joblib import parallel_backend
import anndata as ad
import scanpy as sc
import pytometry as pm
import scanpy.external as sce

warnings.filterwarnings('ignore')

def read_and_clean_data(input_path):
    """Reads the input CSV file and cleans the data."""
    df = pd.read_csv(input_path)
    
    # Filter columns with ': Cell: Median' and extract metadata
    data = df.filter(regex=': Cell: Median').copy()
    meta = df[['Image', 'Object ID', 'Name', 'Class', 'Parent', 'ROI', 'Centroid X µm', 'Centroid Y µm',
               'Nucleus: Area µm^2', 'Nucleus: Length µm', 'Nucleus: Circularity', 'Nucleus: Solidity',
               'Nucleus: Max diameter µm', 'Nucleus: Min diameter µm', 'Cell: Area µm^2', 'Cell: Length µm',
               'Cell: Circularity', 'Cell: Solidity', 'Cell: Max diameter µm', 'Cell: Min diameter µm',
               'Nucleus/Cell area ratio']]
    
    # Clean column names
    data.columns = data.columns.str.replace(": Cell: Median", "")
    meta.rename(columns={'Centroid X µm': "x", 'Centroid Y µm': "y"}, inplace=True)
    
    return data, meta

def create_anndata(data, meta):
    """Creates an AnnData object from data and metadata."""
    adata = ad.AnnData(data)
    adata.obs = meta
    adata.obsm["spatial"] = adata.obs[['x', 'y']].to_numpy()
    adata.obs['Image'] = adata.obs['Image'].str.replace(' ', '-')
    return adata

def filter_by_image_counts(adata, min_cells=500):
    """Filters images based on the minimum number of cells."""
    image_counts = adata.obs['Image'].value_counts()
    adata = adata[adata.obs['Image'].isin(image_counts[image_counts > min_cells].index)].copy()
    return adata

def exclude_cores(adata, excluded_cores):
    """Excludes specified cores from the data."""
    if excluded_cores:
        excluded_core_list = excluded_cores.split(',')
        adata = adata[~adata.obs['Image'].isin(excluded_core_list)].copy()
    return adata

def apply_quality_control(adata):
    """Applies a series of quality control filters to the data."""
    # Filter based on DAPI signal
    DAPI_threshold = np.percentile(adata[:, 'DAPI'].X, 10)
    adata = adata[adata[:, 'DAPI'].X > DAPI_threshold, :]

    # Filter cells based on nucleus area
    adata = adata[(adata.obs['Nucleus: Area µm^2'] >= adata.obs['Nucleus: Area µm^2'].quantile(0.01)) &
                  (adata.obs['Nucleus: Area µm^2'] <= adata.obs['Nucleus: Area µm^2'].quantile(0.99))]
    return adata

def exclude_markers(adata, exclude_marker):
    """Excludes specified markers from the data."""
    if exclude_marker:
        exclude_marker_list = exclude_marker.split(',')
        adata = adata[:, ~adata.var_names.isin(exclude_marker_list)]
    return adata

def normalize_and_scale(adata):
    """Normalizes and scales the data."""
    pm.tl.normalize_arcsinh(adata, cofactor=150, inplace=True)
    adata.X = stats.zscore(adata.X, axis=0)
    adata.X = stats.zscore(adata.X, axis=1)
    return adata

def perform_pca_and_harmony(adata):
    """Performs PCA and applies Harmony for batch correction."""
    sc.tl.pca(adata)
    with parallel_backend('threading', n_jobs=32):
        sc.external.pp.harmony_integrate(adata, key='Image', max_iter_harmony=40)
    return adata

def cluster_and_rank_genes(adata, random_seed=49):
    """Clusters the data and ranks genes."""
    sc.settings.seed = random_seed
    sce.tl.phenograph(adata, clustering_algo='leiden', resolution_parameter=0.5, n_jobs=-1, k=20)
    sc.tl.rank_genes_groups(adata, 'pheno_leiden', method='t-test')
    return adata

def save_anndata(adata, output_path):
    """Saves the AnnData object to an H5AD file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # Clean up any large, unnecessary attributes
    if 'pheno_jaccard_ig' in adata.obsp:
        del adata.obsp["pheno_jaccard_ig"]
    adata.write(filename=output_path)

def main(args):
    data, meta = read_and_clean_data(args.input)
    adata = create_anndata(data, meta)
    adata = filter_by_image_counts(adata)
    adata = exclude_cores(adata, args.exclude_core)
    adata = apply_quality_control(adata)
    adata = exclude_markers(adata, args.exclude_marker)
    adata = normalize_and_scale(adata)
    adata = perform_pca_and_harmony(adata)
    adata = cluster_and_rank_genes(adata)
    save_anndata(adata, args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data with QC and normalization.")
    parser.add_argument("--input", required=True, help="Path to the input CSV file")
    parser.add_argument("--output", required=True, help="Path to save the output H5AD file")
    parser.add_argument("--exclude_core", default="", help="Comma-separated list of cores to exclude during QC")
    parser.add_argument("--exclude_marker", default="", help="Comma-separated list of markers to exclude during QC")
    args = parser.parse_args()
    main(args)

