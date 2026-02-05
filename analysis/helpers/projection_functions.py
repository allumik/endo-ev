import os 
import re
import glob
import warnings
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import squidpy as sq
from typing import Tuple
from scipy.stats import ttest_ind, spearmanr, f_oneway, tukey_hsd, pearsonr, wasserstein_distance
from scipy.spatial.distance import jaccard, euclidean
from sklearn.metrics import auc, jaccard_score # More direct for boolean masks

# Define the key where enrichment values reside
OBSM_KEY = "tangram_ct_pred"


def get_spatial_crop_coords(adata: ad.AnnData) -> tuple[float, float, float, float] | None:
  """
  Calculates the bounding box coordinates for spatial data points.

  Args:
    adata: An AnnData object containing spatial coordinates in adata.obsm['spatial'].

  Returns:
    A tuple containing (min_x, min_y, max_x, max_y) representing the
    bounding box, or None if 'spatial' key is not found in adata.obsm.
  """
  if 'spatial' not in adata.obsm:
    print("Error: 'spatial' key not found in adata.obsm")
    return None

  spatial_coords = adata.obsm['spatial']

  # Ensure spatial_coords is a NumPy array for efficient calculation
  if not isinstance(spatial_coords, np.ndarray):
      try:
          spatial_coords = np.array(spatial_coords)
      except Exception as e:
          print(f"Error converting spatial coordinates to NumPy array: {e}")
          return None

  if spatial_coords.shape[1] < 2:
      print("Error: spatial coordinates need at least 2 columns (x and y).")
      return None

  # Assuming column 0 is x and column 1 is y
  min_x = spatial_coords[:, 0].min()
  max_x = spatial_coords[:, 0].max()
  min_y = spatial_coords[:, 1].min()
  max_y = spatial_coords[:, 1].max()

  crop_coords = (min_x, min_y, max_x, max_y)
  return crop_coords

# Function to load AnnData objects from a directory
def load_spatial_data(data_dir, pattern="c2l_*.h5ad"):
  """Loads AnnData objects matching the pattern."""
  data_dict = {}
  file_paths = glob.glob(os.path.join(data_dir, pattern))
  print(file_paths)
  if not file_paths:
    warnings.warn(f"No files found matching pattern '{pattern}' in directory '{data_dir}'")
    return data_dict

  for file_path in file_paths:
    file_name = os.path.basename(file_path)
    match = re.match(r"c2l_(.*?)_(ev|bio|ref)\.h5ad", file_name)
    if match:
      slide_id = match.group(1)
      ref_type = match.group(2)
      try:
        adata = sc.read_h5ad(file_path, backed="r").to_memory()

        data_dict[(slide_id, ref_type)] = adata
      except Exception as e:
        warnings.warn(f"Could not load or process {file_name}: {e}")
    else:
        warnings.warn(f"Filename {file_name} did not match expected pattern.")

  return data_dict


# define the functions for the comparison metrics
# Function to calculate comparison metrics between two AnnData objects
def calculate_comparison_metrics(adata_ref, adata_comp, jaccard_threshold_quantile=0.75):
  results = []
  # adata_comp = vis_dat_ev
  # adata_ref = vis_dat_ref
  cell_types = np.intersect1d(
    adata_comp.obsm[OBSM_KEY].columns,
    adata_ref.obsm[OBSM_KEY].columns
  )

  # --- Get enrichment data from obsm ---
  try:
    X_ref = np.asarray(adata_ref.obsm[OBSM_KEY][cell_types])
    X_comp = np.asarray(adata_comp.obsm[OBSM_KEY][cell_types])
  except KeyError:
    warnings.warn(f"'{OBSM_KEY}' not found in obsm for comparison. Skipping metrics.")
    return pd.DataFrame() # Return empty dataframe

  coords_ref = adata_ref.obsm['spatial']
  coords_comp = adata_comp.obsm['spatial'] # Should be identical

  # --- Per Cell Type Metrics ---
  for i, ct in enumerate(cell_types):
    vec_ref = X_ref[:, i]
    vec_comp = X_comp[:, i]

    # 1. Mean Absolute Difference per Spot
    abs_diff = np.mean(np.abs(vec_ref - vec_comp))
    results.append({'cell_type': ct, 'metric': 'Mean Absolute Difference', 'value': abs_diff})

    # 2. Spot-wise Correlation (Pearson)
    if np.std(vec_ref) > 1e-6 and np.std(vec_comp) > 1e-6:
        corr, _ = pearsonr(vec_ref, vec_comp)
    else:
        corr = 1.0 if np.allclose(vec_ref, vec_comp) else 0.0
    results.append({'cell_type': ct, 'metric': 'Pearson Correlation', 'value': corr})

    # 3. Centroid Distance of Highly Enriched Spots
    q_ref = np.quantile(vec_ref, jaccard_threshold_quantile)
    q_comp = np.quantile(vec_comp, jaccard_threshold_quantile)
    high_mask_ref = vec_ref > q_ref
    high_mask_comp = vec_comp > q_comp

    if np.any(high_mask_ref) and np.any(high_mask_comp):
      centroid_ref = coords_ref[high_mask_ref, :].mean(axis=0)
      centroid_comp = coords_comp[high_mask_comp, :].mean(axis=0)
      centroid_dist = euclidean(centroid_ref, centroid_comp)
    else:
      centroid_dist = np.nan
    results.append({'cell_type': ct, 'metric': f'Centroid Distance (>{jaccard_threshold_quantile*100:.0f}%)', 'value': centroid_dist})

    # 4. Jaccard Index
    if len(high_mask_ref) == len(high_mask_comp):
      jaccard_val = jaccard_score(high_mask_ref, high_mask_comp)
    else:
      jaccard_val = np.nan
    results.append({'cell_type': ct, 'metric': f'Jaccard Index (>{jaccard_threshold_quantile*100:.0f}%)', 'value': jaccard_val})

    # 5. Wasserstein Distance (1D) --- does not add anything in addition to MAD
    # wasserstein_dist = wasserstein_distance(vec_ref, vec_comp)
    # results.append({'cell_type': ct, 'metric': 'Wasserstein Distance (1D)', 'value': wasserstein_dist})

  return pd.DataFrame(results)



# Function to calculate single-adata metrics (Moran's I and Ripley's L AUC)
def calculate_single_adata_metrics(adata):
  results = []
  cell_types = adata.obsm[OBSM_KEY].columns
  metrics_calc = {} # To store intermediate results

  # Check which expected cell types are actually present as columns in obs
  valid_cts_in_obs = [ct for ct in cell_types if ct in adata.obs.columns]
  missing_cts_in_obs = [ct for ct in cell_types if ct not in adata.obs.columns]
  if missing_cts_in_obs:
    warnings.warn(f"Cell types {missing_cts_in_obs} not found as columns in adata.obs. Moran's I will only be calculated for {valid_cts_in_obs}.")
  
  # Store for Ripley, Moran uses attr='obs'
  enrichment_data = adata.obs[valid_cts_in_obs]

  # --- Moran's I (using attr="obs") ---
  # Calculate Moran's I directly using obs columns
  sq.gr.spatial_autocorr(
    adata,                 # Use original adata
    genes=valid_cts_in_obs,# Pass list of column names in obs
    attr="obs",            # Specify reading from obs
    mode="moran",
    n_perms=1000,           # Use consistent params
    n_jobs=1,              # Use consistent params
  )
  # Results are stored in adata.uns['moranI']
  moran_results = adata.uns['moranI']

  # Store results, handling potentially missing calculations
  for ct in cell_types:
    if ct not in metrics_calc: metrics_calc[ct] = {}
    if ct in moran_results.index:
      metrics_calc[ct]['Moran I'] = moran_results.loc[ct, 'I']
      metrics_calc[ct]['Moran I PValue'] = moran_results.loc[ct, 'pval_norm']
      metrics_calc[ct]['Moran I FDR'] = moran_results.loc[ct, 'pval_sim_fdr_bh']
    else:
      # Could be missing if ct wasn't in valid_cts_in_obs or if calculation failed
      metrics_calc[ct]['Moran I'] = np.nan
      metrics_calc[ct]['Moran I PValue'] = np.nan
      metrics_calc[ct]['Moran I FDR'] = np.nan

  # --- Ripley's Statistic (L-function AUC & Min P-value) ---
  # Assign max enrichment cell type
  max_ct_col = '_max_enrichment_ct_temp' # Use temporary column name
  enrichment_matrix = enrichment_data.to_numpy()
  ct_names_for_max = enrichment_data.columns # Use the actual columns used

  if enrichment_matrix is not None and enrichment_matrix.shape[0] > 0 and enrichment_matrix.shape[1] > 0 :
    max_indices = np.argmax(enrichment_matrix, axis=1)
    max_ct_names = ct_names_for_max[max_indices]
    adata.obs[max_ct_col] = pd.Categorical(max_ct_names)
    valid_categories = adata.obs[max_ct_col].value_counts()
    valid_categories = valid_categories[valid_categories > 0].index

    # Calculate Ripley's L-function using copy=True
    max_coord_diff = np.max(adata.obsm['spatial'].max(axis=0) - adata.obsm['spatial'].min(axis=0))
    ripley_results_dict = sq.gr.ripley(
        adata,
        cluster_key=max_ct_col,
        mode='L',
        max_dist=max_coord_diff / 4,
        n_simulations=1000, # Must be > 0 for pvalues
        copy=True # Return results directly
    )

    l_values_df = ripley_results_dict.get('L_stat', pd.DataFrame()).set_index(max_ct_col)
    p_values = ripley_results_dict.get('pvalues', np.array([]))

    # Calculate AUC and Min P-value for each cell type category
    for ct in cell_types:
      if ct not in metrics_calc: metrics_calc[ct] = {}
      if ct in l_values_df.index: # Check if calculated for this ct
        ith_idx = np.where(l_values_df.index.drop_duplicates() == ct)
        # Calculate AUC by subtracting the distances from the l values
        l_minus_r = l_values_df.loc[ct, "stats"] - l_values_df.loc[ct, "bins"]
        ripley_auc = auc(l_values_df.loc[ct, "bins"], l_minus_r)
        # Extract metrics
        metrics_calc[ct]['Ripley L AUC'] = ripley_auc
        min_p_value = np.nanmax(p_values[ith_idx])
        metrics_calc[ct]['Ripley L Max PValue'] = min_p_value

    # --- Convert the collected metrics to list of dicts ---
    for ct, metrics_dict in metrics_calc.items():
      for metric_name, value in metrics_dict.items():
        results.append({
          'cell_type': ct,
          'metric': metric_name,
          'value': value
        })

  return pd.DataFrame(results)