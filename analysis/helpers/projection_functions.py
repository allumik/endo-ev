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
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import squidpy as sq
from typing import Tuple
from scipy.stats import ttest_ind, spearmanr, f_oneway, tukey_hsd, pearsonr, wasserstein_distance, combine_pvalues
from scipy.spatial.distance import jaccard, euclidean
from sklearn.metrics import auc, jaccard_score # More direct for boolean masks
from statsmodels.stats.multitest import multipletests

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



#### Functions for the biomarker analysis

def calculate_scores(vis_dat, gene_dict):
  """
  Calculates gene set scores for provided gene groups.
  """
  # Iterate over the dictionary to score each group
  for group, genes in gene_dict.items():
    # Intersect with available genes in raw data
    valid_genes = [g for g in genes if g in vis_dat.raw.var_names]
    
    if len(valid_genes) > 0:
      sc.tl.score_genes(
        vis_dat, 
        gene_list=valid_genes, 
        ctrl_as_ref=False, 
        use_raw=True, 
        score_name=group
      )
  return vis_dat

def calculate_correlations(vis_dat, gene_groups, obsm_key, quantile_threshold=0.75, n_perms=1000):
  """
  Calculates Pearson correlation, Jaccard index, and Bivariate Lee's L 
  between gene group scores and cell type abundances.
  
  Returns matrices for the metrics and the Lee's L P-value (uncorrected).
  """
  if obsm_key not in vis_dat.obsm.keys():
    raise KeyError(f"{obsm_key} not found in adata.obsm")

  # --- 1. Prepare Spatial Weights ---
  if "spatial_connectivities" not in vis_dat.obsp:
    sq.gr.spatial_neighbors(vis_dat)
  
  W = vis_dat.obsp["spatial_connectivities"]
  
  # Row-normalize W (Critical for Lee's L)
  row_sums = np.array(W.sum(axis=1)).flatten()
  with np.errstate(divide='ignore', invalid='ignore'):
    W_norm = W.multiply(1 / row_sums[:, np.newaxis])
    W_norm = W_norm.tocsr() 

  cell_types = vis_dat.obsm[obsm_key].columns
  
  # --- 2. Initialize DataFrames ---
  pearson_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  jaccard_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  lees_df    = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  lees_pval_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)

  # --- 3. Iterate and Calculate ---
  for ct in cell_types:
    # Prepare Cell Type Vector (Y)
    y = vis_dat.obsm[obsm_key][ct].values
    
    # Pre-calculate Spatial Lag of Y (W * y)
    # Standardize Y first (Z_y)
    if np.std(y) > 1e-12:
      y_std = (y - np.mean(y)) / np.std(y)
      y_lag = W_norm.dot(y_std) # Lag of standardized Y
    else:
      y_std, y_lag = None, None

    for group in gene_groups:
      if group not in vis_dat.obs.columns:
        continue
      
      # Prepare Gene Group Vector (X)
      x = vis_dat.obs[group].values
      
      # --- Metric 1: Pearson Correlation (Point-wise) ---
      if np.std(y) > 1e-12 and np.std(x) > 1e-12:
        pearson_df.loc[ct, group] = pearsonr(y, x)[0]
      else:
        pearson_df.loc[ct, group] = 0.0
      
      # --- Metric 2: Jaccard Index ---
      q_y = np.quantile(y, quantile_threshold)
      q_x = np.quantile(x, quantile_threshold)
      jaccard_df.loc[ct, group] = jaccard_score(y > q_y, x > q_x)

      # --- Metric 3: Bivariate Lee's L ---
      # Formula for row-standardized W: L = (W*Z_x) . (W*Z_y) / N
      # This essentially correlates the spatial smoothing of X with spatial smoothing of Y
      if y_lag is not None and np.std(x) > 1e-12:
        x_std = (x - np.mean(x)) / np.std(x)
        x_lag = W_norm.dot(x_std) # Lag of standardized X
        
        N = len(x)
        
        # Observed Lee's L
        # Dot product of the two lag vectors divided by N
        obs_L = (x_lag @ y_lag) / N
        lees_df.loc[ct, group] = obs_L
        
        # Permutation Test
        # We shuffle X, re-calculate its Lag, and correlate with fixed Lag of Y
        sim_Ls = np.zeros(n_perms)
        x_perm = x_std.copy()
        
        for p in range(n_perms):
          np.random.shuffle(x_perm)
          # Crucial: Must re-lag the shuffled vector!
          x_perm_lag = W_norm.dot(x_perm) 
          sim_Ls[p] = (x_perm_lag @ y_lag) / N
        
        # Calculate P-value (Pseudo-significance)
        p_val = (np.sum(sim_Ls >= obs_L) + 1) / (n_perms + 1)
        lees_pval_df.loc[ct, group] = p_val
      else:
        lees_df.loc[ct, group] = np.nan
        lees_pval_df.loc[ct, group] = np.nan

  return pearson_df, jaccard_df, lees_df, lees_pval_df

def calculate_correlations_morans(
  vis_dat, gene_groups, obsm_key, 
  quantile_threshold=0.75, n_perms=1000, 
  apply_fdr_correction=False
  ):
  """
  Calculates Pearson correlation, Jaccard index, and Bivariate Moran's I 
  between gene group scores and cell type abundances.
  
  Returns matrices for the metrics and the Bivariate Moran's I P-value (uncorrected).
  """
  if obsm_key not in vis_dat.obsm.keys():
    raise KeyError(f"{obsm_key} not found in adata.obsm")

  # --- 1. Prepare Spatial Weights ---
  # Calculate spatial neighbors if not present
  if "spatial_connectivities" not in vis_dat.obsp:
    sq.gr.spatial_neighbors(vis_dat)
  
  W = vis_dat.obsp["spatial_connectivities"]
  
  # Row-normalize W for Moran's I calculation
  # This ensures the lag is the average of the neighbors
  row_sums = np.array(W.sum(axis=1)).flatten()
  with np.errstate(divide='ignore', invalid='ignore'):
    W_norm = W.multiply(1 / row_sums[:, np.newaxis])
    W_norm = W_norm.tocsr() 

  cell_types = vis_dat.obsm[obsm_key].columns
  
  # --- 2. Initialize DataFrames ---
  pearson_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  jaccard_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  moran_df   = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)
  moran_pval_df = pd.DataFrame(index=cell_types, columns=gene_groups, dtype=float)

  # --- 3. Iterate and Calculate ---
  for ct in cell_types:
    # Prepare Cell Type Vector (Y)
    y = vis_dat.obsm[obsm_key][ct].values
    
    # Pre-calculate Spatial Lag of Y (W * y_std)
    # Standardizing simplifies the Bivariate Moran's I formula
    if np.std(y) > 1e-12:
      y_std = (y - np.mean(y)) / np.std(y)
      y_lag = W_norm.dot(y_std)
    else:
      y_std, y_lag = None, None

    for group in gene_groups:
      if group not in vis_dat.obs.columns:
        continue
      
      # Prepare Gene Group Vector (X)
      x = vis_dat.obs[group].values
      
      # --- Metric 1: Pearson Correlation ---
      if np.std(y) > 1e-12 and np.std(x) > 1e-12:
        pearson_df.loc[ct, group] = pearsonr(y, x)[0]
      else:
        pearson_df.loc[ct, group] = 0.0
      
      # --- Metric 2: Jaccard Index ---
      q_y = np.quantile(y, quantile_threshold)
      q_x = np.quantile(x, quantile_threshold)
      y_mask = y > q_y
      x_mask = x > q_x
      jaccard_df.loc[ct, group] = jaccard_score(y_mask, x_mask)

      # --- Metric 3: Bivariate Moran's I & P-Value ---
      # Formula: I = (Z_x * W * Z_y) / N
      if y_lag is not None and np.std(x) > 1e-12:
        x_std = (x - np.mean(x)) / np.std(x)
        N = len(x)
        
        # Observed I
        obs_I = (x_std @ y_lag) / N
        moran_df.loc[ct, group] = obs_I
        
        # Permutation Test for P-value
        # Shuffle X (Gene Score) while keeping Y (Cell Type Spatial Lag) fixed
        sim_Is = np.zeros(n_perms)
        x_perm = x_std.copy()
        
        for p in range(n_perms):
          np.random.shuffle(x_perm)
          sim_Is[p] = (x_perm @ y_lag) / N
        
        # Calculate P-value (Pseudo-significance)
        # (Number of simulations >= observed) + 1 / (Total simulations + 1)
        p_val = (np.sum(sim_Is >= obs_I) + 1) / (n_perms + 1)
        
        moran_pval_df.loc[ct, group] = p_val
      else:
        moran_df.loc[ct, group] = np.nan
        moran_pval_df.loc[ct, group] = np.nan

  # --- 4. Apply FDR Correction (Benjamini-Hochberg) ---
  # Flatten, correct, and reshape
  if apply_fdr_correction:
    p_values_flat = moran_pval_df.values.flatten()
    mask = ~np.isnan(p_values_flat)
  
    if np.any(mask):
      _, fdr_flat, _, _ = multipletests(p_values_flat[mask], method='fdr_bh')
      p_values_flat[mask] = fdr_flat
  
    moran_pval_df[:] = p_values_flat.reshape(moran_pval_df.shape)

  return pearson_df, jaccard_df, moran_df, moran_pval_df


def plot_correlation_heatmap(data_df, title, metric_name, pval_df=None, sig_threshold=0.05):
  """
  Plots a hierarchical clustering heatmap with masking for non-significant values.
  """
  # 1. Prepare Data and Mask
  # Fill NaNs with 0 for robust clustering (assumes 0 = no correlation/overlap)
  plot_data = data_df.fillna(0).astype(float)
  
  mask = None
  if pval_df is not None:
    # Ensure p-values align with data dimensions
    pval_aligned = pval_df.reindex(index=plot_data.index, columns=plot_data.columns)
    
    # Create Mask: True where P-value >= threshold OR P-value is missing
    # Masked values will appear white/transparent in the heatmap
    mask = (pval_aligned >= sig_threshold) | pval_aligned.isna()
    
    # Also mask where the original data itself was NaN (if any)
    mask = mask | data_df.isna()

  # Determine colormap center based on metric type
  is_diverging = any(x in metric_name for x in ["Correlation", "Lee's L"])
  cmap = "vlag" if is_diverging else "magma"
  center = 0 if is_diverging else None

  # 2. Plot Clustermap
  g = sns.clustermap(
    plot_data,
    mask=mask,               # Apply the significance mask
    cmap=cmap,
    center=center,
    figsize=(10, 8),
    dendrogram_ratio=(.1, .2),
    cbar_pos=(0.02, 0.32, 0.03, 0.2),
    row_cluster=True,
    col_cluster=False
  )
  
  # 3. Aesthetics
  g.ax_heatmap.set_title(f"{title}\n{metric_name}", pad=20)
  
  # Rotate x-axis labels for readability
  g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
  
  # Add subtitle about significance if masking was applied
  if pval_df is not None:
     g.ax_heatmap.set_xlabel(f"Note: Non-significant values (FDR ≥ {sig_threshold}) are masked.", fontsize=9, labelpad=10)

  plt.show()
def aggregate_correlations(slide_ids, adata_collection, gene_groups, obsm_key, quantile_threshold=0.75):
  """
  Calculates metrics for multiple slides and returns the average scores.
  P-values are combined using Fisher's method.
  """
  # Storage for accumulating results
  accumulators = {
    'pearson': [],
    'jaccard': [],
    'lees_l': [],
    'lees_pval': []
  }

  print(f"Processing {len(slide_ids)} slides: {slide_ids}")

  for slide in slide_ids:
    # 1. Load and Prep Slide
    # Handle the specific tuple key format of your adata_collection
    try:
      vis_dat = adata_collection[(slide, "ref")].copy()
    except KeyError:
      print(f"Skipping {slide}: Not found in collection.")
      continue
      
    # 2. Calculate Scores for this specific slide
    # (Ensure your calculate_scores function handles the raw data correctly)
    vis_dat = calculate_scores(vis_dat, gene_groups)

    # 3. Run Correlations
    p, j, m, m_p = calculate_correlations(
      vis_dat, 
      list(gene_groups.keys()), 
      obsm_key, 
      quantile_threshold
    )
    
    accumulators['pearson'].append(p)
    accumulators['jaccard'].append(j)
    accumulators['lees_l'].append(m)
    accumulators['lees_pval'].append(m_p)

  # 4. Aggregate Results
  # Average the coefficients (Pearson, Jaccard, Lee's L)
  avg_pearson = pd.concat(accumulators['pearson']).groupby(level=0).mean()
  avg_jaccard = pd.concat(accumulators['jaccard']).groupby(level=0).mean()
  avg_lees   = pd.concat(accumulators['lees_l']).groupby(level=0).mean()

  # Combine P-values using Fisher's method
  # We need a custom apply because combine_pvalues expects a list of values
  def fisher_agg(series):
    # Filter NaNs
    valid_p = series.dropna()
    if len(valid_p) == 0: return np.nan
    # Fisher's method returns (statistic, pvalue)
    return combine_pvalues(valid_p, method='fisher')[1]

  # Stack all p-value DFs and apply Fisher's method per cell (cell_type x gene_group)
  # This aligns by index (cell_type) and columns (gene_group)
  combined_pval = pd.concat(accumulators['lees_pval']).groupby(level=0).agg(fisher_agg)

  return avg_pearson, avg_jaccard, avg_lees, combined_pval

def plot_ev_phase_heatmap(cycle_data, gene_col="EV-combined", cell_type_list=None, lineage_info=None, fdr_threshold=0.05):
  """
  Plots a heatmap of EV spatial correlation with robust handling for Categorical data.
  """
  
  # 1. Data Consolidation
  lees_records = {}
  fdr_records = {}
  phase_order = ["Proliferative", "Early-Secretory", "Mid-Secretory"]
  
  for phase in phase_order:
    if phase not in cycle_data: continue
    m_df, f_df = cycle_data[phase]
    
    if gene_col in m_df.columns:
      lees_records[phase] = m_df[gene_col]
      fdr_records[phase] = f_df[gene_col]

  plot_df = pd.DataFrame(lees_records).fillna(0)
  fdr_df = pd.DataFrame(fdr_records)
  
  # 2. Filtering
  if cell_type_list is not None:
    valid_cts = [ct for ct in cell_type_list if ct in plot_df.index]
    plot_df = plot_df.loc[valid_cts]
    fdr_df = fdr_df.loc[valid_cts]
    
    if plot_df.empty:
      print("Warning: No data remaining after filtering.")
      return

  # 3. Lineage Annotation (Row Colors)
  row_colors = None
  lut = None 
  
  if lineage_info is not None:
    # Extract Series from inputs
    if isinstance(lineage_info, dict):
      l_series = pd.Series(lineage_info)
    elif isinstance(lineage_info, pd.DataFrame):
      col = 'lineage' if 'lineage' in lineage_info.columns else lineage_info.columns[0]
      l_series = lineage_info[col]
    elif isinstance(lineage_info, pd.Series):
      l_series = lineage_info
    else:
      raise ValueError("lineage_info must be a Dict, Series, or DataFrame")

    # FIX: Convert to object to avoid Categorical errors when adding "Unknown"
    l_series = l_series.astype(object)

    # Deduplicate index
    l_series = l_series[~l_series.index.duplicated(keep='first')]

    # Align lineage info to the current heatmap rows
    current_lineages = l_series.reindex(plot_df.index).fillna("Unknown")
    unique_lineages = current_lineages.unique()
    
    # Create Color Palette
    palette_name = "tab20" if len(unique_lineages) > 10 else "tab10"
    colors = sns.color_palette(palette_name, len(unique_lineages))
    lut = dict(zip(unique_lineages, colors))
    
    # Map colors (using list comprehension to be safe against index types)
    row_colors = pd.Series(
        [lut.get(x, (0.8, 0.8, 0.8)) for x in current_lineages],
        index=plot_df.index,
        name='Lineage'
    )

  # 4. Annotation Matrix
  annot_df = fdr_df.applymap(lambda x: '*' if x < fdr_threshold else '')
  annot_df = annot_df.fillna('')

  # 5. Plotting
  fig_height = max(4, len(plot_df) * 0.25)
  
  g = sns.clustermap(
    plot_df,
    col_cluster=False, 
    row_cluster=True,
    row_colors=row_colors, 
    annot=annot_df,
    fmt='', 
    cmap="vlag",
    center=0,
    figsize=(4.5, fig_height), 
    dendrogram_ratio=(.15, .05),
    linewidths=0.75,
    linecolor='darkgrey',
    cbar_pos=(1, 0.4, 0.03, 0.3)
  )
  
  # 6. Aesthetics & Legends
  g.ax_heatmap.set_title(f"Global EV Spatial Association\n({gene_col})", pad=20)
  g.ax_heatmap.set_xlabel("")
  g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=-45, ha='left')
  
  # Colorbar
  g.ax_cbar.set_ylabel("Bivariate Lee's L", fontsize=10)
  g.ax_cbar.yaxis.set_label_position("right")
  g.ax_cbar.text(
      0.7, -.15, 
      f"* FDR < {fdr_threshold}", 
      ha='left', va='top', fontsize=9, transform=g.ax_cbar.transAxes
  )

  # Lineage Legend 
  if lut is not None:
    handles = [mpatches.Patch(color=c, label=l) for l, c in lut.items() if l != "Unknown"]
    g.figure.legend(
        handles=handles, 
        title="Lineage", 
        bbox_to_anchor=(.96, .9), 
        loc='upper left', 
        frameon=False
    )

  plt.show()



def plot_nested_dot_plot(phase_results, cell_type_order=None, gene_group_order=None, sig_threshold=0.05):
  """
  Plots a nested dot plot with dynamic color scaling and strict axis ordering.
  
  Parameters:
  -----------
  phase_results : dict
      Dictionary { 'PhaseName': (lees_df, pval_df) }
  cell_type_order : list (Optional)
      List of cell types to sort the Y-axis.
  gene_group_order : list (Optional)
      List of gene groups to sort the X-axis.
  """
  
  # --- 1. Data Preparation ---
  long_data = []
  
  phase_styles = {
      'Proliferative': {'offset': -0.25, 'marker': 'o', 'label': 'Prolif (•)'},
      'Early-Secretory':     {'offset':  0.00, 'marker': 's', 'label': 'Early (■)'},
      'Mid-Secretory':       {'offset':  0.25, 'marker': 'D', 'label': 'Mid (♦)'}
  }
  
  # Consolidate data
  for phase_name, (m_df, p_df) in phase_results.items():
    style = next((v for k, v in phase_styles.items() if k in phase_name), None)
    if style is None: continue
    
    # Process Lee's L
    m_melt = m_df.reset_index().melt(id_vars='index', var_name='gene_group', value_name='lees_l')
    m_melt.rename(columns={'index': 'cell_type'}, inplace=True)
    
    # Process P-values
    p_melt = p_df.reset_index().melt(id_vars='index', var_name='gene_group', value_name='pval')
    p_melt.rename(columns={'index': 'cell_type'}, inplace=True)
    
    merged = pd.merge(m_melt, p_melt, on=['cell_type', 'gene_group'])
    merged['phase'] = phase_name
    merged['x_offset'] = style['offset']
    merged['marker'] = style['marker']
    
    long_data.append(merged)

  df_all = pd.concat(long_data, ignore_index=True)

  # --- 2. Sorting and Ordering (X and Y axes) ---
  
  # A. Cell Type Order (Y-axis)
  if cell_type_order is not None:
    df_all = df_all[df_all['cell_type'].isin(cell_type_order)].copy()
    df_all['cell_type'] = pd.Categorical(
      df_all['cell_type'], categories=cell_type_order, ordered=True
    )
  else:
    # Default to alphabetical if not provided
    df_all = df_all.sort_values('cell_type')

  # B. Gene Group Order (X-axis)
  if gene_group_order is not None:
    df_all = df_all[df_all['gene_group'].isin(gene_group_order)].copy()
    df_all['gene_group'] = pd.Categorical(
      df_all['gene_group'], categories=gene_group_order, ordered=True
    )
  # If no order provided, we leave it as is (pandas concat order) or alphabetical
  
  # Apply the sort
  df_all = df_all.sort_values(['cell_type', 'gene_group'])
  
  if df_all.empty:
      print("Warning: No data found after filtering.")
      return

  # --- 3. Calculate Global Limits ---
  # data_min = df_all['lees_l'].min()
  data_min = 0
  data_max = df_all['lees_l'].max()
  vmin, vmax = data_min, data_max

  # Mappings (Categoricals ensure the order is preserved in .unique())
  unique_cells = df_all['cell_type'].unique()
  unique_genes = df_all['gene_group'].unique()
  
  cell_map = {ct: i for i, ct in enumerate(unique_cells)}
  gene_map = {g: i for i, g in enumerate(unique_genes)}
  
  df_all['y_coord'] = df_all['cell_type'].map(cell_map)
  df_all['x_coord'] = df_all['gene_group'].map(gene_map) + df_all['x_offset']
  
  # Size & Alpha
  df_all['size'] = np.where(df_all['pval'] < sig_threshold, 120, 10)
  df_all['alpha'] = np.where(df_all['pval'] < sig_threshold, 1.0, 0.3)

  # --- 4. Plotting ---
  fig, ax = plt.subplots(figsize=(len(unique_genes)*1.2 + 3, len(unique_cells)*0.5 + 1))
  
  sc = None
  for phase_key, style in phase_styles.items():
    subset = df_all[df_all['phase'].str.contains(phase_key, regex=False)]
    if subset.empty: continue
    
    sc = ax.scatter(
      x=subset['x_coord'],
      y=subset['y_coord'],
      c=subset['lees_l'],
      s=subset['size'],
      marker=style['marker'],
      cmap='vlag',
      vmin=vmin, vmax=vmax,
      edgecolor='k',
      linewidth=0.5,
      alpha=subset['alpha'].values
    )

  # --- 5. Formatting ---
  ax.set_yticks(range(len(unique_cells)))
  ax.set_yticklabels(unique_cells, fontsize=10)
  ax.set_xticks(range(len(unique_genes)))
  ax.set_xticklabels(unique_genes, rotation=45, ha='right', fontsize=10)
  
  # Grid
  ax.set_axisbelow(True)
  for y in np.arange(0.5, len(unique_cells) - 0.5):
    ax.axhline(y, color='lightgrey', linestyle='-', linewidth=0.5)
  for x in np.arange(0.5, len(unique_genes) - 0.5):
    ax.axvline(x, color='lightgrey', linestyle='-', linewidth=0.5)

  ax.invert_yaxis()
  ax.set_xlabel("Gene Groups", fontsize=11, labelpad=10)
  ax.set_title("Spatial Correlation by Cycle Phase", fontsize=13, pad=15)

  # --- 6. Legends ---
  cbar = plt.colorbar(sc, ax=ax, fraction=0.02, pad=0.02)
  cbar.set_label("Bivariate Lee's L", fontsize=10)
  
  legend_handles = []
  legend_handles.append(mlines.Line2D([], [], color='none', label=r'$\bf{Cycle\ Phase}$'))
  
  present_phases = df_all['phase'].unique()
  for name, style in phase_styles.items():
    if any(name in p for p in present_phases):
      legend_handles.append(mlines.Line2D([], [], color='k', marker=style['marker'], 
                          linestyle='None', markersize=8, label=style['label']))
  
  legend_handles.append(mlines.Line2D([], [], color='none', label=' '))
  legend_handles.append(mlines.Line2D([], [], color='none', label=r'$\bf{Significance}$'))
  legend_handles.append(mlines.Line2D([], [], color='grey', marker='o', linestyle='None', 
                      markersize=11, label=f'P < {sig_threshold}'))
  legend_handles.append(mlines.Line2D([], [], color='grey', marker='o', linestyle='None', 
                      markersize=4, alpha=0.5, label='Not Sig.'))

  ax.legend(handles=legend_handles, 
        bbox_to_anchor=(1.20, 1.0), 
        loc='upper left', 
        frameon=False)

  plt.tight_layout()
  plt.show()