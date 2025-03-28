import pandas as pd
import numpy as np
import seaborn as sns
import altair as alt
import matplotlib.pyplot as plt
from typing import Tuple
from scipy.stats import spearmanr
from scipy.stats import f_oneway
from scipy.stats import chi2
from scipy.stats import ttest_ind
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.decomposition import PCA

try: from analysis.helpers.dendro_barplotter import *
except ImportError: from helpers.dendro_barplotter import *



# %% Statistical testing functoins
def convert_to_long(fracs, pheno) -> pd.DataFrame:
  # convert to long format
  return (
  fracs
    .reset_index(names="index")
    .melt(id_vars="index", var_name="celltype", value_name="ratio")
    .set_index("index")
    .merge(pheno, left_index=True, right_index=True)
  )

def ttest_table(fracs_long, testing_col="cyclephase", grouping_col="celltype") -> pd.DataFrame:
  ttest_res = {
    grouping_var:
      ttest_ind(*[dat.query(f"`{testing_col}` == @cycle")["ratio"] for cycle in dat[testing_col].unique()])
    for grouping_var, dat in fracs_long.groupby(grouping_col)
  }
  return pd.DataFrame(ttest_res).T.set_axis(["t-stat", "p-val"], axis=1)


# Perform one-way ANOVA and Tukey's test as a function for grouping
def pairwise_tukey_df(dat, testing_col="cyclephase") -> pd.DataFrame:
  tukey_res = pairwise_tukeyhsd(dat["ratio"], dat[testing_col])._results_table
  return pd.DataFrame(
    data=tukey_res.data[1:],
    columns=tukey_res.data[0]
  )

def tukey_table(fracs_long, testing_col="cyclephase", grouping_col="celltype") -> pd.DataFrame:
  tukey_res = {
    celltype:
      pairwise_tukey_df(dat, testing_col=testing_col)
    for celltype, dat in fracs_long.groupby(grouping_col)
  }
  return pd.concat(
    tukey_res, keys=tukey_res.keys(), names=[grouping_col]
    ).reset_index(level=grouping_col)

def anova_table(fracs_long, testing_col="cyclephase", grouping_col="celltype") -> pd.DataFrame:
  anova_res = {
    grouping_var:
      f_oneway(*[dat.query(f"`{testing_col}` == @cycle")["ratio"] for cycle in dat[testing_col].unique()])
    for grouping_var, dat in fracs_long.groupby(grouping_col)
  }
  return pd.DataFrame(anova_res).T.set_axis(["f-stat", "p-val"], axis=1)

def cor_table(fracs_long) -> pd.DataFrame:
  # convert the cyclephase to a dummy var
  fracs_long["dummy_phase"] = pd.Categorical(
    fracs_long.cyclephase,
    categories=["pro", "pre", "rec", "post"],
    ordered=True
    ).codes
  cor_res = {
    celltype:
      spearmanr(dat["dummy_phase"], dat["ratio"])
    for celltype, dat in fracs_long.groupby("celltype")
  }
  return pd.DataFrame(cor_res).T.set_axis(["corr", "p-val"], axis=1)

def mahalanobis_distances(x, refs):
  cov_val = np.cov(refs)
  inv_cov = 1 / cov_val
  diff = x - np.mean(refs)
  mahal = np.sqrt(diff**2 * inv_cov)
  return mahal

# Function to calculate geometric mean, ignoring NaN and non-positive values
def geometric_mean(row):
    valid_values = row[(row > 0) & (abs(row) < np.inf)].dropna()
    if len(valid_values) > 0:
        return np.exp(np.log(valid_values).mean())
    else:
        return np.nan

# Function to calculate statistics for each group
def calculate_stats(group):
  if len(group) == 2:
    ratio_uf = group.iloc[0]
    ratio_biopsy = group.iloc[1]
    correlation, p_value = spearmanr(ratio_uf, ratio_biopsy)
    mean_difference = np.abs(ratio_biopsy.mean() - ratio_uf.mean())
    return pd.Series({'Correlation': correlation, 'P-value': p_value, 'Abs Mean Diff': mean_difference})
  else:
    return None

def get_df_coord(den):
    """
    Get coordinate dataframe.
    """
    # if you view the dendrogram as a collection of upside-down "U" shapes, then
    # we can regard the 4 corners of the upside-down "U" as points 1, 2, 3 and 4.
    cols_xk = ["xk1", "xk2", "xk3", "xk4"]
    cols_yk = ["yk1", "yk2", "yk3", "yk4"]

    df_coord = pd.merge(
        pd.DataFrame(den["icoord"], columns=cols_xk),
        pd.DataFrame(den["dcoord"], columns=cols_yk),
        left_index=True,
        right_index=True
    )
    return df_coord



# %% Functions for plotting
def dendro_barplot(
  fractions: pd.DataFrame,
  phenotype_df: pd.DataFrame,
  legend_phase_bbox: tuple = (1, 0.3),
  legend_class_bbox: tuple = (1, 0.3),
  legend_phase_ncol: int = 1,
  legend_class_ncol: int = 1,
  phenotype_df_col: str = "cyclephase",
  add_dendrogram: bool = True
) -> None:
  """Generates a dendrogram-barplot figure to visualize cell fraction predictions.

  This function creates a combined plot consisting of a dendrogram (optionally)
  and a stacked bar plot. The dendrogram visualizes the hierarchical clustering
  of samples based on their cell fractions. The bar plot shows the predicted
  cell fractions for each sample, with bars stacked according to different
  cell types/states. The rows of the bar plot are colored according to
  phenotype data.

  Args:
      fractions: DataFrame of cell fraction predictions. Rows are samples,
          columns are cell types/states. Values are the predicted fractions.
      phenotype_df: DataFrame of phenotype data. Must have the same index as `fractions`.
      legend_phase_bbox: Bounding box coordinates (x, y) for the phase legend.
          This controls the legend's position.
      legend_class_bbox: Bounding box coordinates (x, y) for the class legend.
      legend_phase_ncol: Number of columns for the phase legend.
      legend_class_ncol: Number of columns for the class legend.
      phenotype_df_col: Column name in `phenotype_df` to use for row
          coloring/labels (e.g., 'cyclephase'). This column determines
          the categories used for coloring the rows of the bar plot.
      add_dendrogram: Whether to include the dendrogram in the plot.
  """

  frac_pred_index = fractions.index

  # 1. Create colormap and assign row colors
  lut_row = create_row_colormap(phenotype_df, phenotype_df_col)
  row_colors = assign_row_colors(phenotype_df, phenotype_df_col, frac_pred_index, lut_row)

  # 2. Setup figure and axes
  fig, (ax_dendro, ax_bar) = setup_figure_axes(add_dendrogram)

  dend_obj = None  # Initialize dend_obj with type hint

  # 3. Plot dendrogram and style axes (conditional)
  if add_dendrogram and ax_dendro is not None:  # Check if ax_dendro is valid
      dend_obj = plot_dendrogram(ax_dendro, fractions, frac_pred_index, phenotype_df, phenotype_df_col)
      style_dendrogram_axes(ax_dendro)

  # 4. Plot barplot and style axes
  plot_barplot(ax_bar, fractions, dend_obj, row_colors)  # Pass dend_obj to barplot
  style_barplot_axes(ax_bar)

  # 5. Create legends
  create_phase_legend(ax_bar, lut_row, phenotype_df_col, legend_phase_bbox, legend_phase_ncol)
  create_class_legend(ax_bar, legend_class_bbox, legend_class_ncol)

  plt.show()



def biplot_fractions(
  fractions_df: pd.DataFrame,
  phenotype_df: pd.DataFrame,
  legend_title: str = "Cycle Phase",
  color_field: str = None, # Added for flexibility
  style_field: str = None, # Added for flexibility
  text_limit = 0.04, # controls on how much labels are shown
  dims = (400, 300)
) -> Tuple[alt.Chart, PCA]:
  """Generates an Altair biplot of Principal Component Analysis (PCA) and its loadings.

  Performs PCA, creates an interactive scatter plot, and overlays loadings. Returns the Altair chart
  and the fitted PCA object.

  Args:
      fractions_df: DataFrame where rows are samples, columns are features (fractions).
      phenotype_df: DataFrame with phenotype data. Must have the same index as `fractions_df`.
      legend_title: Title for the legend.
      color_field: Column in the merged DataFrame to use for color encoding.
      style_field: Column in the merged DataFrame to use for shape encoding.
      text_limit: Controls on how many of loading titles are shown
      dims: The dimensions of the output plot, defaults to (400, 300)

  Returns:
      Tuple[alt.Chart, PCA]: A tuple containing:
          - alt.Chart: The Altair chart object representing the biplot.
          - PCA: The fitted scikit-learn PCA object.
  """

  # PCA calculation
  pca = PCA(n_components=2)
  pca_df = (
    pd.DataFrame(data=pca.fit_transform(fractions_df), columns=["PC1", "PC2"])
    .set_index(fractions_df.index)
    .merge(phenotype_df, left_index=True, right_index=True)
  )
  expl_var = pca.explained_variance_ratio_
  loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
  loadings_df = pd.DataFrame(
    loadings, index=fractions_df.columns, columns=["PC1", "PC2"]
  ).reset_index()  # Reset index for Altair

  # --- Create the base scatter plot ---
  scatter = alt.Chart(pca_df).mark_point(
    size=40, filled=True
  ).encode(
    x=alt.X("PC1:Q", axis=alt.Axis(title=f"PC1 ({expl_var[0]:.2%})", tickCount=1)),
    y=alt.Y("PC2:Q", axis=alt.Axis(title=f"PC2 ({expl_var[1]:.2%})", tickCount=1)),
    color=alt.Color(f"{color_field}:N", title=legend_title) if color_field else "black",
    shape=alt.Shape(f"{style_field}:N") if style_field else "circle"
  )

  # --- Create the loadings plot ---
  loadings_chart = (
    alt.Chart(pd.concat([
      loadings_df, 
      pd.DataFrame({"index": loadings_df["index"], "PC1": 0, "PC2": 0})
      ], ignore_index=True))
    .mark_line(color="red", opacity=0.5)  # Use mark_rule instead of mark_line
    .encode(x='PC1:Q', y='PC2:Q', detail="index")
  )

  # --- Create text labels for loadings (with filtering) ---
  loadings_text = (
    alt.Chart(loadings_df)
      .mark_text(align='center', dx=4, dy=0, color="black", fontSize=10)
      .encode(x="PC1:Q", y="PC2:Q", text='index:N')
      .transform_filter(
        (alt.datum.PC1 > text_limit) |
        (alt.datum.PC1 < -text_limit) |
        (alt.datum.PC2 > text_limit) |
        (alt.datum.PC2 < -text_limit)
      )
  )

  # --- Combine the plots ---
  return (loadings_chart + loadings_text + scatter).properties(
    title="PCA Biplot",
    width=dims[0],  # Adjust as needed
    height=dims[1]   # Adjust as needed
  ), pca



def peruvian_transform(
  fractions: pd.DataFrame,
  phenotype: pd.DataFrame,
  general_cells: pd.DataFrame = None,
  grouping_ids: list = None
  ) -> pd.DataFrame:
  """Transforms and reshapes data for stacked area plots of cell type fractions.

  Prepares a DataFrame for visualization (e.g., with Altair) by:
    1. Optionally grouping and averaging cell type fractions by specified columns.
    2. Reshaping the data from wide to long format (stacking cell types).
    3. Optionally merging cell lineage information based on 'celltype'.

  Args:
      fractions: DataFrame of cell type fraction predictions.  Rows are samples,
          columns are cell types.  Should contain columns to join with `phenotype`.
      phenotype: DataFrame of phenotype data. Must contain the columns specified
          in `grouping_ids`, and index or columns to join with `fractions`.
      general_cells: DataFrame mapping cell types to lineages.
          Requires 'celltype' and 'lineage' columns. If None, no lineage information is added.
      grouping_ids: List of column names in `phenotype` to group by.
          Defaults to ["cyclephase"].  If an empty list is provided (`[]`), no grouping will occur
          and all samples will be present for the plot.

  Returns:
      DataFrame: Transformed DataFrame in long format, ready for plotting.
          Columns include:
          - The columns listed in `grouping_ids`.
          - 'celltype' (the cell type).
          - 'fractions' (the fraction of that cell type).
          - 'lineage' (if `general_cells` is provided).
  """
  # Group and average the samples
  if grouping_ids is not None:
    fractions = fractions.join(phenotype[grouping_ids]).groupby(grouping_ids).mean()

  # Data Reshaping for Altair: Convert to long format
  fractions.columns.names = ["celltype"]
  fractions = (
    fractions
    .stack()
    .reset_index()
    .rename(columns={0: "fractions"})
    .reset_index(drop=True)
  )
  # Add lineage information for the celltype
  if general_cells is not None:  # More concise conditional merge
      return fractions.merge(general_cells, on="celltype")
  else:
      return fractions # Add lineage information for the celltype



def peruvian_sands(
  fractions_df: pd.DataFrame,
  x_axis_col: str = "cyclephase",
  x_order: list = ["pro", "pre", "rec", "post"],
  x_title: str = None,
  global_color_scale: alt.Scale = None,
  lineage_val: str = None,
  grouping_val: str = None,
  is_first_dataset_in_row: bool = True,
  is_first_row: bool = True,
  legend_col_n: int = 1,
  dims: tuple = (200, 150)
) -> alt.Chart:
  """Generates a stacked area chart of cell type fractions across cycle phases.

  Visualizes cell type fractions over cycle phases using Altair-Lite.
  Designed for faceted/concatenated layouts, conditionally displaying
  axis titles, legends, and plot titles based on layout position.

  Args:
      fractions_df: DataFrame with 'cyclephase', 'fractions', and 'celltype' columns.
      cyclephase_order: Order of cycle phases for x-axis.
      lineage_val: Lineage for y-axis title (if first in row).
      grouping_val: Dataset for plot title (if first row).
      is_first_dataset_in_row: If True, display y-axis title and legend.
      is_first_row: If True, display plot title.
      dims: (width, height) of the plot in pixels.

  Returns:
      An Altair-Lite Chart object.
  """
  if global_color_scale is None: global_color_scale = alt.Scale(domain=fractions_df.celltype.unique().tolist())
  return alt.Chart(fractions_df, view=alt.ViewConfig(strokeWidth=0)).mark_area().encode(
    x=alt.X(
      f'{x_axis_col}:O',
      axis=alt.Axis(labelAngle=-45, grid=False, domain=False, title=x_title),
      sort=x_order
      ),
    y=alt.Y(
      'fractions:Q',
      axis=alt.Axis(format='%', grid=False, title=lineage_val) 
        if is_first_dataset_in_row else None
      ),
    color=alt.Color(
      'celltype:N', # Show legend only for the first chart in the row
      scale=global_color_scale,
      legend=alt.Legend(
        title=None,
        values=fractions_df.celltype.unique().tolist(), # Assuming current_row_celltypes is defined in the outer scope or needs to be passed if row-specific
        columns=legend_col_n
        ) if is_first_dataset_in_row else None
      ),
  ).properties( # Title only for first plot, or dataset name as subtitle
    width=dims[0],
    height=dims[1],
    title=f'{grouping_val}' if is_first_row or not None else ""
  )



def peruvian_grouped(
  fractions_df: pd.DataFrame,
  cyclephase_order: list = ["pro", "pre", "rec", "post"],
  grouping_id: str = "dataset",
  global_color_scale: alt.Scale = None,
  dims: tuple = (200, 150)
) -> alt.Chart:
  """Generates a faceted stacked area chart of cell type fractions across cycle phases, faceted by lineage and dataset.

  This function takes cell type fraction predictions, phenotype data, and cell lineage information to create a
  faceted stacked area chart using Altair-Lite. The chart visualizes cell type fractions across different
  cycle phases, with facets arranged by cell lineage (rows) and grouping ID (columns, typically dataset).
  Each row (lineage) has a shared y-axis scale, and the entire grid shares an x-axis scale. A global color
  scheme is used for cell types, and legends are displayed only at the beginning of each row, showing
  cell types present in that lineage.

  Args:
    fractions_df (pd.DataFrame): DataFrame containing cell type fraction predictions.
      Must have columns that will be used to join with `comb_uf_pheno` based on 'sample_id',
      and columns representing cell types (whose names will become the columns of the 'fractions' table
      after grouping and averaging).
    comb_uf_pheno (pd.DataFrame): DataFrame containing phenotype data, including 'cyclephase', and the
      column specified by `grouping_id` (e.g., 'dataset'). Must have a 'sample_id' column to join with `frac_pred`.
    general_cells (pd.DataFrame): DataFrame containing cell type to lineage mapping.
      Must have 'celltype' and 'lineage' columns.
    grouping_id (str, optional): The column name in `comb_uf_pheno` to use for grouping datasets
      horizontally in the faceted plot. Defaults to "dataset".

  Returns:
    Chart: An Altair-Lite Chart object representing the faceted stacked area chart.
  """
  row_charts = [] # List to hold charts for each row (lineage)
  is_first_row = True
  if global_color_scale is None: global_color_scale = alt.Scale(domain=fractions_df.celltype.unique().tolist())

  for lineage_val in fractions_df.lineage.unique():
    dataset_charts_row = [] # List to hold charts for each dataset in the current lineage row
    is_first_dataset_in_row = True # Flag to indicate if it's the first plot in the row (for legend)

    for grouping_val in fractions_df[grouping_id].unique():
      # Create a chart for this dataset and lineage
      current_chart = peruvian_sands(
        fractions_df.query(f"lineage == '{lineage_val}' & {grouping_id} == '{grouping_val}'"),
        x_order=cyclephase_order,
        lineage_val=lineage_val,
        grouping_val=grouping_val,
        is_first_dataset_in_row=is_first_dataset_in_row,
        is_first_row=is_first_row,
        global_color_scale=global_color_scale, # Pass global color scale
        dims=dims
      )
      dataset_charts_row.append(current_chart)
      is_first_dataset_in_row = False
    is_first_row = False

    # Horizontally concatenate charts for the current lineage row
    row_chart = alt.hconcat(*dataset_charts_row, spacing=0).resolve_scale(y='shared')
    row_charts.append(row_chart) # Add the row chart to the list of rows

  # Vertically concatenate all row charts to form the final grid
  return alt.vconcat(*row_charts, spacing=0).resolve_scale( # Adjust vertical spacing as needed
    x='shared',
    y='independent',
    color='independent' # Keep color scale independent - legends are now handled manually
  )


# simple placeholder image generator, use it for composite images where you'll put the image in later
def spacer_with_text(text="[insert plot here]", width=200, height=100, text_size=14, text_color='gray'):
  """Creates a spacer plot with placeholder text centered in the middle.

  Args:
      text: The placeholder text to display.
      width: The width of the spacer plot.
      height: The height of the spacer plot.
      text_size: The size of the placeholder text.
      text_color: The color of the placeholder text.

  Returns:
      An Altair Chart object representing the spacer plot.
  """

  return alt.Chart(
    pd.DataFrame({'x': [0], 'y': [0], 'text': [text]})
  ).mark_text(
    size=text_size,
    color=text_color,
    align='center',
    baseline='middle'
  ).encode(
    x=alt.value(width / 2),  # Center horizontally
    y=alt.value(height / 2), # Center vertically
    text='text'
  ).properties(
    width=width,
    height=height
  )
