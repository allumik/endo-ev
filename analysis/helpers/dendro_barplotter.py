import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

def create_row_colormap(pheno_dat, pheno_dat_col, palette="Set2"):
  """
  Creates a colormap for rows based on unique values in pheno_dat_col.

  Args:
      pheno_dat (pd.DataFrame): DataFrame containing phenotype data.
      pheno_dat_col (str): Column name in pheno_dat to use for coloring rows.
      palette (str): Seaborn color palette name.

  Returns:
      dict: Dictionary mapping unique values in pheno_dat_col to hex colors.
  """
  unique_phases = pheno_dat[pheno_dat_col].unique()
  palette = sns.color_palette(palette, unique_phases.size).as_hex()
  return dict(zip(unique_phases, palette))

def assign_row_colors(pheno_dat, pheno_dat_col, frac_pred_index, lut_row):
  """
  Assigns colors to rows based on pheno_dat_col and a lookup table.

  Args:
      pheno_dat (pd.DataFrame): DataFrame containing phenotype data.
      pheno_dat_col (str): Column name used for coloring.
      frac_pred_index (pd.Index): Index of fraction predictions.
      lut_row (dict): Lookup table mapping values to colors.

  Returns:
      pd.Series: Series of colors for each row, indexed by frac_pred_index.
  """
  return (
      pheno_dat
      .reindex(frac_pred_index.values)
      [pheno_dat_col]
      .astype("object")
      .map(lut_row)
  )

def setup_figure_axes(add_dendro):
  """
  Sets up the figure and axes for the dendrogram and barplot.

  Args:
      add_dendro (bool): Whether to include dendrogram in the plot.

  Returns:
      tuple: Figure and a tuple of axes (ax_dendro, ax_bar).
        If add_dendro is False, ax_dendro will be None.
  """
  fig = plt.figure(figsize=(12, 10))
  if add_dendro:
    gs = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)
    ax_dendro = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])
  else:
    gs = GridSpec(1, 1)
    ax_dendro = None
    ax_bar = fig.add_subplot(gs[0])
  return fig, (ax_dendro, ax_bar)

def plot_dendrogram(ax_dendro, fractions, frac_pred_index, pheno_dat, pheno_dat_col):
  """
  Plots the dendrogram on the provided axes.

  Args:
      ax_dendro (plt.Axes): Axes object for the dendrogram.
      fractions (pd.DataFrame): DataFrame of cell fractions.
      frac_pred_index (pd.Index): Index for fraction predictions.
      pheno_dat (pd.DataFrame): Phenotype data.
      pheno_dat_col (str): Column in pheno_dat for labels.

  Returns:
      dict: Dendrogram object returned by scipy.cluster.hierarchy.dendrogram.
  """
  linkage_mat = linkage(fractions.loc[frac_pred_index], method="complete", optimal_ordering=True)
  dend_obj = dendrogram(
      linkage_mat,
      labels=pheno_dat[pheno_dat_col][frac_pred_index].index,
      color_threshold=0,
      above_threshold_color="#000000",
      ax=ax_dendro
  )
  ax_dendro.set_title('TAPE deconvolution results')
  return dend_obj

def style_dendrogram_axes(ax_dendro):
  """
  Applies minimal styling to the dendrogram axes.

  Args:
      ax_dendro (plt.Axes): Dendrogram axes object.
  """
  for spine in ax_dendro.spines.keys():
    ax_dendro.spines[spine].set_visible(False)
  ax_dendro.tick_params(left=False, bottom=False)
  ax_dendro.set_xticks([])
  ax_dendro.set_yticks([])

def plot_barplot(ax_bar, fractions, dend_obj, row_colors):
  """
  Plots the stacked barplot on the provided axes.

  Args:
      ax_bar (plt.Axes): Axes object for the barplot.
      fractions (pd.DataFrame): DataFrame of cell fractions.
      dend_obj (dict): Dendrogram object from plot_dendrogram.
      row_colors (pd.Series): Series of colors for sample names.
  """
  if dend_obj is not None:
    fractions.iloc[dend_obj["leaves"]].plot(kind='bar', stacked=True, figsize=(12, 6), ax=ax_bar)
  else:
    fractions.plot(kind='bar', stacked=True, figsize=(12, 6), ax=ax_bar) # No reordering if dendrogram is absent
  ax_bar.set_xlabel('Sample')
  ax_bar.set_ylabel('Cell Fraction')
  # color the samplenames
  for label in ax_bar.get_xmajorticklabels():
    label.set_color(row_colors.to_dict()[label.get_text()])

def style_barplot_axes(ax_bar):
  """
  Styles the barplot axes (spines, grid, ticks, labels).

  Args:
      ax_bar (plt.Axes): Barplot axes object.
  """
  for spine in ax_bar.spines.keys():
    ax_bar.spines[spine].set_visible(False)
  ax_bar.grid(False)
  ax_bar.set_xticklabels(ax_bar.get_xticklabels(), rotation=90, fontsize=10)

def create_phase_legend(ax_bar, lut_row, pheno_dat_col, legend_phase_bbox, legend_phase_ncol):
  """
  Creates and adds the legend for sample phases to the barplot.

  Args:
      ax_bar (plt.Axes): Barplot axes object.
      lut_row (dict): Lookup table for row colors (phases).
      pheno_dat_col (str): Column name for phase information.
      legend_phase_bbox (tuple): Bounding box for phase legend.
      legend_phase_ncol (int): Number of columns for phase legend.
  """
  legend_elements_phases = [
      Patch(facecolor=color, edgecolor='black', label=phase) for phase, color in lut_row.items()
  ]
  legend_phases = ax_bar.legend(
      handles=legend_elements_phases,
      title=pheno_dat_col,
      bbox_to_anchor=legend_phase_bbox,
      loc='upper left',
      ncol=legend_phase_ncol
  )
  ax_bar.add_artist(legend_phases) # Add to ensure both legends are shown

def create_class_legend(ax_bar, legend_class_bbox, legend_class_ncol):
  """
  Creates and adds the legend for cell classes (bar colors) to the barplot.

  Args:
      ax_bar (plt.Axes): Barplot axes object.
      legend_class_bbox (tuple): Bounding box for class legend.
      legend_class_ncol (int): Number of columns for class legend.
  """
  handles, labels = ax_bar.get_legend_handles_labels()
  ax_bar.legend(
      handles=handles,
      labels=labels,
      title='Cell Types',
      bbox_to_anchor=legend_class_bbox,
      loc='lower left',
      ncol=legend_class_ncol
  )