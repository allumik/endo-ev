# %% Setup and load some preliminaries

## NB! run this on a gpu node with the endo_ev_preproc environment
# theres a conflict with the torch version and the omicverse version
# so we have endo_ev and endo_ev_preproc environments separately

import pandas as pd
import scanpy as sc
import omicverse as ov
import anndata as an
import numpy as np
from pathlib import Path
from os import getenv
from dotenv import load_dotenv

## set the random seed for torch
RANDO_SEED = 13
import torch
torch.manual_seed(RANDO_SEED)
import random
random.seed(RANDO_SEED)
np.random.seed(RANDO_SEED)


## load the environment variables from the .env file
load_dotenv()
if getenv("DATA_FOLDER") is None:
  load_dotenv(Path.cwd() / ".env")

anndata_folder = getenv("ANNDATA_FOLDER")
atlas_folder = getenv("ATLAS_FOLDER")
data_folder = getenv("DATA_FOLDER")
st_folder = getenv("ST_FOLDER")
model_folder = Path(data_folder).expanduser() / "saved_models"
raw_data_folder = getenv("RAW_DATA_FOLDER")

# some prelims for data standardisation
replace_dict = {
  "Proliferative Early": "EP",
  "Proliferative Late": "LP",
  "Proliferative": "LP",
  "Secretory Early": "ES",
  "Secretory Mid": "MS",
  "Secretory Late": "LS",
  "Secretory Early-Mid": "MS",
  # "Secretory": "Secretory Mid",
  "Mid-Secretory": "MS",
  "Proliferative Disordered": "EP"
}

replace_immu = {
  "Lymphatic": ["uNK1", "uNK1_cycling", "uNK2", "uNK3", "ILC3", "Peripheral_lymphocyte", "Lymphatic", "Immune_Lymphoid"],
  "Myeloid": ["eM1", "eM2", "cDC1", "cDC2", "pDC", "Monocyte", "Immune_Myeloid"],
  "B-cells": ["B_cell", "Plasma_B_cell"],
  "T-cells": ["T_Reg", "T_cell_CD8", "T_cell_CD4", "T_cell_cycling"],
  "Glandular_secretory": ["Glandular_secretory_FGF7"] # also replace the FGF7 cells
}
inv_replace_immu = {value: key for key, values in replace_immu.items() for value in values}
cyclephase_to_include = ["LS", "LP", "EP", "ES", "MS", "Menstrual"]



# %% Do the scRNAseq transformation
snapshot_an_loc = Path(data_folder) / "sc_deconv_snapshot.h5ad"
if not snapshot_an_loc.exists():
  ## Load the main HECA dataset if not already formatted
  sc_dat = sc.read_h5ad(
    Path(raw_data_folder).expanduser() / "endometriumAtlasV2_cells_with_counts.h5ad",
    backed = "r"
  )

  sc_dat_immu = sc.read_h5ad(
    Path(raw_data_folder).expanduser() / "endometriumAtlasV2_cells_immune.h5ad",
    backed = "r"
  )

  sc_dat.obs = (
    sc_dat.obs
    .assign(
      cyclephase=lambda x: x.Stage.replace(replace_dict).astype("category"),
      celltype=lambda x: 
        x.celltype
        .astype("object")
        .mask(
          x.index.isin(sc_dat_immu.obs.celltype.index), 
          sc_dat_immu.obs.celltype.astype("object")
        ),
      label_long=lambda x:
        x.label_long
        .astype("object")
        .mask(
          x.index.isin(sc_dat_immu.obs.label_long.index), 
          sc_dat_immu.obs.label_long.astype("object")
        )
        .astype("category")
    )
    .rename(columns={"sample": "samplename"})
    # now generalise the cell types for better deconvolution classes
    .assign(celltype=lambda x: x.celltype.replace(inv_replace_immu).astype("category"))
  )

  del sc_dat_immu

  ## also exclude the endometriosis samples
  the_condition = "(celltype not in @sc_dat.var.index) & (Endometrial_pathology == 'C') & (cyclephase in @cyclephase_to_include)"
  sc_dat = sc_dat[sc_dat.obs.query(the_condition).index, :].to_memory()

  ## exclude some celltypes, such as she hormones cell types and poorly represented samples
  count_tmp = sc_dat.obs.celltype.value_counts()
  no_go_celltype = list(count_tmp[count_tmp < 100].index) + ['eHormones', 'sHormones', 'dHormones']
  sc_dat = sc_dat[sc_dat.obs.query("celltype not in @no_go_celltype").index, :]

  sc_dat = ov.pp.preprocess(sc_dat.to_memory(), mode="pearson|pearson", n_HVGs=3000)

  ## switch the expression matrix to raw counts if the matrix is preprocessed
  sc_dat.layers["scaled"] = sc_dat.X
  sc_dat.X = sc_dat.layers["counts"]
  del sc_dat.layers["counts"]

  # save the processed dataset
  sc_dat.write(snapshot_an_loc)
else:
  sc_dat = an.read_h5ad(snapshot_an_loc)



# %% load in the UF and biopsy bulk datasets

## load in the EV data
comb_all = pd.read_feather(Path(data_folder).expanduser() / "combined" / "comb_all_batch.feather").set_index("gene_id").iloc[:-1]
comb_all_raw = pd.read_feather(Path(data_folder).expanduser() / "combined" / "comb_all_raw.feather").set_index("external_gene_name").iloc[:-1]
comb_uf = pd.read_feather(Path(data_folder).expanduser() / "combined" / "comb_uf_batch.feather").set_index("gene_id").iloc[:-1]
## filter out some samples
terminator = ["HUT26_UF", "HUT26_biopsy"]
comb_all = comb_all.drop(columns=terminator)
comb_all_raw = comb_all_raw.drop(columns=terminator)
comb_all_pheno = (
  pd.read_table(Path(data_folder).expanduser() / "combined" / "comb_all_pheno.tsv")
  .set_index("samplename")
  .query("samplename not in @terminator")
  .assign(dataset = lambda x: np.where(x.dataset == "HUT", x.dataset, "Vigano"))
)
comb_uf_pheno = comb_all_pheno.query("cyclephase in ['rec', 'pre'] and group == 'UF'")
comb_uf = comb_uf[comb_uf.columns.intersection(comb_uf_pheno.index)]

## CCHT clinical data
clin_de_genes = pd.read_table(Path(data_folder).expanduser() / "clin_de_genes.tsv").locus.to_list()
uf_bio_genes = pd.read_table(Path(data_folder).expanduser() / "combined" / "uf_bio_genes.tsv").locus.to_list()

## load in the SCRATCH clinical bio data
scratch_mat = pd.read_feather(Path(raw_data_folder).expanduser() / "star_mat_pheno.feather").set_index("sample")
cols_to_keep = ["HID", "ages", "Timing_LH", "status_rif", "no_lb"]
scratch_pheno = scratch_mat[cols_to_keep]
scratch_mat = (
  scratch_mat
  .drop(columns=cols_to_keep)
  .transpose()
  .rename_axis("gene_id")
  .infer_objects()
  .iloc[:-1]
)

## load in the EV CCHT only data
ccht_uf_raw = pd.read_feather(Path(data_folder).expanduser() / "filtered" / "annot_raw.feather").set_index("gene_id").iloc[:-1]
ccht_uf_pheno = (
  pd.read_table(Path(data_folder).expanduser() / "filtered" / "phenotype.tsv")
  .set_index("samplename")
  .query("samplename in @ccht_uf_raw.columns")
  .assign(
    cyclephase=lambda x: 
      pd.Categorical(x.cyclephase, categories=["pro", "pre", "rec", "post"], ordered=True)
    )
)
terminator = ["HUT26_UF", "HUT17_UF", "HUT53_UF", "HUT71_UF"]
ccht_uf_pheno = ccht_uf_pheno.query("samplename not in @terminator")
ccht_uf_raw = ccht_uf_raw.drop(columns=terminator)

## load in clinical samples
clin_uf_raw = pd.read_feather(Path(data_folder).expanduser() / "clin_counts_raw.feather").iloc[:-1]
clin_uf_pheno = pd.read_table(Path(data_folder).expanduser() / "phenotype_clin.tsv").set_index("samplename")
clin_uf_pheno = clin_uf_pheno.query("samplename in @clin_uf_raw.columns")
## change the gene_id for the clinical samples
annot = pd.read_csv(Path(raw_data_folder).expanduser() / "annot_table.csv")
clin_uf_raw = (
  clin_uf_raw
  .merge(annot, left_on='gene_id', right_on='ensembl_gene_id')
  .drop(["ensembl_gene_id", "description", "gene_id"], axis=1)
  .groupby("external_gene_name")
  .sum()
)
clin_uf_raw.index.name = "gene_id"
clin_uf_raw = clin_uf_raw.loc[:, clin_uf_pheno.index]



# %% Parameters for the celltype models
model_params = {
  "celltype_key": "celltype",
  "top_marker_num": 250, # reduce it from the default parameter of 500
  # get a half of all the single ~ladies~ cells due to the size of the input dataset
  "max_single_cells": round(len(sc_dat.obs.index) / 8),
  "ratio_num": 1,
  "gpu": 0
}

frac_params = {
  "batch_size": 512,
  "epochs": 1000, # looking at loss plot then 500 seems to be already enough
  "method": "tape",
  "scaler": "ss", # seems to harmonise distributions better
  "mode": "high-resolution"
}

# for the projection part
vae_params = {
  "batch_size": 512,
  "hidden_size": 256,
  "epoch_num": 100 # looking at loss plot then 500 seems to be already enough
}



# %% Run for the celltype models
fractions_folder = Path(data_folder) / "tape_fractions"

## process the batch corrected counts
fractions_file = fractions_folder / "comb_fracs.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=comb_all,
  single_data=sc_dat,
  **model_params
)
frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_pred.to_csv(fractions_file, sep="\t")

## quick check of the model reproducibility with the seed
# frac_pred_old = frac_pred.copy()
# frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
# print(np.subtract(frac_pred, frac_pred_old))

## process the combined raw counts
fractions_file = fractions_folder / "comb_fracs_raw.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=comb_all_raw,
  single_data=sc_dat,
  **model_params
)
frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_pred.to_csv(fractions_file, sep="\t")

# subset only uf samples
fractions_file = fractions_folder / "comb_uf_fracs.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=comb_uf,
  single_data=sc_dat,
  **model_params
)
frac_tmp = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_tmp.to_csv(fractions_file, sep="\t")

## process only the non-batch-corrected ccht UF and biopsy samples
## use this model later for the spatial projection too
fractions_file = fractions_folder / "ccht_fracs.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=ccht_uf_raw,
  single_data=sc_dat,
  **model_params
)
frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_pred.to_csv(fractions_file, sep = "\t")

## the SCRaTCH RIF dataset
fractions_file = fractions_folder / "scratch_fracs.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=scratch_mat,
  single_data=sc_dat,
  **model_params
)
frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_pred.to_csv(fractions_file, sep = "\t")

## the clinical CCHT EV samples
fractions_file = fractions_folder / "clin_fracs.tsv"
model_obj = ov.bulk2single.Bulk2Single(
  bulk_data=clin_uf_raw,
  single_data=sc_dat,
  **model_params
)
frac_pred = model_obj.predicted_fraction(**frac_params, seed=RANDO_SEED)
frac_pred.to_csv(fractions_file, sep = "\t")
