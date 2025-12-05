#!/usr/bin/env python3

from modules.TCRpreparation import TCRpreparation, SlurmTCRModel2 
from modules.pMHCpreparation import pMHCpreparation, SlurmColabFold
#from modules.renumber_pMHC import renumber_pmhc
from modules.gather_results import find_tcr_models, find_pmhc_models, gather_pmhc_models
import os
import pandas as pd
from datetime import datetime


# Get current date for file naming
cdate = datetime.now()
DATE_STRING = f"{cdate.year}{cdate.month:02d}{cdate.day:02d}"

# Define paths and create directories
pipeline_name = "vdjdb_score2_wojust10x"
raw_data_path = f"../databases/VDJdb/02_processed/fullSeqs_dataset_score2_wojust10x_{DATE_STRING}.csv"

# Create Run environment directories
runs_dir     = f"../runs/{pipeline_name}"
id_path      = f"{runs_dir}/01_identification"
pred_path    = f"{runs_dir}/02_predictions"
results_path = f"{runs_dir}/03_results"
assets_path  = f"{runs_dir}/04_assets"
logs_path    = f"{runs_dir}/05_logs"

os.makedirs(runs_dir, exist_ok=True)
os.makedirs(id_path, exist_ok=True)
os.makedirs(pred_path, exist_ok=True)
os.makedirs(results_path, exist_ok=True)
os.makedirs(assets_path, exist_ok=True)
os.makedirs(logs_path, exist_ok=True)

# Remove just10x references from VDJdb score2 dataset
vdjdb_score2_path = "../databases/VDJdb/02_processed/fullSeqs_dataset_score2_20251021.csv"
just10x_ref_path = "assets/immrep_just10x_exclusion.csv"

raw_df = pd.read_csv(vdjdb_score2_path)
reference = pd.read_csv(just10x_ref_path)
just10x_df = reference[reference['Exclude'] == 1]

just10x_df['PMID'] = just10x_df['PMID'].replace({'PMID: ': ''})

filtered_df = raw_df[~raw_df['Reference'].isin(just10x_df['PMID'])]
filtered_df = filtered_df[~filtered_df['Reference'].isin(just10x_df['VDJdb_alias'])]
filtered_df.reset_index(drop=True, inplace=True)
filtered_df.to_csv(raw_data_path, index=False)

# statistics
print(f"Original dataset size: {raw_df.shape[0]}")
print(f"Filtered dataset size: {filtered_df.shape[0]}")

# Run TCRmodel2
TCRpreparation(
    input_path=raw_data_path,
    run_dir=runs_dir)

SlurmTCRModel2(
    pipeline_name=pipeline_name,
    run_dir=runs_dir,
    input_model_path=os.path.join(id_path, f'unique_tcr_samples_{DATE_STRING}.csv'))

# Run ColabFold for pMHCs
pMHCpreparation(
    input_path=raw_data_path,
    run_dir=runs_dir)

SlurmColabFold(
    pipeline_name=pipeline_name,
    run_dir=runs_dir,
    input_model_path=os.path.join(pred_path, f'unique_pmhc_samples_{DATE_STRING}.csv'),
    colabfold_path='/shared/groups/proteindesign/singularity_containers/colabfold_155.sif')
