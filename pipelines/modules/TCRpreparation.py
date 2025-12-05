#!/usr/bin/env python3

import os
import pandas as pd
import subprocess
from datetime import datetime

cdate = datetime.now()
DATE_STRING = f"{cdate.year}{cdate.month:02d}{cdate.day:02d}"

class TCRpreparation:
    def __init__(self, input_path:str, run_dir:str):
        self.input_path = input_path
        self.run_dir = run_dir 

        self.preparation()

    def preparation(self):
        
        # Read the dataset
        raw_dataset = pd.read_csv(self.input_path)
        raw_samples_num = raw_dataset.shape[0]

        # Process the dataset to find unique samples
        dataset = raw_dataset.drop_duplicates(subset=['TRA', 'TRB'], inplace=False)
        dataset.reset_index(drop=True, inplace=True)
        dataset['unique_tcr_id'] = [f'TCR{i:04d}' for i in range(len(dataset))]
        self.tcr_dataset = dataset[['unique_tcr_id', 'TRA', 'TRB']]
        unique_samples_num = dataset.shape[0]

        print(f"Number of TCR in raw samples: {raw_samples_num}")
        print(f"Number of unique TCR samples: {unique_samples_num}")

        # assign tcr sample to raw
        raw_dataset = raw_dataset.merge(self.tcr_dataset, on=['TRA', 'TRB'], how='left')
        self.match_ids = raw_dataset[['TCR_ID', 'unique_tcr_id']]

        match_path = os.path.join(self.run_dir, '01_identification', f'match_ids_tcr_{DATE_STRING}.csv')
        self.match_ids.to_csv(match_path, index=False)

        tcr_path = os.path.join(self.run_dir, '01_identification', f'unique_tcr_samples_{DATE_STRING}.csv')
        self.tcr_dataset.to_csv(tcr_path, index=False)

class SlurmTCRModel2:
    def __init__(self, 
                 pipeline_name:str,
                 run_dir:str, 
                 input_model_path:str, 
                 tcrmodel_path:str = '',
                 ignore_pdbs_string:str = '', 
                 relax_structures:str = "True", 
                 max_template_date:str = "2100-01-01", 
                 databases_path:str = "/opt/tcrmodel2/data/databases"):
        
        self.pipeline_name = pipeline_name
        self.run_dir = run_dir
        self.input_model_path = input_model_path
        self.tcrmodel_path = tcrmodel_path
        self.ignore_pdbs_string = ignore_pdbs_string
        self.relax_structures = relax_structures
        self.max_template_date = max_template_date
        self.databases_path = databases_path


        os.makedirs(os.path.join(self.run_dir, '05_logs'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '05_logs', 'TCR'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '02_predictions'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '02_predictions', 'TCR'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '02_predictions', 'TCR', 'TCR_raw'), exist_ok=True)

        self.generate_slurm_template()

        # Save the template to a file
        with open(f'{self.run_dir}/tcrmodel2_{self.pipeline_name}.sh', 'w') as file:
            file.write(self.template)

    def generate_slurm_template(self, dir_slurm:str=None):

        self.template = f'''
#!/bin/bash
#SBATCH --array=1-800%1
#SBATCH --job-name={self.pipeline_name}_tcrmodel2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=short-gpu-big
#SBATCH --gres=gpu:a100:1
#SBATCH --output=05_logs/TCR/slurm.{self.pipeline_name}.%A_%a.out
#SBATCH --error=05_logs/TCR/slurm.{self.pipeline_name}.%A_%a.error

ml load tcrmodel2

echo $SLURM_ARRAY_TASK_ID

value=$((SLURM_ARRAY_TASK_ID))

readarray -t table < <(tail -n +2 01_identification/unique_tcr_samples_{DATE_STRING}.csv | cut -d, -f1-4)

# assuming $value starts at 1
line="${{table[$((value-1))]}}"

sample_id=$(echo "$line" | cut -d, -f1)
TCRA=$(echo "$line" | cut -d, -f2)
TCRB=$(echo "$line" | cut -d, -f3)

echo "Processing sample ID: $sample_id"
echo "TCRA Sequence: $TCRA"
echo "TCRB Sequence: $TCRB"

tcrmodel2_ub_tcr \
--job_id=$sample_id \
--output_dir=/02_predictions/TCR/TCR_raw \
--tcra_seq=$TCRA \
--tcrb_seq=$TCRB \
--ori_db=/database/ \
--tp_db={self.databases_path} \
--relax_structures={self.relax_structures} \
--max_template_date={self.max_template_date}'''

        return self.template
    