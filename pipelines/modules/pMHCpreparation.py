#!/usr/bin/env python3

import os
import pandas as pd
import subprocess
from datetime import datetime

cdate = datetime.now()
DATE_STRING = f"{cdate.year}{cdate.month:02d}{cdate.day:02d}"

class pMHCpreparation:
    def __init__(self, input_path:str,
                 run_dir:str):
        
        self.input_path = input_path
        self.run_dir = run_dir

        print('Create fasta directory')
        self.fastas_dir_path = os.path.join(self.run_dir, '04_assets', 'pMHC_fastas')
        os.makedirs(self.fastas_dir_path, exist_ok=True)

        self.preparation()
        self.create_fastas()

    def preparation(self):
        # Read the dataset
        raw_dataset = pd.read_csv(self.input_path)
        raw_samples_num = raw_dataset.shape[0]

        # Process the dataset to find unique samples
        dataset = raw_dataset.drop_duplicates(subset=['MHCseq', 'peptide'], inplace=False)
        dataset = dataset[~dataset['MHCseq'].isna() & dataset['MHCseq'].notnull() & (dataset['MHCseq'] != '') & (dataset['MHCseq'] != 'nan')]
        dataset.reset_index(drop=True, inplace=True)
        
        dataset['unique_pmhc_id'] = dataset.apply(lambda x: x['assigned_allele'].replace('*','').replace(':', '') + f"_{x['peptide']}", axis=1)

        non_mhc = dataset[dataset['MHCseq'].isna() | (dataset['MHCseq'].isnull()) | (dataset['MHCseq'] == '') | (dataset['MHCseq'] == 'nan')]

        print(f"Number of samples without MHC sequence: {non_mhc.shape[0]}")
        
        self.mhc_dataset = dataset[['unique_pmhc_id', 'MHCseq', 'peptide']]
        
        unique_samples_num = dataset.shape[0]

        print(f"Number of raw samples: {raw_samples_num}")
        print(f"Number of unique samples: {unique_samples_num}")

        # assign tcr sample to raw
        raw_dataset = raw_dataset.merge(self.mhc_dataset, on=['MHCseq', 'peptide'], how='left')
        match_ids = raw_dataset[['TCR_ID', 'unique_pmhc_id']]
        self.match_path = os.path.join(self.run_dir, '01_identification', f'match_ids_pmhc_{DATE_STRING}.csv')
        match_ids.to_csv(self.match_path, index=False)

        self.mhc_path = os.path.join(self.run_dir, '01_identification', f'unique_pmhc_samples_{DATE_STRING}.csv')
        self.mhc_dataset.to_csv(self.mhc_path, index=False)

    def create_fastas(self):
        samples = pd.read_csv(self.mhc_path)

        for index, row in samples.iterrows():
            fasta_path = os.path.join(self.fastas_dir_path, f"{row['unique_pmhc_id']}_{index:03d}.fasta")
            with open(fasta_path, 'w') as fasta_file:
                fasta_file.write(f">{row['unique_pmhc_id']}\n{row['MHCseq']}:{row['peptide']}\n")
            fasta_file.close()

class SlurmColabFold:
    def __init__(self, 
                 pipeline_name:str,
                 run_dir:str, 
                 input_model_path:str, 
                 colabfold_path:str = '/shared/groups/proteindesign/singularity_containers/colabfold_155.sif'):
        
        self.pipeline_name = pipeline_name
        self.run_dir = run_dir
        self.input_model_path = input_model_path
        self.colabfold_path = colabfold_path

        os.makedirs(os.path.join(self.run_dir, '05_logs'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '05_logs', 'pMHC'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '02_predictions'), exist_ok=True)
        os.makedirs(os.path.join(self.run_dir, '02_predictions', 'pMHC'), exist_ok=True)

        self.generate_slurm_template()

        # Save the template to a file
        with open(f'{self.run_dir}/colabfold_{self.pipeline_name}.sh', 'w') as file:
            file.write(self.template)
    
    def generate_slurm_template(self, dir_slurm:str=None):

        self.template = f"""
#!/bin/bash
#SBATCH --array=1-800%1
#SBATCH --job-name={self.pipeline_name}_colabfold
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short-gpu-small
#SBATCH --mem-per-cpu=8G
#SBATCH --gres=gpu:1g.5gb:1
#SBATCH --output=05_logs/pMHC/slurm.{self.pipeline_name}.%A_%a.out
#SBATCH --error=05_logs/pMHC/slurm.{self.pipeline_name}.%A_%a.error

echo "Job ID Array: $SLURM_ARRAY_JOB_ID"

COLABFOLD_SIF={self.colabfold_path}
FASTA_DIR=04_assets/pMHC_fastas

# essa variável aponta para o arquivo fasta (MUDE PARA O SEU ARQUIVO)
FASTA_FILE=$(find "$FASTA_DIR" -type f -name "*_$SLURM_ARRAY_TASK_ID.fasta" | head -n 1)
basename=$(basename "$FASTA_FILE" .fasta)

# nome da pasta onde os modelos e resultados serão salvos (PODE MUDAR PARA UM NOME QUE ESCOLHER)
OUTPUT_DIR=/02_predictions/pMHC/$basename

echo "Processing $FASTA_FILE"
echo "Basename: $basename"
echo "Output directory: $OUTPUT_DIR"
echo "Using ColabFold SIF: $COLABFOLD_SIF"

# run alphafold
singularity run --nv $COLABFOLD_SIF colabfold_batch \
    $FASTA_FILE \
    $OUTPUT_DIR \
    --templates \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 20 \
    --amber"""
    
        return self.template
